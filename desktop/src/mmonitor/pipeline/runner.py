#!/usr/bin/env python3
"""
Unified Snakemake runner for MMonitor pipelines.

Handles running Snakemake workflows natively (Linux/macOS) or via Docker (Windows).
Used by both the CLI and GUI.
"""
import os
import sys
import platform
import subprocess
import tempfile
import logging
import shutil
from pathlib import Path
from typing import Optional, Callable

import yaml

logger = logging.getLogger(__name__)


def _find_workflows_dir() -> Path:
    """Locate the workflows directory.

    Search order:
    1. Inside the installed package: mmonitor/workflows/Snakefile
    2. Walk up the filesystem (dev-checkout fallback)
    """
    # 1. Package-internal (pip install / editable install)
    pkg_workflows = Path(__file__).resolve().parent.parent / 'workflows'
    if (pkg_workflows / 'Snakefile').exists():
        return pkg_workflows

    # 2. Walk up from this file (legacy dev tree: desktop/workflows/)
    current = Path(__file__).resolve().parent
    while current != current.parent:
        for candidate in [current / 'workflows', current / 'desktop' / 'workflows']:
            if (candidate / 'Snakefile').exists():
                return candidate
        current = current.parent

    raise FileNotFoundError(
        "Cannot locate workflows/Snakefile. "
        "Ensure you are running from the MMonitor project tree or that "
        "the mmonitor package is properly installed."
    )


DOCKER_IMAGE = "mmonitor/pipelines:latest"


class SnakemakeRunner:
    """
    Unified runner that invokes Snakemake targets for all MMonitor pipelines.

    Supports:
    - Native execution with conda environments (Linux/macOS)
    - Docker-based execution (Windows or explicit --docker flag)
    - Progress callbacks for GUI integration
    """

    def __init__(self, workflow_dir: Optional[Path] = None,
                 use_docker: Optional[bool] = None):
        """
        Args:
            workflow_dir: Path to the workflows/ directory containing Snakefile.
                          Auto-detected if not provided.
            use_docker: Force Docker mode. If None, auto-detect (True on Windows).
        """
        if workflow_dir is not None:
            self.workflow_dir = Path(workflow_dir)
        else:
            self.workflow_dir = _find_workflows_dir()

        if use_docker is not None:
            self.use_docker = use_docker
        else:
            self.use_docker = self._should_use_docker()

    @staticmethod
    def _should_use_docker() -> bool:
        """Auto-detect whether to use Docker (always on Windows)."""
        return platform.system() == 'Windows'

    @staticmethod
    def docker_available() -> bool:
        """Check if Docker is installed and accessible."""
        try:
            result = subprocess.run(
                ['docker', '--version'],
                capture_output=True, text=True, timeout=10
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    @staticmethod
    def snakemake_available() -> bool:
        """Check if Snakemake is installed and accessible."""
        try:
            result = subprocess.run(
                ['snakemake', '--version'],
                capture_output=True, text=True, timeout=10
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    def run(self, target: str, config: dict, threads: int = 4,
            dry_run: bool = False, verbose: bool = False,
            callback: Optional[Callable[[str], None]] = None) -> int:
        """
        Run a Snakemake target.

        Args:
            target: Snakemake target rule name (e.g. 'taxonomy_16s', 'assembly')
            config: Configuration dictionary matching the Snakemake config schema
            threads: Number of cores to use
            dry_run: Only show what would be done
            verbose: Enable verbose Snakemake output
            callback: Optional function called with each line of stdout/stderr
                      for progress reporting (GUI integration)

        Returns:
            Exit code (0 for success)
        """
        if self.use_docker:
            return self._run_docker(target, config, threads, dry_run, verbose, callback)
        return self._run_native(target, config, threads, dry_run, verbose, callback)

    def _run_native(self, target: str, config: dict, threads: int,
                    dry_run: bool, verbose: bool,
                    callback: Optional[Callable[[str], None]]) -> int:
        """Run Snakemake natively with conda environments."""
        if not self.snakemake_available():
            msg = (
                "Snakemake is not installed or not in PATH. "
                "Install with: pip install snakemake"
            )
            logger.error(msg)
            if callback:
                callback(f"ERROR: {msg}")
            return 1

        config_path = self._write_config(config)
        try:
            cmd = [
                'snakemake',
                '--snakefile', str(self.workflow_dir / 'Snakefile'),
                '--configfile', config_path,
                '--cores', str(threads),
                target,
                '--use-conda',
                '--conda-frontend', 'mamba',
            ]

            if dry_run:
                cmd.append('--dry-run')
            if verbose:
                cmd.append('--verbose')

            logger.info("Running: %s", ' '.join(cmd))

            return self._execute(cmd, callback)
        finally:
            self._cleanup_config(config_path)

    def _run_docker(self, target: str, config: dict, threads: int,
                    dry_run: bool, verbose: bool,
                    callback: Optional[Callable[[str], None]]) -> int:
        """Run Snakemake inside a Docker container."""
        if not self.docker_available():
            msg = (
                "Docker is not installed or not running. "
                "Install Docker Desktop from https://www.docker.com/products/docker-desktop/"
            )
            logger.error(msg)
            if callback:
                callback(f"ERROR: {msg}")
            return 1

        # Resolve host paths for volume mounts
        output_dir = Path(config.get('output_dir', 'results')).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write config with container-internal paths
        docker_config = config.copy()
        docker_config['output_dir'] = '/data/output'

        # Map input files to container paths
        input_fastq = config.get('input_fastq', '')
        volume_mounts = []
        if isinstance(input_fastq, list) and input_fastq:
            input_dir = str(Path(input_fastq[0]).resolve().parent)
            volume_mounts.append(f"{input_dir}:/data/input:ro")
            docker_config['input_fastq'] = [
                f"/data/input/{Path(f).name}" for f in input_fastq
            ]
        elif isinstance(input_fastq, str) and input_fastq:
            input_dir = str(Path(input_fastq).resolve().parent)
            volume_mounts.append(f"{input_dir}:/data/input:ro")
            docker_config['input_fastq'] = f"/data/input/{Path(input_fastq).name}"

        volume_mounts.append(f"{output_dir}:/data/output")

        # Mount database directories
        databases = config.get('databases', {})
        db_mounts = {}
        for db_name, db_path in databases.items():
            if db_path and Path(db_path).exists():
                container_path = f"/data/databases/{db_name}"
                volume_mounts.append(f"{Path(db_path).resolve()}:{container_path}:ro")
                db_mounts[db_name] = container_path
        docker_config['databases'] = db_mounts

        # Also handle tool-specific db paths (e.g. config['emu']['database'])
        for tool in ('emu', 'centrifuger', 'checkm2', 'gtdbtk', 'bakta'):
            tool_cfg = docker_config.get(tool, {})
            if isinstance(tool_cfg, dict) and tool_cfg.get('database'):
                host_path = Path(tool_cfg['database']).resolve()
                if host_path.exists():
                    container_path = f"/data/databases/{tool}"
                    volume_mounts.append(f"{host_path}:{container_path}:ro")
                    tool_cfg['database'] = container_path

        config_path = self._write_config(docker_config)
        try:
            cmd = ['docker', 'run', '--rm']
            for mount in volume_mounts:
                cmd.extend(['-v', mount])
            cmd.extend(['-v', f"{config_path}:/opt/mmonitor/config.yaml:ro"])

            cmd.extend([
                DOCKER_IMAGE,
                '--configfile', '/opt/mmonitor/config.yaml',
                '--use-conda',
                '--cores', str(threads),
            ])

            if dry_run:
                cmd.append('--dry-run')
            if verbose:
                cmd.append('--verbose')

            cmd.append(target)

            logger.info("Running Docker: %s", ' '.join(cmd))

            return self._execute(cmd, callback)
        finally:
            self._cleanup_config(config_path)

    def _execute(self, cmd: list,
                 callback: Optional[Callable[[str], None]]) -> int:
        """Execute a command, optionally streaming output to a callback."""
        if callback:
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1
            )
            for line in proc.stdout:
                line = line.rstrip('\n')
                callback(line)
                logger.debug(line)
            proc.wait()
            return proc.returncode
        else:
            result = subprocess.run(cmd, check=False)
            return result.returncode

    @staticmethod
    def _write_config(config: dict) -> str:
        """Write config dict to a temporary YAML file. Returns the path."""
        fd, path = tempfile.mkstemp(suffix='.yaml', prefix='mmonitor_')
        with os.fdopen(fd, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
        return path

    @staticmethod
    def _cleanup_config(path: str):
        """Remove temporary config file."""
        try:
            os.unlink(path)
        except OSError:
            pass
