#!/usr/bin/env python3
"""
MMonitor Database Manager

Handles downloading, building, and managing databases for all analysis tools.
Implements FAIR principles with versioning, checksums, and manifests.

Supported databases:
- EMU (16S rRNA taxonomy)
- Centrifuger (WGS taxonomy)
- GTDB-TK (MAG taxonomy)
- CheckM2 (MAG quality)
- Bakta (gene annotation)
"""
import os
import sys
import json
import hashlib
import shutil
import tarfile
import tempfile
import subprocess
import logging
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List, Any
from dataclasses import dataclass, asdict
from urllib.parse import urlparse

import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)

# ============ Conda Environment Mapping ============
# Maps tools to their conda environment names
TOOL_CONDA_ENVS = {
    'bakta': 'mmonitor-smk-annotation',
    'bakta_db': 'mmonitor-smk-annotation',
    'centrifuger': 'mmonitor-smk-centrifuger',
    'centrifuger-download': 'mmonitor-smk-centrifuger',
    'centrifuger-build': 'mmonitor-smk-centrifuger',
    'checkm2': 'mmonitor-smk-checkm2',
    'gtdbtk': 'mmonitor-smk-gtdbtk',
    'minimap2': 'mmonitor-smk-taxonomy',
    'emu': 'mmonitor-smk-taxonomy',
}

# Path to conda environment yaml files (relative to workflows dir)
CONDA_ENV_YAMLS = {
    'mmonitor-smk-annotation': 'envs/annotation.yaml',
    'mmonitor-smk-taxonomy': 'envs/taxonomy.yaml',
    'mmonitor-smk-centrifuger': 'envs/centrifuger.yaml',
    'mmonitor-smk-checkm2': 'envs/checkm2.yaml',
    'mmonitor-smk-gtdbtk': 'envs/gtdbtk.yaml',
}

# ============ Database Presets ============

EMU_PRESETS = {
    'default': {
        'name': 'EMU Default',
        'description': 'Default EMU database (hosted on OSF)',
        'osf_project': '56uf7',
        'osf_path': 'osfstorage/emu-prebuilt/emu.tar',
        'version': '1.0',
        'size_mb': 100,
    },
    'silva138': {
        'name': 'SILVA 138.2',
        'description': 'SILVA 138.2 Small Subunit rRNA database (hosted on OSF)',
        'osf_project': '56uf7',
        'osf_path': 'osfstorage/emu-prebuilt/silva-138.2.tar',
        'version': '138.2',
        'size_mb': 500,
    },
    'rdp': {
        'name': 'RDP',
        'description': 'Ribosomal Database Project database (hosted on OSF)',
        'osf_project': '56uf7',
        'osf_path': 'osfstorage/emu-prebuilt/rdp.tar',
        'version': 'latest',
        'size_mb': 200,
    },
}

# Prebuilt indices from Dropbox (ready to use, no build step)
# Check https://github.com/mourisl/centrifuger for newer versions
CENTRIFUGER_PREBUILT = {
    'gtdb': {
        'name': 'GTDB r226 (prebuilt)',
        'description': 'GTDB r226 prebuilt index (2025-05-20). Check https://github.com/mourisl/centrifuger for updates.',
        'url': 'https://www.dropbox.com/scl/fo/xjp5r81jxkzxest9ijxul/ADfYFKoxIyl0hrICeEI63QM?rlkey=5lij0ocrbre165pa52mavux5z&dl=1',
        'filename': 'gtdb_r226.zip',
        'version': 'r226',
        'size_gb': 176,
    },
    'nt': {
        'name': 'NCBI core nt (prebuilt)',
        'description': 'NCBI core nt 2025 prebuilt index (2025-05-20). Check https://github.com/mourisl/centrifuger for updates.',
        'url': 'https://www.dropbox.com/scl/fo/f1mbf7nf893pisoruanb4/AHS06LaJr9EN0Pg7hbifWn8?rlkey=7fgtj6pi53l2iwrjw1k6xq8o8&dl=1',
        'filename': 'core_nt_2025.zip',
        'version': '2025',
        'size_gb': 212,
    },
}

# Build-from-sequences presets (download + build, slower)
CENTRIFUGER_PRESETS = {
    'bacteria': {
        'name': 'NCBI RefSeq Bacteria',
        'description': 'RefSeq complete bacterial genomes (download + build)',
        'domain': 'bacteria',
        'version': 'latest',
        'size_gb': 50,
    },
    'archaea': {
        'name': 'NCBI RefSeq Archaea',
        'description': 'RefSeq complete archaeal genomes (download + build)',
        'domain': 'archaea',
        'version': 'latest',
        'size_gb': 5,
    },
    'viral': {
        'name': 'NCBI RefSeq Viral',
        'description': 'RefSeq complete viral genomes (download + build)',
        'domain': 'viral',
        'version': 'latest',
        'size_gb': 2,
    },
    'fungi': {
        'name': 'NCBI RefSeq Fungi',
        'description': 'RefSeq fungal genomes (download + build)',
        'domain': 'fungi',
        'version': 'latest',
        'size_gb': 10,
    },
}

GTDBTK_RELEASES = {
    'r220': {
        'name': 'GTDB R220',
        'url': 'https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz',
        'version': 'R220',
        'size_gb': 85,
    },
    'r214': {
        'name': 'GTDB R214',
        'url': 'https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz',
        'version': 'R214',
        'size_gb': 80,
    },
}

CHECKM2_RELEASES = {
    'latest': {
        'name': 'CheckM2 DIAMOND DB',
        'url': 'https://zenodo.org/record/5571251/files/checkm2_database.tar.gz',
        'version': '1.0',
        'size_gb': 3,
    }
}

BAKTA_RELEASES = {
    'full': {
        'name': 'Bakta Full Database',
        'description': 'Full Bakta annotation database',
        'version': '5.1',
        'size_gb': 60,
        # Bakta has its own download command
    },
    'light': {
        'name': 'Bakta Light Database',
        'description': 'Light Bakta annotation database (smaller)',
        'version': '5.1',
        'size_gb': 3,
    }
}


# ============ Data Classes ============

@dataclass
class DatabaseManifest:
    """FAIR-compliant database manifest."""
    name: str
    tool: str
    version: str
    created: str
    source: str
    description: str
    files: List[Dict[str, Any]]
    total_size: int
    checksum: str  # SHA256 of all files combined
    doi: Optional[str] = None
    build_command: Optional[str] = None

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> 'DatabaseManifest':
        data = json.loads(json_str)
        return cls(**data)

    def save(self, path: Path):
        with open(path, 'w') as f:
            f.write(self.to_json())

    @classmethod
    def load(cls, path: Path) -> 'DatabaseManifest':
        with open(path) as f:
            return cls.from_json(f.read())


@dataclass
class DatabaseInfo:
    """Information about an installed database."""
    name: str
    tool: str
    path: str
    version: str
    installed: str
    size_bytes: int
    verified: bool = False


# ============ Database Manager ============

class DatabaseManager:
    """
    Unified database management for all MMonitor tools.

    Handles downloading, building, verifying, and tracking databases.
    """

    def __init__(self, base_dir: Optional[Path] = None):
        """
        Initialize database manager.

        Args:
            base_dir: Base directory for databases. Defaults to ~/mmonitor_databases
        """
        self.base_dir = Path(base_dir) if base_dir else Path.home() / 'mmonitor_databases'
        self.registry_file = self.base_dir / 'registry.json'
        self._registry = None

    @property
    def registry(self) -> Dict[str, List[DatabaseInfo]]:
        """Load or return cached registry."""
        if self._registry is None:
            self._registry = self._load_registry()
        return self._registry

    def _load_registry(self) -> Dict[str, List[DatabaseInfo]]:
        """Load database registry from file."""
        if self.registry_file.exists():
            with open(self.registry_file) as f:
                data = json.load(f)
                # Convert to DatabaseInfo objects
                registry = {}
                for tool, dbs in data.items():
                    registry[tool] = [DatabaseInfo(**db) for db in dbs]
                return registry
        return {}

    def _save_registry(self):
        """Save registry to file."""
        self.base_dir.mkdir(parents=True, exist_ok=True)
        data = {}
        for tool, dbs in self.registry.items():
            data[tool] = [asdict(db) for db in dbs]
        with open(self.registry_file, 'w') as f:
            json.dump(data, f, indent=2)

    def _register_database(self, info: DatabaseInfo):
        """Register a database in the registry."""
        if info.tool not in self.registry:
            self.registry[info.tool] = []

        # Remove existing entry with same name
        self.registry[info.tool] = [
            db for db in self.registry[info.tool]
            if db.name != info.name
        ]
        self.registry[info.tool].append(info)
        self._save_registry()

    def list_databases(self) -> Dict[str, List[DatabaseInfo]]:
        """List all installed databases."""
        return self.registry

    def get_database(self, tool: str, name: str = None) -> Optional[DatabaseInfo]:
        """Get a specific database."""
        if tool not in self.registry:
            return None

        dbs = self.registry[tool]
        if not dbs:
            return None

        if name:
            for db in dbs:
                if db.name == name:
                    return db
            return None

        # Return first/default
        return dbs[0]

    # ============ Download Utilities ============

    def _download_file(self, url: str, dest: Path, desc: str = None) -> bool:
        """
        Download a file with progress bar.

        Args:
            url: URL to download
            dest: Destination path
            desc: Description for progress bar

        Returns:
            True if successful
        """
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()

            total_size = int(response.headers.get('content-length', 0))

            dest.parent.mkdir(parents=True, exist_ok=True)

            with open(dest, 'wb') as f:
                with tqdm(total=total_size, unit='B', unit_scale=True, desc=desc or dest.name) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                        pbar.update(len(chunk))

            return True

        except Exception as e:
            logger.error(f"Download failed: {e}")
            return False

    def _extract_archive(self, archive: Path, dest: Path) -> bool:
        """Extract tar/tar.gz/tar.bz2 archive (auto-detects format)."""
        try:
            dest.mkdir(parents=True, exist_ok=True)

            with tarfile.open(archive, 'r:*') as tar:
                tar.extractall(dest)

            return True

        except Exception as e:
            logger.error(f"Extraction failed: {e}")
            return False

    def _calculate_checksum(self, path: Path) -> str:
        """Calculate SHA256 checksum of file or directory."""
        sha256 = hashlib.sha256()

        if path.is_file():
            with open(path, 'rb') as f:
                for chunk in iter(lambda: f.read(8192), b''):
                    sha256.update(chunk)
        else:
            # For directories, hash all files sorted by name
            for file_path in sorted(path.rglob('*')):
                if file_path.is_file():
                    sha256.update(str(file_path.relative_to(path)).encode())
                    with open(file_path, 'rb') as f:
                        for chunk in iter(lambda: f.read(8192), b''):
                            sha256.update(chunk)

        return sha256.hexdigest()

    def _get_dir_size(self, path: Path) -> int:
        """Get total size of directory in bytes."""
        total = 0
        for file_path in path.rglob('*'):
            if file_path.is_file():
                total += file_path.stat().st_size
        return total

    def _get_workflows_dir(self) -> Path:
        """Get the workflows directory path."""
        # Try to find workflows dir relative to this file
        this_file = Path(__file__).resolve()
        # Go up from cli/database_manager.py to desktop/
        desktop_dir = this_file.parent.parent.parent.parent
        workflows_dir = desktop_dir / 'workflows'
        if workflows_dir.exists():
            return workflows_dir
        # Fallback: try relative to cwd
        cwd_workflows = Path.cwd() / 'workflows'
        if cwd_workflows.exists():
            return cwd_workflows
        return workflows_dir

    def _ensure_conda_env(self, env_name: str) -> bool:
        """
        Ensure a conda environment exists, creating it if necessary.

        Args:
            env_name: Name of the conda environment

        Returns:
            True if environment exists or was created successfully
        """
        # Check if environment already exists
        try:
            result = subprocess.run(
                ['conda', 'env', 'list'],
                capture_output=True, text=True, check=True
            )
            if env_name in result.stdout:
                logger.debug(f"Conda environment '{env_name}' already exists")
                return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.warning(f"Could not check conda environments: {e}")
            # Try mamba as fallback
            try:
                result = subprocess.run(
                    ['mamba', 'env', 'list'],
                    capture_output=True, text=True, check=True
                )
                if env_name in result.stdout:
                    return True
            except (subprocess.CalledProcessError, FileNotFoundError):
                pass

        # Environment doesn't exist, try to create it
        if env_name in CONDA_ENV_YAMLS:
            yaml_path = self._get_workflows_dir() / CONDA_ENV_YAMLS[env_name]
            if yaml_path.exists():
                logger.info(f"Creating conda environment '{env_name}' from {yaml_path}...")
                try:
                    # Try mamba first (faster), fallback to conda
                    for pkg_manager in ['mamba', 'conda']:
                        try:
                            subprocess.run(
                                [pkg_manager, 'env', 'create', '-f', str(yaml_path)],
                                check=True
                            )
                            logger.info(f"Successfully created environment '{env_name}'")
                            return True
                        except FileNotFoundError:
                            continue
                        except subprocess.CalledProcessError as e:
                            logger.error(f"Failed to create environment with {pkg_manager}: {e}")
                            continue
                except Exception as e:
                    logger.error(f"Failed to create conda environment: {e}")
                    return False
            else:
                logger.warning(f"Environment yaml not found: {yaml_path}")

        return False

    def _run_in_conda_env(self, cmd: List[str], env_name: str,
                          check: bool = True, **kwargs) -> subprocess.CompletedProcess:
        """
        Run a command inside a conda environment.

        Args:
            cmd: Command and arguments as a list
            env_name: Name of the conda environment
            check: Whether to raise exception on non-zero exit
            **kwargs: Additional arguments passed to subprocess.run

        Returns:
            CompletedProcess instance
        """
        # Ensure the environment exists
        if not self._ensure_conda_env(env_name):
            logger.warning(f"Conda environment '{env_name}' may not exist, attempting to run anyway...")

        # Build the command to run inside conda environment
        # Using conda run which properly activates the environment
        full_cmd = ['conda', 'run', '-n', env_name, '--no-capture-output'] + cmd

        logger.debug(f"Running in conda env '{env_name}': {' '.join(cmd)}")

        try:
            return subprocess.run(full_cmd, check=check, **kwargs)
        except FileNotFoundError:
            # Try mamba as fallback
            full_cmd = ['mamba', 'run', '-n', env_name, '--no-capture-output'] + cmd
            return subprocess.run(full_cmd, check=check, **kwargs)

    def _get_conda_env_for_tool(self, tool: str) -> Optional[str]:
        """Get the conda environment name for a tool."""
        return TOOL_CONDA_ENVS.get(tool)

    def _create_manifest(self, db_path: Path, tool: str, name: str,
                         version: str, source: str, description: str,
                         build_command: str = None, doi: str = None) -> DatabaseManifest:
        """Create a FAIR-compliant manifest for a database."""
        files = []
        for file_path in sorted(db_path.rglob('*')):
            if file_path.is_file():
                files.append({
                    'name': str(file_path.relative_to(db_path)),
                    'size': file_path.stat().st_size,
                    'sha256': self._calculate_checksum(file_path)
                })

        total_size = sum(f['size'] for f in files)
        overall_checksum = self._calculate_checksum(db_path)

        manifest = DatabaseManifest(
            name=name,
            tool=tool,
            version=version,
            created=datetime.now().isoformat(),
            source=source,
            description=description,
            files=files,
            total_size=total_size,
            checksum=overall_checksum,
            doi=doi,
            build_command=build_command
        )

        # Save manifest
        manifest.save(db_path / 'manifest.json')

        return manifest

    # ============ EMU Database ============

    def download_emu_database(self, preset: str, output_dir: Path = None,
                              force: bool = False) -> Optional[Path]:
        """
        Download EMU database via osfclient.

        Args:
            preset: Database preset (default, silva138, rdp)
            output_dir: Output directory
            force: Overwrite existing

        Returns:
            Path to database directory
        """
        if preset not in EMU_PRESETS:
            logger.error(f"Unknown EMU preset: {preset}")
            logger.info(f"Available presets: {', '.join(EMU_PRESETS.keys())}")
            return None

        info = EMU_PRESETS[preset]

        if output_dir is None:
            output_dir = self.base_dir / 'emu' / preset
        else:
            output_dir = Path(output_dir)

        if output_dir.exists() and not force:
            # Check if it actually has database files
            if (output_dir / 'species_taxid.fasta').exists():
                logger.warning(f"Database already exists at {output_dir}")
                logger.info("Use --force to overwrite")
                return output_dir

        logger.info(f"Downloading EMU {info['name']} database...")
        logger.info(f"Size: ~{info['size_mb']} MB")

        # Ensure osfclient is available
        try:
            subprocess.run(['osf', '--help'], capture_output=True, check=True)
        except (FileNotFoundError, subprocess.CalledProcessError):
            logger.info("Installing osfclient...")
            subprocess.run([sys.executable, '-m', 'pip', 'install', 'osfclient'],
                           check=True)

        output_dir.mkdir(parents=True, exist_ok=True)

        # Download via osfclient
        tar_file = output_dir / 'emu.tar'
        try:
            logger.info(f"Fetching from OSF project {info['osf_project']}...")
            result = subprocess.run(
                ['osf', '-p', info['osf_project'], 'fetch',
                 info['osf_path'], str(tar_file)],
                capture_output=True, text=True, cwd=str(output_dir)
            )
            if result.returncode != 0:
                logger.error(f"osf fetch failed: {result.stderr}")
                return None

            # Extract
            logger.info("Extracting database...")
            if not self._extract_archive(tar_file, output_dir):
                return None

        finally:
            # Clean up tar file
            if tar_file.exists():
                tar_file.unlink()

        # Create manifest
        manifest = self._create_manifest(
            db_path=output_dir,
            tool='emu',
            name=f"emu-{preset}",
            version=info['version'],
            source=f"osf.io/{info['osf_project']}",
            description=info['description']
        )

        # Register
        db_info = DatabaseInfo(
            name=f"emu-{preset}",
            tool='emu',
            path=str(output_dir),
            version=info['version'],
            installed=datetime.now().isoformat(),
            size_bytes=manifest.total_size,
            verified=True
        )
        self._register_database(db_info)

        logger.info(f"EMU database installed at: {output_dir}")
        return output_dir

    def build_emu_database(self, fasta: Path, taxonomy: Path,
                           output_dir: Path, name: str,
                           threads: int = 4) -> Optional[Path]:
        """
        Build custom EMU database.

        Args:
            fasta: Input FASTA file with sequences
            taxonomy: TSV file with taxonomy mapping
            output_dir: Output directory
            name: Database name
            threads: Number of threads for indexing

        Returns:
            Path to database directory
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Building EMU database '{name}'...")

        # Copy taxonomy file
        shutil.copy(taxonomy, output_dir / 'taxonomy.tsv')

        # Copy and rename FASTA
        shutil.copy(fasta, output_dir / 'species_taxid.fasta')

        # Build minimap2 index - run in conda environment
        logger.info("Building minimap2 index...")
        env_name = self._get_conda_env_for_tool('minimap2')
        try:
            cmd = [
                'minimap2', '-d',
                str(output_dir / 'species_taxid.fasta.mmi'),
                str(output_dir / 'species_taxid.fasta'),
                '-t', str(threads)
            ]
            if env_name:
                self._run_in_conda_env(cmd, env_name)
            else:
                subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to build minimap2 index: {e}")
            return None
        except FileNotFoundError:
            logger.error("minimap2 not found. Please ensure the conda environment is set up.")
            logger.info("Try: mamba env create -f workflows/envs/taxonomy.yaml")
            return None

        # Create manifest
        build_cmd = f"mmonitor database emu build -f {fasta} -t {taxonomy} -o {output_dir} -n {name}"
        manifest = self._create_manifest(
            db_path=output_dir,
            tool='emu',
            name=name,
            version='custom',
            source=str(fasta),
            description=f'Custom EMU database built from {fasta.name}',
            build_command=build_cmd
        )

        # Register
        db_info = DatabaseInfo(
            name=name,
            tool='emu',
            path=str(output_dir),
            version='custom',
            installed=datetime.now().isoformat(),
            size_bytes=manifest.total_size,
            verified=True
        )
        self._register_database(db_info)

        logger.info(f"Custom EMU database created at: {output_dir}")
        return output_dir

    def restrict_emu_database(self, source_db: Path, taxids: List[int],
                              output_dir: Path, name: str) -> Optional[Path]:
        """
        Create restricted EMU database for specific taxa.

        Args:
            source_db: Source EMU database directory
            taxids: List of taxonomy IDs to include
            output_dir: Output directory
            name: Database name

        Returns:
            Path to restricted database
        """
        from Bio import SeqIO
        import pandas as pd

        source_db = Path(source_db)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Creating restricted EMU database with {len(taxids)} taxa...")

        # Load taxonomy
        tax_file = source_db / 'taxonomy.tsv'
        if not tax_file.exists():
            logger.error(f"Taxonomy file not found: {tax_file}")
            return None

        tax_df = pd.read_csv(tax_file, sep='\t')

        # Filter taxonomy
        taxids_set = set(str(t) for t in taxids)
        filtered_tax = tax_df[tax_df['tax_id'].astype(str).isin(taxids_set)]

        if len(filtered_tax) == 0:
            logger.error("No matching taxa found in database")
            return None

        logger.info(f"Found {len(filtered_tax)} matching taxa")

        # Save filtered taxonomy
        filtered_tax.to_csv(output_dir / 'taxonomy.tsv', sep='\t', index=False)

        # Filter FASTA
        fasta_file = source_db / 'species_taxid.fasta'
        if not fasta_file.exists():
            logger.error(f"FASTA file not found: {fasta_file}")
            return None

        kept_seqs = []
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Extract tax_id from header (format: >seqid|taxid|...)
            parts = record.id.split('|')
            if len(parts) >= 2 and parts[1] in taxids_set:
                kept_seqs.append(record)

        logger.info(f"Kept {len(kept_seqs)} sequences")

        # Write filtered FASTA
        SeqIO.write(kept_seqs, output_dir / 'species_taxid.fasta', 'fasta')

        # Build index - run in conda environment
        logger.info("Building minimap2 index...")
        env_name = self._get_conda_env_for_tool('minimap2')
        try:
            cmd = [
                'minimap2', '-d',
                str(output_dir / 'species_taxid.fasta.mmi'),
                str(output_dir / 'species_taxid.fasta')
            ]
            if env_name:
                self._run_in_conda_env(cmd, env_name)
            else:
                subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to build index: {e}")
            return None

        # Create manifest
        manifest = self._create_manifest(
            db_path=output_dir,
            tool='emu',
            name=name,
            version='restricted',
            source=str(source_db),
            description=f'Restricted EMU database with {len(taxids)} taxa'
        )

        # Register
        db_info = DatabaseInfo(
            name=name,
            tool='emu',
            path=str(output_dir),
            version='restricted',
            installed=datetime.now().isoformat(),
            size_bytes=manifest.total_size,
            verified=True
        )
        self._register_database(db_info)

        logger.info(f"Restricted database created at: {output_dir}")
        return output_dir

    # ============ Centrifuger Database ============

    def download_centrifuger_database(self, preset: str, output_dir: Path = None,
                                       threads: int = 4, force: bool = False) -> Optional[Path]:
        """
        Download and build Centrifuger database.

        Uses centrifuger-download to fetch genomes from NCBI.

        Args:
            preset: Database preset (bacteria, archaea, viral, fungi)
            output_dir: Output directory
            threads: Download threads
            force: Overwrite existing

        Returns:
            Path to database
        """
        if preset not in CENTRIFUGER_PRESETS:
            logger.error(f"Unknown Centrifuger preset: {preset}")
            logger.info(f"Available presets: {', '.join(CENTRIFUGER_PRESETS.keys())}")
            return None

        info = CENTRIFUGER_PRESETS[preset]

        if output_dir is None:
            output_dir = self.base_dir / 'centrifuger' / preset
        else:
            output_dir = Path(output_dir)

        if output_dir.exists() and not force:
            logger.warning(f"Database already exists at {output_dir}")
            logger.info("Use --force to overwrite")
            return output_dir

        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Downloading Centrifuger {info['name']} database...")
        logger.info(f"This will download ~{info['size_gb']} GB of data")
        logger.warning("This process may take several hours depending on your connection")

        # Get conda environment for centrifuger
        env_name = self._get_conda_env_for_tool('centrifuger-download')

        # Step 1: Download taxonomy
        logger.info("Step 1/3: Downloading NCBI taxonomy...")
        try:
            cmd = ['centrifuger-download', '-o', str(output_dir / 'taxonomy'), '-d', 'taxonomy']
            if env_name:
                self._run_in_conda_env(cmd, env_name)
            else:
                subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to download taxonomy: {e}")
            return None
        except FileNotFoundError:
            logger.error("centrifuger-download not found. Please ensure the conda environment is set up.")
            logger.info("Try: mamba env create -f workflows/envs/centrifuger.yaml")
            return None

        # Step 2: Download genomes
        logger.info(f"Step 2/3: Downloading {info['domain']} genomes...")
        try:
            cmd = [
                'centrifuger-download', '-o', str(output_dir / 'library'),
                '-d', info['domain'],
                '-t', str(threads)
            ]
            if env_name:
                self._run_in_conda_env(cmd, env_name)
            else:
                subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to download genomes: {e}")
            return None

        # Step 3: Build index
        logger.info("Step 3/3: Building Centrifuger index...")
        db_prefix = output_dir / f"centrifuger_{preset}"
        try:
            cmd = [
                'centrifuger-build',
                '-p', str(threads),
                '--taxonomy-tree', str(output_dir / 'taxonomy' / 'nodes.dmp'),
                '--name-table', str(output_dir / 'taxonomy' / 'names.dmp'),
                '-r', str(output_dir / 'library'),
                str(db_prefix)
            ]
            if env_name:
                self._run_in_conda_env(cmd, env_name)
            else:
                subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to build index: {e}")
            return None

        # Create manifest
        manifest = self._create_manifest(
            db_path=output_dir,
            tool='centrifuger',
            name=f"centrifuger-{preset}",
            version=datetime.now().strftime('%Y%m%d'),
            source=f"NCBI RefSeq {info['domain']}",
            description=info['description'],
            build_command=f"mmonitor database centrifuger download --preset {preset}"
        )

        # Register
        db_info = DatabaseInfo(
            name=f"centrifuger-{preset}",
            tool='centrifuger',
            path=str(db_prefix),
            version=datetime.now().strftime('%Y%m%d'),
            installed=datetime.now().isoformat(),
            size_bytes=manifest.total_size,
            verified=True
        )
        self._register_database(db_info)

        logger.info(f"Centrifuger database built at: {db_prefix}")
        return db_prefix

    def download_centrifuger_prebuilt(self, preset: str, output_dir: Path = None,
                                       force: bool = False) -> Optional[Path]:
        """
        Download prebuilt Centrifuger index from Dropbox.

        Much faster than download+build since the index is already compiled.

        Args:
            preset: Prebuilt preset (gtdb, nt)
            output_dir: Output directory
            force: Overwrite existing

        Returns:
            Path to database directory
        """
        if preset not in CENTRIFUGER_PREBUILT:
            logger.error(f"Unknown prebuilt preset: {preset}")
            logger.info(f"Available prebuilt: {', '.join(CENTRIFUGER_PREBUILT.keys())}")
            return None

        info = CENTRIFUGER_PREBUILT[preset]

        if output_dir is None:
            output_dir = self.base_dir / 'centrifuger' / preset
        else:
            output_dir = Path(output_dir)

        if output_dir.exists() and not force:
            # Check for .cfr files
            cfr_files = list(output_dir.glob('*.cfr'))
            if cfr_files:
                logger.warning(f"Database already exists at {output_dir} ({len(cfr_files)} .cfr files)")
                logger.info("Use --force to overwrite")
                return output_dir

        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Downloading prebuilt {info['name']}...")
        logger.info(f"Size: ~{info['size_gb']} GB")
        logger.info("Note: check https://github.com/mourisl/centrifuger for newer prebuilt indices")

        zip_file = output_dir / info['filename']

        # Download with wget (supports resume with -c)
        try:
            logger.info(f"Downloading {info['filename']}...")
            cmd = [
                'wget', '-c', '-O', str(zip_file),
                '--progress=dot:giga',
                info['url']
            ]
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Download failed: {e}")
            return None
        except FileNotFoundError:
            # Fallback to curl if wget not available
            logger.info("wget not found, trying curl...")
            try:
                cmd = ['curl', '-L', '-C', '-', '-o', str(zip_file), info['url']]
                subprocess.run(cmd, check=True)
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                logger.error(f"Download failed: {e}")
                return None

        # Extract
        logger.info("Extracting prebuilt index...")
        try:
            import zipfile
            with zipfile.ZipFile(zip_file, 'r') as zf:
                zf.extractall(output_dir)
            zip_file.unlink()
        except Exception as e:
            logger.error(f"Extraction failed: {e}")
            return None

        # Find the database prefix (look for .cfr files)
        cfr_files = list(output_dir.rglob('*.1.cfr'))
        if cfr_files:
            # Database prefix is the path without the .1.cfr suffix
            db_prefix = str(cfr_files[0]).replace('.1.cfr', '')
            logger.info(f"Database prefix: {db_prefix}")
        else:
            db_prefix = str(output_dir)
            logger.warning("Could not find .cfr index files. Check the extracted contents.")

        # Create manifest
        manifest = self._create_manifest(
            db_path=output_dir,
            tool='centrifuger',
            name=f"centrifuger-{preset}-prebuilt",
            version=info['version'],
            source=info['url'],
            description=info['description']
        )

        # Register
        db_info = DatabaseInfo(
            name=f"centrifuger-{preset}-prebuilt",
            tool='centrifuger',
            path=db_prefix,
            version=info['version'],
            installed=datetime.now().isoformat(),
            size_bytes=manifest.total_size,
            verified=True
        )
        self._register_database(db_info)

        logger.info(f"Prebuilt Centrifuger index installed at: {db_prefix}")
        return Path(db_prefix)

    def build_centrifuger_database(self, fasta_dir: Path, taxonomy_dir: Path,
                                    output_prefix: Path, name: str,
                                    threads: int = 4) -> Optional[Path]:
        """
        Build custom Centrifuger database.

        Args:
            fasta_dir: Directory with genome FASTA files
            taxonomy_dir: Directory with NCBI taxonomy files (nodes.dmp, names.dmp)
            output_prefix: Output database prefix
            name: Database name
            threads: Number of threads

        Returns:
            Path to database prefix
        """
        logger.info(f"Building custom Centrifuger database '{name}'...")

        env_name = self._get_conda_env_for_tool('centrifuger-build')
        try:
            cmd = [
                'centrifuger-build',
                '-p', str(threads),
                '--taxonomy-tree', str(taxonomy_dir / 'nodes.dmp'),
                '--name-table', str(taxonomy_dir / 'names.dmp'),
                '-r', str(fasta_dir),
                str(output_prefix)
            ]
            if env_name:
                self._run_in_conda_env(cmd, env_name)
            else:
                subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to build index: {e}")
            return None
        except FileNotFoundError:
            logger.error("centrifuger-build not found. Please ensure the conda environment is set up.")
            logger.info("Try: mamba env create -f workflows/envs/centrifuger.yaml")
            return None

        # Register
        db_info = DatabaseInfo(
            name=name,
            tool='centrifuger',
            path=str(output_prefix),
            version='custom',
            installed=datetime.now().isoformat(),
            size_bytes=self._get_dir_size(output_prefix.parent),
            verified=True
        )
        self._register_database(db_info)

        logger.info(f"Custom Centrifuger database built at: {output_prefix}")
        return output_prefix

    # ============ GTDB-TK Database ============

    def download_gtdbtk_database(self, release: str = 'r220',
                                  output_dir: Path = None,
                                  force: bool = False) -> Optional[Path]:
        """
        Download GTDB-TK database.

        Args:
            release: GTDB release (r220, r214)
            output_dir: Output directory
            force: Overwrite existing

        Returns:
            Path to database directory
        """
        if release not in GTDBTK_RELEASES:
            logger.error(f"Unknown GTDB-TK release: {release}")
            logger.info(f"Available releases: {', '.join(GTDBTK_RELEASES.keys())}")
            return None

        info = GTDBTK_RELEASES[release]

        if output_dir is None:
            output_dir = self.base_dir / 'gtdbtk' / release
        else:
            output_dir = Path(output_dir)

        if output_dir.exists() and not force:
            logger.warning(f"Database already exists at {output_dir}")
            return output_dir

        logger.info(f"Downloading GTDB-TK {info['name']} database...")
        logger.info(f"Size: ~{info['size_gb']} GB")
        logger.warning("This is a large download and may take several hours")

        with tempfile.TemporaryDirectory() as tmpdir:
            archive = Path(tmpdir) / 'gtdbtk.tar.gz'

            if not self._download_file(info['url'], archive, f"GTDB-TK {release}"):
                return None

            logger.info("Extracting database (this may take a while)...")
            output_dir.mkdir(parents=True, exist_ok=True)

            if not self._extract_archive(archive, output_dir):
                return None

        # Register
        db_info = DatabaseInfo(
            name=f"gtdbtk-{release}",
            tool='gtdbtk',
            path=str(output_dir),
            version=info['version'],
            installed=datetime.now().isoformat(),
            size_bytes=self._get_dir_size(output_dir),
            verified=True
        )
        self._register_database(db_info)

        logger.info(f"GTDB-TK database installed at: {output_dir}")
        logger.info(f"Set environment variable: export GTDBTK_DATA_PATH={output_dir}")

        return output_dir

    # ============ CheckM2 Database ============

    def download_checkm2_database(self, output_dir: Path = None,
                                   force: bool = False) -> Optional[Path]:
        """
        Download CheckM2 database.

        Args:
            output_dir: Output directory
            force: Overwrite existing

        Returns:
            Path to database directory
        """
        info = CHECKM2_RELEASES['latest']

        if output_dir is None:
            output_dir = self.base_dir / 'checkm2'
        else:
            output_dir = Path(output_dir)

        if output_dir.exists() and not force:
            logger.warning(f"Database already exists at {output_dir}")
            return output_dir

        logger.info(f"Downloading CheckM2 database...")
        logger.info(f"Size: ~{info['size_gb']} GB")

        # CheckM2 has its own download command - run in conda environment
        env_name = self._get_conda_env_for_tool('checkm2')
        try:
            cmd = ['checkm2', 'database', '--download', '--path', str(output_dir)]
            if env_name:
                self._run_in_conda_env(cmd, env_name)
            else:
                subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            # Fallback to direct download
            logger.info("Falling back to direct download...")
            with tempfile.TemporaryDirectory() as tmpdir:
                archive = Path(tmpdir) / 'checkm2_db.tar.gz'

                if not self._download_file(info['url'], archive, "CheckM2"):
                    return None

                output_dir.mkdir(parents=True, exist_ok=True)
                if not self._extract_archive(archive, output_dir):
                    return None
        except FileNotFoundError:
            logger.error("checkm2 not found. Please ensure the conda environment is set up.")
            logger.info("Try: mamba env create -f workflows/envs/checkm2.yaml")
            return None

        # Register
        db_info = DatabaseInfo(
            name='checkm2',
            tool='checkm2',
            path=str(output_dir),
            version=info['version'],
            installed=datetime.now().isoformat(),
            size_bytes=self._get_dir_size(output_dir),
            verified=True
        )
        self._register_database(db_info)

        logger.info(f"CheckM2 database installed at: {output_dir}")
        return output_dir

    # ============ Bakta Database ============

    def download_bakta_database(self, db_type: str = 'light',
                                 output_dir: Path = None,
                                 force: bool = False) -> Optional[Path]:
        """
        Download Bakta database.

        Args:
            db_type: Database type (full, light)
            output_dir: Output directory
            force: Overwrite existing

        Returns:
            Path to database directory
        """
        if db_type not in BAKTA_RELEASES:
            logger.error(f"Unknown Bakta database type: {db_type}")
            return None

        info = BAKTA_RELEASES[db_type]

        if output_dir is None:
            output_dir = self.base_dir / 'bakta' / db_type
        else:
            output_dir = Path(output_dir)

        if output_dir.exists() and not force:
            logger.warning(f"Database already exists at {output_dir}")
            return output_dir

        logger.info(f"Downloading Bakta {info['name']}...")
        logger.info(f"Size: ~{info['size_gb']} GB")

        # Use Bakta's download command inside conda environment
        try:
            cmd = ['bakta_db', 'download', '--output', str(output_dir)]
            if db_type == 'light':
                cmd.append('--type')
                cmd.append('light')

            env_name = self._get_conda_env_for_tool('bakta_db')
            if env_name:
                self._run_in_conda_env(cmd, env_name)
            else:
                subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to download Bakta database: {e}")
            return None
        except FileNotFoundError:
            logger.error("bakta_db not found. Please ensure the conda environment is set up.")
            logger.info("Try: mamba env create -f workflows/envs/annotation.yaml")
            return None

        # Register
        db_info = DatabaseInfo(
            name=f'bakta-{db_type}',
            tool='bakta',
            path=str(output_dir),
            version=info['version'],
            installed=datetime.now().isoformat(),
            size_bytes=self._get_dir_size(output_dir),
            verified=True
        )
        self._register_database(db_info)

        logger.info(f"Bakta database installed at: {output_dir}")
        return output_dir

    # ============ Verification ============

    def verify_database(self, tool: str, name: str = None) -> bool:
        """
        Verify database integrity using manifest checksums.

        Args:
            tool: Tool name (emu, centrifuger, etc.)
            name: Database name (optional)

        Returns:
            True if verification passes
        """
        db = self.get_database(tool, name)
        if not db:
            logger.error(f"Database not found: {tool}/{name}")
            return False

        db_path = Path(db.path)
        manifest_file = db_path / 'manifest.json'

        if not manifest_file.exists():
            logger.warning(f"No manifest found for {tool}/{name}")
            return False

        manifest = DatabaseManifest.load(manifest_file)

        logger.info(f"Verifying {tool}/{name}...")

        # Check overall checksum
        current_checksum = self._calculate_checksum(db_path)
        if current_checksum != manifest.checksum:
            logger.error("Database checksum mismatch!")
            return False

        logger.info("Database verification passed")
        return True
