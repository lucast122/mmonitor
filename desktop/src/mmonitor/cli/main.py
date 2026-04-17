#!/usr/bin/env python3
"""
MMonitor CLI v0.3
New click-based command line interface with Snakemake integration.

Usage:
    mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p project1
    mmonitor run taxonomy-16s --sample-sheet samples.csv -p project1
    mmonitor run assembly -i reads.fastq -s sample1
    mmonitor database emu download --preset silva138
    mmonitor config show
    mmonitor status --project project1
"""
import os
import sys
import csv
import json
import time
import subprocess
import click
import keyring
import logging
from pathlib import Path
from datetime import datetime

import yaml
import requests as req_lib
from mmonitor.pipeline.runner import SnakemakeRunner

# Keyring service name for storing MMonitor server password
KEYRING_SERVICE = 'mmonitor'

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('mmonitor')

# Version
__version__ = '0.3.0'


class Config:
    """CLI configuration management.

    Uses the same nested Snakemake-compatible YAML schema as the GUI
    (PipelineConfig). The config file is ``~/.mmonitor/config.yaml`` (YAML).
    A legacy ``config.json`` is migrated on first load.
    """

    def __init__(self):
        self.config_dir = Path.home() / '.mmonitor'
        self.config_file = self.config_dir / 'config.yaml'
        self._legacy_json = self.config_dir / 'config.json'
        self.credentials_file = self.config_dir / 'credentials.json'
        self._config = None

    @property
    def config(self) -> dict:
        if self._config is None:
            self._config = self.load()
        return self._config

    def load(self) -> dict:
        """Load configuration from YAML file, migrating from JSON if needed."""
        if self.config_file.exists():
            with open(self.config_file) as f:
                return yaml.safe_load(f) or self.defaults()
        # Migrate legacy JSON if it exists
        if self._legacy_json.exists():
            try:
                with open(self._legacy_json) as f:
                    legacy = json.load(f)
                # Best-effort migration: merge legacy values into defaults
                cfg = self.defaults()
                cfg = self._deep_merge(cfg, legacy)
                self.save(cfg)
                return cfg
            except Exception:
                pass
        return self.defaults()

    def save(self, config: dict):
        """Save configuration to YAML file."""
        self.config_dir.mkdir(parents=True, exist_ok=True)
        with open(self.config_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        self._config = config

    def load_config_file(self, path: str) -> dict:
        """Load an external YAML config file (e.g. exported from GUI)."""
        with open(path) as f:
            return yaml.safe_load(f) or {}

    @staticmethod
    def _deep_merge(base: dict, override: dict) -> dict:
        """Recursively merge override into base."""
        result = dict(base)
        for k, v in override.items():
            if k in result and isinstance(result[k], dict) and isinstance(v, dict):
                result[k] = Config._deep_merge(result[k], v)
            else:
                result[k] = v
        return result

    @staticmethod
    def defaults() -> dict:
        """Default configuration matching the Snakemake config schema."""
        return {
            'threads': os.cpu_count() or 4,
            'memory_gb': 16,
            'output_dir': str(Path.home() / 'mmonitor_results'),
            'db_base_dir': '',  # Empty = ~/mmonitor_databases; set to e.g. /mnt/disk2/db

            'server': {
                'url': 'http://localhost:8000',
                'username': '',
                'password': '',  # leave empty; use 'mmonitor config set-password' to store securely
                'upload_results': True,
            },

            'filtlong': {
                'enabled': True,
                'min_length': 1000,
                'min_mean_q': 7.0,
                'keep_percent': 90.0,
                'target_bases': None,
                'max_length': None,
            },

            'emu': {
                'database': '',
                'min_abundance': 0.0001,
            },

            'centrifuger': {
                'database': '',
                'min_hitlen': 22,
            },

            'flye': {
                'mode': 'nano-raw',
                'meta': True,
                'min_overlap': 1000,
            },

            'medaka': {
                'model': 'r1041_e82_400bps_hac_v5.0.0',
            },

            'metabat2': {
                'min_contig': 2500,
            },

            'checkm2': {
                'database': '',
            },

            'gtdbtk': {
                'database': '',
            },

            'bakta': {
                'database': '',
                'min_contig_length': 1,
            },

            'eggnog': {
                'database': '',
            },

            'realtime': {
                'min_files_for_auto_analysis': 5,
            },
        }


# Global config instance
config = Config()


# ============ Sample Sheet Parsing ============
def parse_sample_sheet(path: str) -> dict:
    """Parse a sample sheet CSV/TSV with columns: sample_name, file_path.

    Returns dict: {sample_name: {"fastq": [file_path, ...]}, ...}
    Multiple rows with the same sample_name aggregate their file paths.
    """
    samples = {}
    with open(path) as f:
        sample_text = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample_text, delimiters=',\t')
        except csv.Error:
            dialect = 'excel'  # default to comma-separated
        reader = csv.DictReader(f, dialect=dialect)
        for row in reader:
            name = row.get('sample_name', '').strip()
            fpath = row.get('file_path', '').strip()
            if not name or not fpath:
                continue
            if not os.path.exists(fpath):
                raise click.BadParameter(
                    f"File not found for sample '{name}': {fpath}",
                    param_hint="--sample-sheet"
                )
            if name in samples:
                samples[name]['fastq'].append(fpath)
            else:
                samples[name] = {'fastq': [fpath]}
    if not samples:
        raise click.BadParameter(
            f"No valid samples found in {path}. "
            "Expected columns: sample_name, file_path",
            param_hint="--sample-sheet"
        )
    return samples


def _build_samples(sample_sheet, sample, input_files):
    """Build samples dict from CLI args. Raises click.UsageError on invalid input."""
    if sample_sheet:
        if sample or input_files:
            raise click.UsageError(
                "Cannot use --sample-sheet together with -s/--sample or -i/--input"
            )
        return parse_sample_sheet(sample_sheet)
    if not sample or not input_files:
        raise click.UsageError(
            "Either --sample-sheet or both -s/--sample and -i/--input are required"
        )
    return {sample: {"fastq": list(input_files)}}


# ============ Credential Management ============
def _get_password(username: str) -> str:
    """Get password from keyring, or prompt interactively if not stored."""
    password = keyring.get_password(KEYRING_SERVICE, username)
    if password:
        return password
    click.echo(f"No stored password found for user '{username}'.")
    password = click.prompt("Server password", hide_input=True)
    return password


def _verify_credentials(server_url: str, username: str, password: str) -> bool:
    """Verify credentials against the MMonitor server. Returns True on success."""
    if not server_url.startswith('http'):
        server_url = f"https://{server_url}"

    login_url = f"{server_url}/api/v1/auth/login/"
    try:
        resp = req_lib.post(login_url, json={
            'username': username,
            'password': password,
        }, timeout=10)
        if resp.status_code == 200:
            return True
        click.echo(click.style(
            f"Login failed (HTTP {resp.status_code}): {resp.text}", fg='red'), err=True)
        return False
    except req_lib.ConnectionError:
        click.echo(click.style(
            f"Cannot connect to server at {server_url}", fg='red'), err=True)
        return False
    except req_lib.Timeout:
        click.echo(click.style(
            f"Server at {server_url} timed out", fg='red'), err=True)
        return False


def _ensure_credentials(cfg: dict) -> dict:
    """Ensure valid server credentials are present in cfg.

    - Reads username from cfg['server']['username'] (prompts if empty).
    - Reads password from keyring (prompts if not stored).
    - Verifies credentials against the server.
    - On success, stores the password in keyring and sets cfg['server']['password'].
    - On failure, exits.
    """
    server_cfg = cfg.setdefault('server', {})
    server_url = server_cfg.get('url', 'http://localhost:8000')
    upload = server_cfg.get('upload_results', True)

    if not upload:
        click.echo("Upload disabled, skipping credential check.")
        return cfg

    username = server_cfg.get('username', '')
    if not username:
        username = click.prompt("MMonitor server username")
        server_cfg['username'] = username

    password = _get_password(username)

    click.echo(f"Verifying credentials for '{username}' at {server_url}...")
    if not _verify_credentials(server_url, username, password):
        # Give one retry with manual entry
        click.echo("Please try again.")
        password = click.prompt("Server password", hide_input=True)
        if not _verify_credentials(server_url, username, password):
            click.echo(click.style("Authentication failed. Aborting.", fg='red'), err=True)
            sys.exit(1)

    # Store validated password in keyring for future use
    keyring.set_password(KEYRING_SERVICE, username, password)
    click.echo(click.style("Credentials verified.", fg='green'))

    # Pass password to Snakemake config so upload scripts can use it
    server_cfg['password'] = password
    return cfg


# ============ Main CLI Group ============
@click.group()
@click.version_option(__version__, prog_name='mmonitor')
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose output')
@click.option('-q', '--quiet', is_flag=True, help='Suppress non-error output')
@click.pass_context
def cli(ctx, verbose, quiet):
    """
    MMonitor - Metagenomic Monitoring Tool

    Analyze nanopore sequencing data for taxonomy, assembly, and functional analysis.

    \b
    Quick start (single sample):
      mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p myproject

    \b
    Quick start (multiple samples via sample sheet):
      mmonitor run taxonomy-16s --sample-sheet samples.csv -p myproject

    \b
    Sample sheet format (CSV or TSV):
      sample_name,file_path
      sample1,/path/to/reads1.fastq
      sample2,/path/to/reads2.fastq

    \b
    For more information on a command:
      mmonitor run --help
      mmonitor database --help
    """
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    ctx.obj['quiet'] = quiet

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    elif quiet:
        logging.getLogger().setLevel(logging.ERROR)


# ============ Run Command Group ============
@cli.group()
@click.pass_context
def run(ctx):
    """
    Run analysis pipelines.

    \b
    Available pipelines:
      taxonomy-16s   16S rRNA taxonomy using EMU
      taxonomy-wgs   Whole genome taxonomy using Centrifuger
      assembly       Assembly + binning + taxonomy + annotation + upload
      assembly-full  Assembly + full functional analysis + upload

    \b
    Input modes:
      Single sample:  -s sample1 -i reads.fastq
      Sample sheet:   --sample-sheet samples.csv
    """
    pass


# Common options for run commands
def common_run_options(f):
    """Common options for all run commands."""
    f = click.option('-i', '--input', 'input_files', multiple=True,
                     type=click.Path(exists=True),
                     help='Input FASTQ file(s). Use with -s for single sample.')(f)
    f = click.option('-s', '--sample', default='',
                     help='Sample name. Use with -i for single sample.')(f)
    f = click.option('--sample-sheet', type=click.Path(exists=True),
                     help='CSV/TSV sample sheet with columns: sample_name, file_path. '
                          'Alternative to -s/-i for processing multiple samples.')(f)
    f = click.option('-p', '--project', required=True,
                     help='Project name')(f)
    f = click.option('-u', '--subproject', default='',
                     help='Subproject name')(f)
    f = click.option('-d', '--date', default='',
                     help='Sample date (YYYY-MM-DD)')(f)
    f = click.option('-o', '--output', 'output_dir',
                     type=click.Path(),
                     help='Output directory (default: ~/mmonitor_results)')(f)
    f = click.option('-t', '--threads', type=int,
                     help='Number of threads')(f)
    f = click.option('-c', '--config-file', type=click.Path(exists=True),
                     help='Path to config file')(f)
    f = click.option('--dry-run', is_flag=True,
                     help='Show what would be done without executing')(f)
    f = click.option('--docker', is_flag=True, default=False,
                     help='Run pipeline inside Docker (auto-enabled on Windows)')(f)
    return f


def run_snakemake(target: str, config_dict: dict, dry_run: bool = False,
                  cores: int = 4, verbose: bool = False,
                  use_docker: bool = False) -> int:
    """
    Run a Snakemake workflow via the unified SnakemakeRunner.

    Args:
        target: Snakemake target rule
        config_dict: Configuration dictionary
        dry_run: If True, only show what would be done
        cores: Number of cores to use
        verbose: Enable verbose output
        use_docker: Force Docker mode

    Returns:
        Exit code (0 for success)
    """
    runner = SnakemakeRunner(use_docker=use_docker if use_docker else None)
    return runner.run(
        target, config_dict,
        threads=cores, dry_run=dry_run, verbose=verbose,
    )


def _build_cfg(config_file, sample_sheet, sample, input_files, project,
               subproject, date, output_dir, threads):
    """Build base config dict with samples, shared across all run commands."""
    if config_file:
        cfg = config.load_config_file(config_file)
    else:
        cfg = config.config.copy()

    # Build samples dict
    samples = _build_samples(sample_sheet, sample, input_files)
    cfg['samples'] = samples

    # Overlay CLI options
    t = threads or cfg.get('threads', 4)
    cfg['project_name'] = project
    cfg['subproject_name'] = subproject
    cfg['sample_date'] = date or datetime.now().strftime('%Y-%m-%d')
    cfg['output_dir'] = output_dir or cfg.get('output_dir', 'results')
    cfg['threads'] = t

    # Report what we're running
    sample_names = list(samples.keys())
    if len(sample_names) == 1:
        click.echo(f"Sample: {sample_names[0]}")
    else:
        click.echo(f"Samples ({len(sample_names)}): {', '.join(sample_names)}")

    # Verify server credentials before running the pipeline
    cfg = _ensure_credentials(cfg)

    return cfg


@run.command('taxonomy-16s')
@common_run_options
@click.option('--emu-db', type=click.Path(exists=True),
              help='Path to EMU database')
@click.option('--min-abundance', type=float, default=0.0001,
              help='Minimum abundance threshold')
@click.option('--use-filtlong/--no-filtlong', default=True,
              help='Enable/disable Filtlong pre-filtering')
@click.pass_context
def run_taxonomy_16s(ctx, input_files, sample, sample_sheet, project, subproject, date,
                     output_dir, threads, config_file, dry_run, docker,
                     emu_db, min_abundance, use_filtlong):
    """
    Run 16S rRNA taxonomy analysis using EMU.

    \b
    Single sample:
      mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p project1
      mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p project1 --emu-db /path/to/db

    \b
    Multiple samples:
      mmonitor run taxonomy-16s --sample-sheet samples.csv -p project1
    """
    click.echo("Running 16S taxonomy analysis")

    cfg = _build_cfg(config_file, sample_sheet, sample, input_files,
                     project, subproject, date, output_dir, threads)

    t = cfg['threads']

    if 'filtlong' not in cfg:
        cfg['filtlong'] = {}
    cfg['filtlong']['enabled'] = use_filtlong

    if 'emu' not in cfg:
        cfg['emu'] = {}
    if emu_db:
        cfg['emu']['database'] = emu_db
    cfg['emu']['min_abundance'] = min_abundance
    cfg['emu']['threads'] = t

    # Validate EMU database
    if not cfg['emu'].get('database'):
        click.echo("Error: EMU database path required. Use --emu-db or configure in ~/.mmonitor/config.yaml", err=True)
        click.echo("To download a database: mmonitor database emu download --preset silva138", err=True)
        sys.exit(1)

    # Run Snakemake
    exit_code = run_snakemake(
        'taxonomy_16s',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False),
        use_docker=docker,
    )

    if exit_code == 0:
        click.echo(click.style("Analysis complete!", fg='green'))
    else:
        click.echo(click.style("Analysis failed!", fg='red'), err=True)

    sys.exit(exit_code)


@run.command('taxonomy-wgs')
@common_run_options
@click.option('--centrifuger-db', type=str,
              help='Path to Centrifuger database prefix (e.g. /path/to/cfr_gtdb_r226)')
@click.option('--min-hitlen', type=int, default=22,
              help='Minimum hit length')
@click.option('--use-filtlong/--no-filtlong', default=True,
              help='Enable/disable Filtlong pre-filtering')
@click.pass_context
def run_taxonomy_wgs(ctx, input_files, sample, sample_sheet, project, subproject, date,
                     output_dir, threads, config_file, dry_run, docker,
                     centrifuger_db, min_hitlen, use_filtlong):
    """
    Run whole-genome taxonomy analysis using Centrifuger.

    \b
    Single sample:
      mmonitor run taxonomy-wgs -i reads.fastq -s sample1 -p project1

    \b
    Multiple samples:
      mmonitor run taxonomy-wgs --sample-sheet samples.csv -p project1
    """
    click.echo("Running WGS taxonomy analysis")

    cfg = _build_cfg(config_file, sample_sheet, sample, input_files,
                     project, subproject, date, output_dir, threads)

    t = cfg['threads']

    if 'filtlong' not in cfg:
        cfg['filtlong'] = {}
    cfg['filtlong']['enabled'] = use_filtlong

    if 'centrifuger' not in cfg:
        cfg['centrifuger'] = {}
    if centrifuger_db:
        cfg['centrifuger']['database'] = centrifuger_db
    cfg['centrifuger']['min_hitlen'] = min_hitlen
    cfg['centrifuger']['threads'] = t

    # Validate Centrifuger database
    if not cfg['centrifuger'].get('database'):
        click.echo("Error: Centrifuger database path required. Use --centrifuger-db or configure in ~/.mmonitor/config.yaml", err=True)
        click.echo("To download a database: mmonitor database centrifuger download --preset bacteria", err=True)
        sys.exit(1)

    # Run Snakemake
    exit_code = run_snakemake(
        'taxonomy_wgs',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False),
        use_docker=docker,
    )

    if exit_code == 0:
        click.echo(click.style("Analysis complete!", fg='green'))
    else:
        click.echo(click.style("Analysis failed!", fg='red'), err=True)

    sys.exit(exit_code)


@run.command('assembly')
@common_run_options
@click.option('--mode', type=click.Choice(['nano-raw', 'nano-corr', 'nano-hq']),
              default='nano-raw', help='Flye assembly mode')
@click.option('--medaka-model', default='r1041_e82_400bps_hac_v5.0.0',
              help='Medaka polishing model')
@click.option('--meta/--no-meta', default=True,
              help='Enable/disable metagenomic mode')
@click.option('--gtdbtk-db', type=click.Path(exists=True),
              help='Path to GTDB-TK database')
@click.option('--bakta-db', type=click.Path(exists=True),
              help='Path to Bakta database')
@click.option('--checkm2-db', type=click.Path(exists=True),
              help='Path to CheckM2 database')
@click.pass_context
def run_assembly(ctx, input_files, sample, sample_sheet, project, subproject, date,
                 output_dir, threads, config_file, dry_run, docker,
                 mode, medaka_model, meta, gtdbtk_db, bakta_db, checkm2_db):
    """
    Assembly + binning + taxonomy + annotation pipeline.

    \b
    Runs: Flye assembly -> Medaka polishing -> MetaBAT2 binning ->
          CheckM2 quality -> GTDB-TK taxonomy -> Bakta annotation -> upload

    \b
    Single sample:
      mmonitor run assembly -i reads.fastq -s sample1 -p project1

    \b
    Multiple samples:
      mmonitor run assembly --sample-sheet samples.csv -p project1
    """
    click.echo("Running assembly pipeline")

    cfg = _build_cfg(config_file, sample_sheet, sample, input_files,
                     project, subproject, date, output_dir, threads)

    t = cfg['threads']

    if 'flye' not in cfg:
        cfg['flye'] = {}
    cfg['flye']['mode'] = mode
    cfg['flye']['meta'] = meta
    cfg['flye']['threads'] = t

    if 'medaka' not in cfg:
        cfg['medaka'] = {}
    cfg['medaka']['model'] = medaka_model
    cfg['medaka']['threads'] = t

    for tool, cli_db in [('checkm2', checkm2_db), ('gtdbtk', gtdbtk_db),
                         ('bakta', bakta_db)]:
        if tool not in cfg:
            cfg[tool] = {}
        if cli_db:
            cfg[tool]['database'] = cli_db
        cfg[tool].setdefault('database', '')
        cfg[tool]['threads'] = t

    # Run Snakemake
    exit_code = run_snakemake(
        'assembly',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False),
        use_docker=docker,
    )

    if exit_code == 0:
        click.echo(click.style("Assembly pipeline complete!", fg='green'))
    else:
        click.echo(click.style("Assembly pipeline failed!", fg='red'), err=True)

    sys.exit(exit_code)


@run.command('assembly-full')
@common_run_options
@click.option('--mode', type=click.Choice(['nano-raw', 'nano-corr', 'nano-hq']),
              default='nano-raw', help='Flye assembly mode')
@click.option('--medaka-model', default='r1041_e82_400bps_hac_v5.0.0',
              help='Medaka polishing model')
@click.option('--meta/--no-meta', default=True,
              help='Enable/disable metagenomic mode')
@click.option('--gtdbtk-db', type=click.Path(exists=True),
              help='Path to GTDB-TK database')
@click.option('--bakta-db', type=click.Path(exists=True),
              help='Path to Bakta database')
@click.option('--checkm2-db', type=click.Path(exists=True),
              help='Path to CheckM2 database')
@click.option('--eggnog-db', type=click.Path(exists=True),
              help='Path to eggNOG database')
@click.pass_context
def run_assembly_full(ctx, input_files, sample, sample_sheet, project, subproject, date,
                      output_dir, threads, config_file, dry_run, docker,
                      mode, medaka_model, meta, gtdbtk_db, bakta_db,
                      checkm2_db, eggnog_db):
    """
    Full assembly + functional analysis pipeline.

    \b
    Runs everything in 'assembly' plus:
    eggNOG-mapper (COG/KEGG) -> InterProScan (PFAM) -> upload MAGs

    \b
    Single sample:
      mmonitor run assembly-full -i reads.fastq -s sample1 -p project1

    \b
    Multiple samples:
      mmonitor run assembly-full --sample-sheet samples.csv -p project1
    """
    click.echo("Running full assembly + functional pipeline")

    cfg = _build_cfg(config_file, sample_sheet, sample, input_files,
                     project, subproject, date, output_dir, threads)

    t = cfg['threads']

    if 'flye' not in cfg:
        cfg['flye'] = {}
    cfg['flye']['mode'] = mode
    cfg['flye']['meta'] = meta
    cfg['flye']['threads'] = t

    if 'medaka' not in cfg:
        cfg['medaka'] = {}
    cfg['medaka']['model'] = medaka_model
    cfg['medaka']['threads'] = t

    for tool, cli_db in [('checkm2', checkm2_db), ('gtdbtk', gtdbtk_db),
                         ('bakta', bakta_db), ('eggnog', eggnog_db)]:
        if tool not in cfg:
            cfg[tool] = {}
        if cli_db:
            cfg[tool]['database'] = cli_db
        cfg[tool].setdefault('database', '')
        if tool != 'eggnog':
            cfg[tool]['threads'] = t

    # Run Snakemake
    exit_code = run_snakemake(
        'assembly_full',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False),
        use_docker=docker,
    )

    if exit_code == 0:
        click.echo(click.style("Full assembly pipeline complete!", fg='green'))
    else:
        click.echo(click.style("Full assembly pipeline failed!", fg='red'), err=True)

    sys.exit(exit_code)


# ============ Database Command Group ============

# Initialize database manager
from mmonitor.cli.database_manager import DatabaseManager, EMU_PRESETS, CENTRIFUGER_PRESETS, CENTRIFUGER_PREBUILT
db_manager = DatabaseManager()


@cli.group()
def database():
    """
    Manage analysis databases.

    \b
    Commands:
      emu          Manage EMU databases
      centrifuger  Manage Centrifuger databases
      gtdbtk       Manage GTDB-TK database
      checkm2      Manage CheckM2 database
      bakta        Manage Bakta database
      list         List installed databases
      verify       Verify database integrity
    """
    pass


# ============ EMU Database Commands ============
@database.group('emu')
def database_emu():
    """Manage EMU databases for 16S analysis."""
    pass


@database_emu.command('download')
@click.option('--preset', type=click.Choice(['silva138', 'rdp', 'default']),
              required=True, help='Download preset database')
@click.option('--output', '-o', type=click.Path(),
              help='Output directory')
@click.option('--force', '-f', is_flag=True,
              help='Overwrite existing database')
def emu_download(preset, output, force):
    """
    Download EMU database.

    \b
    Available presets:
      default   - Default EMU database (~100 MB)
      silva138  - SILVA 138.2 SSU database (~500 MB)
      rdp       - RDP database (~200 MB)

    \b
    Example:
      mmonitor database emu download --preset silva138
      mmonitor database emu download --preset silva138 -o /data/databases/emu
    """
    output_dir = Path(output) if output else None

    result = db_manager.download_emu_database(preset, output_dir, force)

    if result:
        click.echo(click.style(f"EMU database installed at: {result}", fg='green'))

        # Update config
        cfg = config.config
        cfg.setdefault('emu', {})['database'] = str(result)
        config.save(cfg)
        click.echo(f"Configuration updated with database path")
    else:
        click.echo(click.style("Download failed", fg='red'), err=True)
        sys.exit(1)


# ============ Centrifuger Database Commands ============
@database.group('centrifuger')
def database_centrifuger():
    """
    Manage Centrifuger databases for WGS analysis.

    \b
    Commands:
      download        Download prebuilt index (recommended, fast)
      download-build  Download sequences from NCBI and build index (slow)
      build           Build index from local FASTA files
    """
    pass


@database_centrifuger.command('download')
@click.option('--preset', type=click.Choice(['gtdb', 'nt']),
              required=True, help='Prebuilt index to download')
@click.option('--output', '-o', type=click.Path(),
              help='Output directory')
@click.option('--force', '-f', is_flag=True,
              help='Overwrite existing database')
def centrifuger_download(preset, output, force):
    """
    Download prebuilt Centrifuger index (recommended).

    Downloads a ready-to-use prebuilt index from Dropbox. Much faster than
    building from sequences.

    \b
    Available presets:
      gtdb  - GTDB r226 (~176 GB, 2025-05-20)
      nt    - NCBI core nt (~212 GB, 2025-05-20)

    \b
    Note: Check https://github.com/mourisl/centrifuger for newer prebuilt indices.

    \b
    Example:
      mmonitor database centrifuger download --preset gtdb
      mmonitor database centrifuger download --preset nt -o /data/databases/centrifuger
    """
    info = CENTRIFUGER_PREBUILT[preset]
    output_dir = Path(output) if output else None

    click.echo(f"This will download the prebuilt {info['name']} index.")
    click.echo(f"Size: ~{info['size_gb']} GB")
    click.echo(f"Note: check https://github.com/mourisl/centrifuger for newer prebuilt indices.")

    if not click.confirm("Continue?"):
        click.echo("Aborted.")
        return

    result = db_manager.download_centrifuger_prebuilt(preset, output_dir, force)

    if result:
        click.echo(click.style(f"Prebuilt Centrifuger index installed at: {result}", fg='green'))

        # Update config
        cfg = config.config
        cfg.setdefault('centrifuger', {})['database'] = str(result)
        config.save(cfg)
        click.echo("Configuration updated with database path")
    else:
        click.echo(click.style("Download failed", fg='red'), err=True)
        sys.exit(1)


@database_centrifuger.command('download-build')
@click.option('--preset', type=click.Choice(['bacteria', 'archaea', 'viral', 'fungi']),
              required=True, help='Download and build preset database')
@click.option('--output', '-o', type=click.Path(),
              help='Output directory')
@click.option('--threads', '-t', type=int, default=4,
              help='Download threads')
@click.option('--force', '-f', is_flag=True,
              help='Overwrite existing database')
def centrifuger_download_build(preset, output, threads, force):
    """
    Download sequences from NCBI and build Centrifuger index.

    This downloads raw genomes and builds the index locally. For most users,
    'mmonitor database centrifuger download' (prebuilt) is faster and easier.

    \b
    Available presets:
      bacteria  - RefSeq bacterial genomes (~50 GB)
      archaea   - RefSeq archaeal genomes (~5 GB)
      viral     - RefSeq viral genomes (~2 GB)
      fungi     - RefSeq fungal genomes (~10 GB)

    \b
    Warning: This downloads genomes from NCBI and builds an index.
    It may take several hours depending on your connection.

    \b
    Example:
      mmonitor database centrifuger download-build --preset bacteria
      mmonitor database centrifuger download-build --preset bacteria -t 8
    """
    output_dir = Path(output) if output else None

    click.echo(f"This will download and build the {preset} database.")
    click.echo("This process may take several hours.")
    click.echo("Tip: use 'mmonitor database centrifuger download' for prebuilt indices instead.")

    if not click.confirm("Continue?"):
        click.echo("Aborted.")
        return

    result = db_manager.download_centrifuger_database(preset, output_dir, threads, force)

    if result:
        click.echo(click.style(f"Centrifuger database built at: {result}", fg='green'))

        # Update config
        cfg = config.config
        cfg.setdefault('centrifuger', {})['database'] = str(result)
        config.save(cfg)
        click.echo("Configuration updated with database path")
    else:
        click.echo(click.style("Download/build failed", fg='red'), err=True)
        sys.exit(1)


@database_centrifuger.command('build')
@click.option('--fasta-dir', '-f', required=True, type=click.Path(exists=True),
              help='Directory with genome FASTA files')
@click.option('--taxonomy-dir', '-t', required=True, type=click.Path(exists=True),
              help='Directory with NCBI taxonomy files (nodes.dmp, names.dmp)')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output database prefix')
@click.option('--name', '-n', required=True,
              help='Database name')
@click.option('--threads', type=int, default=4,
              help='Threads for building')
def centrifuger_build(fasta_dir, taxonomy_dir, output, name, threads):
    """
    Build custom Centrifuger database.

    \b
    Example:
      mmonitor database centrifuger build -f /data/genomes -t /data/taxonomy -o /data/mydb -n custom
    """
    result = db_manager.build_centrifuger_database(
        Path(fasta_dir), Path(taxonomy_dir), Path(output), name, threads
    )

    if result:
        click.echo(click.style(f"Centrifuger database built at: {result}", fg='green'))
    else:
        click.echo(click.style("Build failed", fg='red'), err=True)
        sys.exit(1)


# ============ GTDB-TK Database Commands ============
@database.group('gtdbtk')
def database_gtdbtk():
    """Manage GTDB-TK database for MAG taxonomy."""
    pass


@database_gtdbtk.command('download')
@click.option('--release', type=click.Choice(['r220', 'r214']),
              default='r220', help='GTDB release version')
@click.option('--output', '-o', type=click.Path(),
              help='Output directory')
@click.option('--force', '-f', is_flag=True,
              help='Overwrite existing database')
def gtdbtk_download(release, output, force):
    """
    Download GTDB-TK database.

    \b
    Warning: This is a large download (~85 GB for R220).

    \b
    Example:
      mmonitor database gtdbtk download --release r220
    """
    output_dir = Path(output) if output else None

    click.echo(f"Downloading GTDB-TK {release} database...")
    click.echo("This is a large download (~85 GB) and may take several hours.")

    if not click.confirm("Continue?"):
        click.echo("Aborted.")
        return

    result = db_manager.download_gtdbtk_database(release, output_dir, force)

    if result:
        click.echo(click.style(f"GTDB-TK database installed at: {result}", fg='green'))
        click.echo(f"\nSet environment variable:")
        click.echo(f"  export GTDBTK_DATA_PATH={result}")

        # Update config
        cfg = config.config
        cfg.setdefault('gtdbtk', {})['database'] = str(result)
        config.save(cfg)
    else:
        click.echo(click.style("Download failed", fg='red'), err=True)
        sys.exit(1)


# ============ CheckM2 Database Commands ============
@database.group('checkm2')
def database_checkm2():
    """Manage CheckM2 database for MAG quality."""
    pass


@database_checkm2.command('download')
@click.option('--output', '-o', type=click.Path(),
              help='Output directory')
@click.option('--force', '-f', is_flag=True,
              help='Overwrite existing database')
def checkm2_download(output, force):
    """
    Download CheckM2 database.

    \b
    Example:
      mmonitor database checkm2 download
    """
    output_dir = Path(output) if output else None

    result = db_manager.download_checkm2_database(output_dir, force)

    if result:
        click.echo(click.style(f"CheckM2 database installed at: {result}", fg='green'))

        # Update config
        cfg = config.config
        cfg.setdefault('checkm2', {})['database'] = str(result)
        config.save(cfg)
    else:
        click.echo(click.style("Download failed", fg='red'), err=True)
        sys.exit(1)


# ============ Bakta Database Commands ============
@database.group('bakta')
def database_bakta():
    """Manage Bakta database for gene annotation."""
    pass


@database_bakta.command('download')
@click.option('--type', 'db_type', type=click.Choice(['full', 'light']),
              default='light', help='Database type')
@click.option('--output', '-o', type=click.Path(),
              help='Output directory')
@click.option('--force', '-f', is_flag=True,
              help='Overwrite existing database')
def bakta_download(db_type, output, force):
    """
    Download Bakta database.

    \b
    Database types:
      light - Smaller database (~3 GB)
      full  - Complete database (~60 GB)

    \b
    Example:
      mmonitor database bakta download --type light
    """
    output_dir = Path(output) if output else None

    result = db_manager.download_bakta_database(db_type, output_dir, force)

    if result:
        click.echo(click.style(f"Bakta database installed at: {result}", fg='green'))

        # Update config
        cfg = config.config
        cfg.setdefault('bakta', {})['database'] = str(result)
        config.save(cfg)
    else:
        click.echo(click.style("Download failed", fg='red'), err=True)
        sys.exit(1)


# ============ Database List & Verify ============
@database.command('list')
def database_list():
    """List all installed databases."""
    click.echo("\n" + "=" * 60)
    click.echo("MMonitor Installed Databases")
    click.echo("=" * 60)

    registry = db_manager.list_databases()

    if not registry:
        click.echo("\nNo databases installed yet.")
        click.echo("\nTo install databases:")
        click.echo("  mmonitor database emu download --preset silva138")
        click.echo("  mmonitor database centrifuger download --preset bacteria")
        click.echo()
        return

    for tool, databases in registry.items():
        click.echo(f"\n{tool.upper()}:")
        click.echo("-" * 40)

        for db in databases:
            status = click.style("OK", fg='green') if db.verified else click.style("?", fg='yellow')
            size_mb = db.size_bytes / (1024 * 1024)
            click.echo(f"  {status} {db.name}")
            click.echo(f"      Version:   {db.version}")
            click.echo(f"      Path:      {db.path}")
            click.echo(f"      Size:      {size_mb:.1f} MB")
            click.echo(f"      Installed: {db.installed[:10]}")

    click.echo()


@database.command('verify')
@click.option('--tool', '-t', required=True,
              type=click.Choice(['emu', 'centrifuger', 'gtdbtk', 'checkm2', 'bakta']),
              help='Tool database to verify')
@click.option('--name', '-n', help='Database name (optional)')
def database_verify(tool, name):
    """
    Verify database integrity using checksums.

    \b
    Example:
      mmonitor database verify -t emu
      mmonitor database verify -t emu -n emu-silva138
    """
    click.echo(f"Verifying {tool} database...")

    if db_manager.verify_database(tool, name):
        click.echo(click.style("Database verification passed", fg='green'))
    else:
        click.echo(click.style("Database verification failed", fg='red'), err=True)
        sys.exit(1)


# ============ Config Command Group ============
@cli.group('config')
def config_group():
    """Manage MMonitor configuration."""
    pass


@config_group.command('show')
def config_show():
    """Show current configuration (passwords are masked)."""
    import copy
    cfg = copy.deepcopy(config.config)
    # Mask password if present
    if 'server' in cfg and cfg['server'].get('password'):
        cfg['server']['password'] = '********'
    # Show keyring status
    username = cfg.get('server', {}).get('username', '')
    if username:
        has_keyring_pw = keyring.get_password(KEYRING_SERVICE, username) is not None
        cfg['server']['_password_in_keyring'] = has_keyring_pw
    click.echo(yaml.dump(cfg, default_flow_style=False, sort_keys=False))


@config_group.command('set')
@click.argument('key')
@click.argument('value')
def config_set(key, value):
    """
    Set a configuration value.

    \b
    Example:
      mmonitor config set threads 8
      mmonitor config set emu.database /path/to/emu/db
    """
    cfg = config.config

    # Handle nested keys
    keys = key.split('.')
    current = cfg
    for k in keys[:-1]:
        if k not in current:
            current[k] = {}
        current = current[k]

    # Try to parse value as JSON, otherwise use as string
    try:
        parsed_value = json.loads(value)
    except json.JSONDecodeError:
        parsed_value = value

    current[keys[-1]] = parsed_value
    config.save(cfg)

    click.echo(f"Set {key} = {parsed_value}")


@config_group.command('init')
@click.option('--force', is_flag=True, help='Overwrite existing config')
def config_init(force):
    """Initialize configuration with defaults."""
    if config.config_file.exists() and not force:
        click.echo(f"Config already exists at {config.config_file}")
        click.echo("Use --force to overwrite")
        return

    config.save(config.defaults())
    click.echo(f"Created config at {config.config_file}")


@config_group.command('set-password')
@click.option('-u', '--username', default='',
              help='Username (if not set, reads from config)')
def config_set_password(username):
    """
    Store server password securely in the system keyring.

    The password is never saved in config files. It is stored in your
    OS keyring (e.g. GNOME Keyring, macOS Keychain, Windows Credential Locker).

    \b
    Example:
      mmonitor config set-password
      mmonitor config set-password -u myuser
    """
    cfg = config.config
    if not username:
        username = cfg.get('server', {}).get('username', '')
    if not username:
        username = click.prompt("Username")
        # Save username to config
        cfg.setdefault('server', {})['username'] = username
        config.save(cfg)

    password = click.prompt("Password", hide_input=True, confirmation_prompt=True)

    server_url = cfg.get('server', {}).get('url', 'http://localhost:8000')
    click.echo(f"Verifying credentials against {server_url}...")
    if _verify_credentials(server_url, username, password):
        keyring.set_password(KEYRING_SERVICE, username, password)
        click.echo(click.style(
            f"Password stored in system keyring for user '{username}'.", fg='green'))
    else:
        click.echo(click.style("Password NOT stored (verification failed).", fg='red'), err=True)
        sys.exit(1)


# ============ Status Command ============
@cli.command()
@click.option('-p', '--project', help='Filter by project')
@click.option('-s', '--sample', help='Filter by sample')
def status(project, sample):
    """Check status of running or completed analyses."""
    click.echo("Analysis Status")
    click.echo("-" * 50)

    # Check output directory for results
    output_dir = Path(config.config.get('output_dir', 'results'))

    if not output_dir.exists():
        click.echo("No analyses found.")
        return

    for sample_dir in sorted(output_dir.iterdir()):
        if sample_dir.is_dir():
            if sample and sample != sample_dir.name:
                continue

            click.echo(f"\nSample: {sample_dir.name}")

            # Check for pipeline outputs
            pipelines = {
                'taxonomy_16s': sample_dir / 'taxonomy_16s' / 'emu_results.json',
                'taxonomy_wgs': sample_dir / 'taxonomy_wgs' / 'centrifuger_results.json',
                'assembly': sample_dir / 'assembly' / 'polished' / 'consensus.fasta',
                'binning': sample_dir / 'binning' / 'checkm2' / 'quality_report.tsv',
                'annotation': sample_dir / 'annotation' / 'annotation_complete.json'
            }

            for pipeline, output_file in pipelines.items():
                if output_file.exists():
                    status = click.style("complete", fg='green')
                else:
                    status = click.style("not run", fg='yellow')
                click.echo(f"  {pipeline:15} {status}")


# ============ Watch Command ============
@cli.command()
@click.argument('folder', type=click.Path(exists=True))
@click.option('-p', '--project', required=True, help='Project name')
@click.option('-u', '--subproject', default='', help='Subproject name')
@click.option('-a', '--analysis', 'analysis_type',
              type=click.Choice(['taxonomy-16s', 'taxonomy-wgs', 'assembly', 'assembly-full']),
              default='taxonomy-16s', help='Analysis pipeline to run')
@click.option('--interval', type=int, default=300,
              help='Seconds between analysis checks (default: 300)')
@click.option('--min-files', type=int, default=5,
              help='Minimum new files per sample before triggering analysis (default: 5)')
@click.option('-o', '--output', 'output_dir', type=click.Path(),
              help='Output directory (default: ~/mmonitor_results)')
@click.option('-t', '--threads', type=int, help='Number of threads')
@click.option('-c', '--config-file', type=click.Path(exists=True),
              help='Path to config file')
@click.option('--preload-index', is_flag=True, default=False,
              help='Preload centrifuger index into page cache at start (recommended for taxonomy-wgs)')
@click.option('--emu-db', type=click.Path(exists=True), help='Path to EMU database')
@click.option('--centrifuger-db', type=str, help='Path to Centrifuger database')
@click.option('--use-filtlong/--no-filtlong', default=True,
              help='Enable/disable Filtlong pre-filtering')
@click.pass_context
def watch(ctx, folder, project, subproject, analysis_type, interval, min_files,
          output_dir, threads, config_file, preload_index, emu_db,
          centrifuger_db, use_filtlong):
    """
    Watch a folder for new FASTQ files and analyze in real-time.

    Monitors FOLDER recursively for new sequencing reads. Files are grouped
    by sample (from barcode subdirectory names). When enough new files
    accumulate for a sample, they are concatenated and the analysis pipeline
    is triggered automatically.

    \b
    The index is loaded once per analysis batch, not per file. For large
    databases (e.g. GTDB), use --preload-index to warm the page cache at
    start for near-instant subsequent loads.

    \b
    Examples:
      mmonitor watch /data/sequencing/run1/fastq_pass -p myproject -a taxonomy-16s
      mmonitor watch /data/run1/fastq_pass -p myproject -a taxonomy-wgs --preload-index --interval 120
      mmonitor watch /data/run1/fastq_pass -p myproject --min-files 10 --interval 600
    """
    from mmonitor.pipeline.watcher import SequencingWatcher, preload_index as do_preload

    # Build base config
    if config_file:
        cfg = config.load_config_file(config_file)
    else:
        cfg = config.config.copy()

    t = threads or cfg.get('threads', 4)
    cfg['project_name'] = project
    cfg['subproject_name'] = subproject
    cfg['sample_date'] = datetime.now().strftime('%Y-%m-%d')
    cfg['output_dir'] = output_dir or cfg.get('output_dir', str(Path.home() / 'mmonitor_results'))
    cfg['threads'] = t

    if 'filtlong' not in cfg:
        cfg['filtlong'] = {}
    cfg['filtlong']['enabled'] = use_filtlong

    # Tool-specific config
    if analysis_type == 'taxonomy-16s':
        if 'emu' not in cfg:
            cfg['emu'] = {}
        if emu_db:
            cfg['emu']['database'] = emu_db
        if not cfg['emu'].get('database'):
            click.echo("Error: EMU database required. Use --emu-db or configure in ~/.mmonitor/config.yaml", err=True)
            sys.exit(1)
        target = 'taxonomy_16s'

    elif analysis_type == 'taxonomy-wgs':
        if 'centrifuger' not in cfg:
            cfg['centrifuger'] = {}
        if centrifuger_db:
            cfg['centrifuger']['database'] = centrifuger_db
        if not cfg['centrifuger'].get('database'):
            click.echo("Error: Centrifuger database required. Use --centrifuger-db or configure in ~/.mmonitor/config.yaml", err=True)
            sys.exit(1)
        target = 'taxonomy_wgs'

    elif analysis_type == 'assembly':
        target = 'assembly'
    elif analysis_type == 'assembly-full':
        target = 'assembly_full'
    else:
        target = analysis_type.replace('-', '_')

    # Verify credentials before starting the long-running watch
    cfg = _ensure_credentials(cfg)

    verbose = ctx.obj.get('verbose', False)

    # Preload centrifuger index into page cache if requested
    if preload_index and analysis_type == 'taxonomy-wgs':
        db_path = cfg.get('centrifuger', {}).get('database', '')
        if db_path:
            click.echo(f"Preloading centrifuger index into page cache...")
            do_preload(db_path)
            click.echo("Index preloaded.")

    # Analysis callback: triggered by the watcher when samples are ready
    def run_analysis(samples_dict: dict, concat_dir: str):
        """Run Snakemake for pending samples."""
        run_cfg = dict(cfg)
        run_cfg['samples'] = samples_dict

        click.echo(f"\n--- Analysis triggered at {datetime.now().strftime('%H:%M:%S')} ---")
        sample_names = list(samples_dict.keys())
        click.echo(f"Samples: {', '.join(sample_names)}")

        exit_code = run_snakemake(
            target, run_cfg,
            cores=t, verbose=verbose,
        )

        if exit_code != 0:
            raise RuntimeError(f"Snakemake exited with code {exit_code}")

    # Status callback
    def on_status(msg: str):
        timestamp = datetime.now().strftime('%H:%M:%S')
        click.echo(f"[{timestamp}] {msg}")

    # Create and start watcher
    concat_dir = str(Path(cfg['output_dir']) / '_concat')
    watcher = SequencingWatcher(
        watch_dir=folder,
        analysis_callback=run_analysis,
        min_files=min_files,
        interval=interval,
        concat_dir=concat_dir,
        on_status=on_status,
    )

    click.echo(f"MMonitor Real-Time Watch")
    click.echo(f"  Folder:    {folder}")
    click.echo(f"  Project:   {project}")
    click.echo(f"  Pipeline:  {analysis_type}")
    click.echo(f"  Interval:  {interval}s")
    click.echo(f"  Min files: {min_files}")
    click.echo(f"  Output:    {cfg['output_dir']}")
    click.echo(f"Press Ctrl+C to stop.\n")

    watcher.start()

    try:
        # Keep main thread alive
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        click.echo("\nStopping watcher...")
        watcher.stop()
        click.echo("Done.")


# ============ GUI Command ============
@cli.command()
def gui():
    """Launch the MMonitor graphical interface."""
    try:
        from mmonitor.userside.view import GUI
    except ImportError as e:
        click.echo(click.style(
            f"Error: GUI dependencies not installed. Install with: pip install mmonitor[gui]",
            fg='red'), err=True)
        click.echo(f"Details: {e}", err=True)
        sys.exit(1)

    app = GUI()
    app.mainloop()


# ============ Serve Command ============
@cli.command()
@click.option('--host', default='127.0.0.1', help='Host to bind to')
@click.option('--port', default=8000, type=int, help='Port to bind to')
@click.option('--no-browser', is_flag=True, help='Do not open browser automatically')
@click.pass_context
def serve(ctx, host, port, no_browser):
    """Start the MMonitor local server and open the dashboard in a browser."""
    import signal
    import time

    # Find server directory by walking up from this file
    _here = Path(__file__).resolve().parent
    server_dir = None
    for parent in [_here] + list(_here.parents):
        candidate = parent / 'server'
        if (candidate / 'manage.py').exists():
            server_dir = candidate
            break
    if server_dir is None:
        click.echo(click.style("Error: Server directory not found", fg='red'), err=True)
        sys.exit(1)

    manage_py = server_dir / 'manage.py'
    if not manage_py.exists():
        click.echo(click.style(f"Error: manage.py not found at {manage_py}", fg='red'), err=True)
        sys.exit(1)

    dashboard_url = f'http://{host}:{port}/dashboard/'
    click.echo(f"Starting MMonitor server at http://{host}:{port}/")
    click.echo(f"Dashboard: {dashboard_url}")
    click.echo("Press Ctrl+C to stop the server")

    env = os.environ.copy()
    env['PYTHONUNBUFFERED'] = '1'
    env['DJANGO_SETTINGS_MODULE'] = 'MMonitor.settings'

    proc = subprocess.Popen(
        [sys.executable, str(manage_py), 'runserver', f'{host}:{port}'],
        cwd=str(server_dir),
        env=env,
    )

    # Wait for server to be ready
    import requests as req_lib
    ready = False
    for attempt in range(30):
        try:
            resp = req_lib.get(f'http://{host}:{port}/', timeout=1)
            ready = True
            break
        except Exception:
            time.sleep(1)

    if ready and not no_browser:
        import webbrowser
        webbrowser.open(dashboard_url)
    elif not ready:
        click.echo(click.style("Warning: Server may not have started yet", fg='yellow'))

    # Wait for the process to finish (Ctrl+C will send SIGINT)
    try:
        proc.wait()
    except KeyboardInterrupt:
        click.echo("\nShutting down server...")
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()


# ============ Entry Point ============
def main():
    """Main entry point."""
    cli()


if __name__ == '__main__':
    main()
