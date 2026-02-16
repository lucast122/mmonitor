#!/usr/bin/env python3
"""
MMonitor CLI v0.3
New click-based command line interface with Snakemake integration.

Usage:
    mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p project1
    mmonitor run assembly -i reads.fastq -s sample1
    mmonitor database emu download --preset silva138
    mmonitor config show
    mmonitor status --project project1
"""
import os
import sys
import json
import subprocess
import click
import logging
from pathlib import Path
from datetime import datetime

import yaml
from mmonitor.pipeline.runner import SnakemakeRunner

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
    Quick start:
      mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p myproject
      mmonitor run assembly -i reads.fastq -s sample1 -p myproject

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
      assembly       Metagenomic assembly with Flye + Medaka
      binning        MAG binning with MetaBAT2 + CheckM2
      annotation     MAG annotation with GTDB-TK + Bakta
      functional     Full functional analysis pipeline
    """
    pass


# Common options for run commands
def common_run_options(f):
    """Common options for all run commands."""
    f = click.option('-i', '--input', 'input_files', required=True, multiple=True,
                     type=click.Path(exists=True),
                     help='Input FASTQ file(s)')(f)
    f = click.option('-s', '--sample', required=True,
                     help='Sample name')(f)
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


@run.command('taxonomy-16s')
@common_run_options
@click.option('--emu-db', type=click.Path(exists=True),
              help='Path to EMU database')
@click.option('--min-abundance', type=float, default=0.0001,
              help='Minimum abundance threshold')
@click.option('--use-filtlong/--no-filtlong', default=True,
              help='Enable/disable Filtlong pre-filtering')
@click.pass_context
def run_taxonomy_16s(ctx, input_files, sample, project, subproject, date,
                     output_dir, threads, config_file, dry_run, docker,
                     emu_db, min_abundance, use_filtlong):
    """
    Run 16S rRNA taxonomy analysis using EMU.

    \b
    Example:
      mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p project1
      mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p project1 --emu-db /path/to/db
    """
    click.echo(f"Running 16S taxonomy analysis for sample: {sample}")

    # Build config: start with config-file or defaults
    if config_file:
        cfg = config.load_config_file(config_file)
    else:
        cfg = config.config.copy()

    # Overlay CLI options
    t = threads or cfg.get('threads', 4)
    cfg['sample_name'] = sample
    cfg['project_name'] = project
    cfg['subproject_name'] = subproject
    cfg['sample_date'] = date or datetime.now().strftime('%Y-%m-%d')
    cfg['input_fastq'] = list(input_files)
    cfg['output_dir'] = output_dir or cfg.get('output_dir', 'results')
    cfg['threads'] = t

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
def run_taxonomy_wgs(ctx, input_files, sample, project, subproject, date,
                     output_dir, threads, config_file, dry_run, docker,
                     centrifuger_db, min_hitlen, use_filtlong):
    """
    Run whole-genome taxonomy analysis using Centrifuger.

    \b
    Example:
      mmonitor run taxonomy-wgs -i reads.fastq -s sample1 -p project1
      mmonitor run taxonomy-wgs -i reads.fastq -s sample1 -p project1 --centrifuger-db /path/to/db
    """
    click.echo(f"Running WGS taxonomy analysis for sample: {sample}")

    # Build config: start with config-file or defaults
    if config_file:
        cfg = config.load_config_file(config_file)
    else:
        cfg = config.config.copy()

    # Overlay CLI options
    t = threads or cfg.get('threads', 4)
    cfg['sample_name'] = sample
    cfg['project_name'] = project
    cfg['subproject_name'] = subproject
    cfg['sample_date'] = date or datetime.now().strftime('%Y-%m-%d')
    cfg['input_fastq'] = list(input_files)
    cfg['output_dir'] = output_dir or cfg.get('output_dir', 'results')
    cfg['threads'] = t

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
@click.pass_context
def run_assembly(ctx, input_files, sample, project, subproject, date,
                 output_dir, threads, config_file, dry_run, docker,
                 mode, medaka_model, meta):
    """
    Run metagenomic assembly with Flye and Medaka polishing.

    \b
    Example:
      mmonitor run assembly -i reads.fastq -s sample1 -p project1
      mmonitor run assembly -i reads.fastq -s sample1 -p project1 --mode nano-hq
    """
    click.echo(f"Running assembly for sample: {sample}")

    # Build config: start with config-file or defaults
    if config_file:
        cfg = config.load_config_file(config_file)
    else:
        cfg = config.config.copy()

    # Overlay CLI options
    t = threads or cfg.get('threads', 4)
    cfg['sample_name'] = sample
    cfg['project_name'] = project
    cfg['subproject_name'] = subproject
    cfg['sample_date'] = date or datetime.now().strftime('%Y-%m-%d')
    cfg['input_fastq'] = list(input_files)
    cfg['output_dir'] = output_dir or cfg.get('output_dir', 'results')
    cfg['threads'] = t

    if 'flye' not in cfg:
        cfg['flye'] = {}
    cfg['flye']['mode'] = mode
    cfg['flye']['meta'] = meta
    cfg['flye']['threads'] = t

    if 'medaka' not in cfg:
        cfg['medaka'] = {}
    cfg['medaka']['model'] = medaka_model
    cfg['medaka']['threads'] = t

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
        click.echo(click.style("Assembly complete!", fg='green'))
    else:
        click.echo(click.style("Assembly failed!", fg='red'), err=True)

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
def run_assembly_full(ctx, input_files, sample, project, subproject, date,
                      output_dir, threads, config_file, dry_run, docker,
                      mode, medaka_model, meta, gtdbtk_db, bakta_db,
                      checkm2_db, eggnog_db):
    """
    Run the full WGS assembly pipeline: assembly + binning + annotation + functional + upload.

    Chains: Flye assembly -> Medaka polishing -> MetaBAT2 binning -> CheckM2 quality ->
    GTDB-TK taxonomy -> Bakta annotation -> eggNOG-mapper (COG/KEGG) ->
    InterProScan (PFAM) -> Upload to server.

    \b
    Example:
      mmonitor run assembly-full -i reads.fastq -s sample1 -p project1
      mmonitor run assembly-full -i reads.fastq -s sample1 -p project1 --mode nano-hq
    """
    click.echo(f"Running full assembly pipeline for sample: {sample}")

    # Build config: start with config-file or defaults
    if config_file:
        cfg = config.load_config_file(config_file)
    else:
        cfg = config.config.copy()

    # Overlay CLI options
    t = threads or cfg.get('threads', 4)
    cfg['sample_name'] = sample
    cfg['project_name'] = project
    cfg['subproject_name'] = subproject
    cfg['sample_date'] = date or datetime.now().strftime('%Y-%m-%d')
    cfg['input_fastq'] = list(input_files)
    cfg['output_dir'] = output_dir or cfg.get('output_dir', 'results')
    cfg['threads'] = t

    if 'flye' not in cfg:
        cfg['flye'] = {}
    cfg['flye']['mode'] = mode
    cfg['flye']['meta'] = meta
    cfg['flye']['threads'] = t

    if 'medaka' not in cfg:
        cfg['medaka'] = {}
    cfg['medaka']['model'] = medaka_model
    cfg['medaka']['threads'] = t

    # Override database paths only if CLI flags are given
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


@run.command('binning')
@common_run_options
@click.option('--min-contig', type=int, default=2500,
              help='Minimum contig length for binning')
@click.option('--checkm2-db', type=click.Path(exists=True),
              help='Path to CheckM2 database')
@click.pass_context
def run_binning(ctx, input_files, sample, project, subproject, date,
                output_dir, threads, config_file, dry_run, docker,
                min_contig, checkm2_db):
    """
    Run MAG binning with MetaBAT2 and quality assessment with CheckM2.

    Note: Requires assembly output. Run 'mmonitor run assembly' first.

    \b
    Example:
      mmonitor run binning -i reads.fastq -s sample1 -p project1
    """
    click.echo(f"Running binning for sample: {sample}")

    # Build config: start with config-file or defaults
    if config_file:
        cfg = config.load_config_file(config_file)
    else:
        cfg = config.config.copy()

    # Overlay CLI options
    t = threads or cfg.get('threads', 4)
    cfg['sample_name'] = sample
    cfg['project_name'] = project
    cfg['subproject_name'] = subproject
    cfg['sample_date'] = date or datetime.now().strftime('%Y-%m-%d')
    cfg['input_fastq'] = list(input_files)
    cfg['output_dir'] = output_dir or cfg.get('output_dir', 'results')
    cfg['threads'] = t

    if 'metabat2' not in cfg:
        cfg['metabat2'] = {}
    cfg['metabat2']['min_contig'] = min_contig

    if 'checkm2' not in cfg:
        cfg['checkm2'] = {}
    if checkm2_db:
        cfg['checkm2']['database'] = checkm2_db
    cfg['checkm2'].setdefault('database', '')
    cfg['checkm2']['threads'] = t

    # Run Snakemake
    exit_code = run_snakemake(
        'binning',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False),
        use_docker=docker,
    )

    if exit_code == 0:
        click.echo(click.style("Binning complete!", fg='green'))
    else:
        click.echo(click.style("Binning failed!", fg='red'), err=True)

    sys.exit(exit_code)


@run.command('functional')
@common_run_options
@click.option('--gtdbtk-db', type=click.Path(exists=True),
              help='Path to GTDB-TK database')
@click.option('--bakta-db', type=click.Path(exists=True),
              help='Path to Bakta database')
@click.pass_context
def run_functional(ctx, input_files, sample, project, subproject, date,
                   output_dir, threads, config_file, dry_run, docker,
                   gtdbtk_db, bakta_db):
    """
    Run full functional analysis pipeline.

    Includes: Assembly → Binning → Quality → Taxonomy → Annotation

    \b
    Example:
      mmonitor run functional -i reads.fastq -s sample1 -p project1
    """
    click.echo(f"Running functional analysis for sample: {sample}")

    # Build config: start with config-file or defaults
    if config_file:
        cfg = config.load_config_file(config_file)
    else:
        cfg = config.config.copy()

    # Overlay CLI options
    t = threads or cfg.get('threads', 4)
    cfg['sample_name'] = sample
    cfg['project_name'] = project
    cfg['subproject_name'] = subproject
    cfg['sample_date'] = date or datetime.now().strftime('%Y-%m-%d')
    cfg['input_fastq'] = list(input_files)
    cfg['output_dir'] = output_dir or cfg.get('output_dir', 'results')
    cfg['threads'] = t

    if 'gtdbtk' not in cfg:
        cfg['gtdbtk'] = {}
    if gtdbtk_db:
        cfg['gtdbtk']['database'] = gtdbtk_db
    cfg['gtdbtk'].setdefault('database', '')
    cfg['gtdbtk']['threads'] = t

    if 'bakta' not in cfg:
        cfg['bakta'] = {}
    if bakta_db:
        cfg['bakta']['database'] = bakta_db
    cfg['bakta'].setdefault('database', '')
    cfg['bakta']['threads'] = t

    # Run Snakemake
    exit_code = run_snakemake(
        'functional',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False),
        use_docker=docker,
    )

    if exit_code == 0:
        click.echo(click.style("Functional analysis complete!", fg='green'))
    else:
        click.echo(click.style("Functional analysis failed!", fg='red'), err=True)

    sys.exit(exit_code)


# ============ Database Command Group ============

# Initialize database manager
from mmonitor.cli.database_manager import DatabaseManager, EMU_PRESETS, CENTRIFUGER_PRESETS
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
      silva138  - SILVA 138.1 SSU database (~500 MB)
      rdp       - RDP 11.5 database (~200 MB)
      default   - Default EMU database (~100 MB)

    \b
    Example:
      mmonitor database emu download --preset silva138
      mmonitor database emu download --preset silva138 -o /data/databases/emu
    """
    output_dir = Path(output) if output else None

    result = db_manager.download_emu_database(preset, output_dir, force)

    if result:
        click.echo(click.style(f"✓ EMU database installed at: {result}", fg='green'))

        # Update config
        cfg = config.config
        cfg.setdefault('emu', {})['database'] = str(result)
        config.save(cfg)
        click.echo(f"Configuration updated with database path")
    else:
        click.echo(click.style("✗ Download failed", fg='red'), err=True)
        sys.exit(1)


# ============ Centrifuger Database Commands ============
@database.group('centrifuger')
def database_centrifuger():
    """Manage Centrifuger databases for WGS analysis."""
    pass


@database_centrifuger.command('download')
@click.option('--preset', type=click.Choice(['bacteria', 'archaea', 'viral', 'fungi']),
              required=True, help='Download preset database')
@click.option('--output', '-o', type=click.Path(),
              help='Output directory')
@click.option('--threads', '-t', type=int, default=4,
              help='Download threads')
@click.option('--force', '-f', is_flag=True,
              help='Overwrite existing database')
def centrifuger_download(preset, output, threads, force):
    """
    Download and build Centrifuger database from NCBI.

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
      mmonitor database centrifuger download --preset bacteria
      mmonitor database centrifuger download --preset bacteria -t 8
    """
    output_dir = Path(output) if output else None

    click.echo(f"This will download and build the {preset} database.")
    click.echo("This process may take several hours.")

    if not click.confirm("Continue?"):
        click.echo("Aborted.")
        return

    result = db_manager.download_centrifuger_database(preset, output_dir, threads, force)

    if result:
        click.echo(click.style(f"✓ Centrifuger database built at: {result}", fg='green'))

        # Update config
        cfg = config.config
        cfg.setdefault('centrifuger', {})['database'] = str(result)
        config.save(cfg)
        click.echo(f"Configuration updated with database path")
    else:
        click.echo(click.style("✗ Download/build failed", fg='red'), err=True)
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
        click.echo(click.style(f"✓ Centrifuger database built at: {result}", fg='green'))
    else:
        click.echo(click.style("✗ Build failed", fg='red'), err=True)
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
        click.echo(click.style(f"✓ GTDB-TK database installed at: {result}", fg='green'))
        click.echo(f"\nSet environment variable:")
        click.echo(f"  export GTDBTK_DATA_PATH={result}")

        # Update config
        cfg = config.config
        cfg.setdefault('gtdbtk', {})['database'] = str(result)
        config.save(cfg)
    else:
        click.echo(click.style("✗ Download failed", fg='red'), err=True)
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
        click.echo(click.style(f"✓ CheckM2 database installed at: {result}", fg='green'))

        # Update config
        cfg = config.config
        cfg.setdefault('checkm2', {})['database'] = str(result)
        config.save(cfg)
    else:
        click.echo(click.style("✗ Download failed", fg='red'), err=True)
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
        click.echo(click.style(f"✓ Bakta database installed at: {result}", fg='green'))

        # Update config
        cfg = config.config
        cfg.setdefault('bakta', {})['database'] = str(result)
        config.save(cfg)
    else:
        click.echo(click.style("✗ Download failed", fg='red'), err=True)
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
            status = click.style("✓", fg='green') if db.verified else click.style("?", fg='yellow')
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
        click.echo(click.style("✓ Database verification passed", fg='green'))
    else:
        click.echo(click.style("✗ Database verification failed", fg='red'), err=True)
        sys.exit(1)


# ============ Config Command Group ============
@cli.group('config')
def config_group():
    """Manage MMonitor configuration."""
    pass


@config_group.command('show')
def config_show():
    """Show current configuration."""
    click.echo(yaml.dump(config.config, default_flow_style=False, sort_keys=False))


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
                    status = click.style("✓ complete", fg='green')
                else:
                    status = click.style("○ not run", fg='yellow')
                click.echo(f"  {pipeline:15} {status}")


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
