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
import shutil
import subprocess
import click
import logging
from pathlib import Path
from datetime import datetime

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('mmonitor')

# Version
__version__ = '0.3.0'

# Find project paths
def get_project_root() -> Path:
    """Find the MMonitor project root directory."""
    # Try to find relative to this file
    current = Path(__file__).resolve().parent
    while current != current.parent:
        if (current / 'workflows' / 'Snakefile').exists():
            return current
        if (current / 'desktop' / 'workflows' / 'Snakefile').exists():
            return current / 'desktop'
        current = current.parent
    # Fallback
    return Path(__file__).resolve().parent.parent.parent

PROJECT_ROOT = get_project_root()
WORKFLOWS_DIR = PROJECT_ROOT / 'workflows'
CONFIG_DIR = PROJECT_ROOT / 'workflows' / 'config'


class Config:
    """CLI configuration management."""

    def __init__(self):
        self.config_dir = Path.home() / '.mmonitor'
        self.config_file = self.config_dir / 'config.json'
        self.credentials_file = self.config_dir / 'credentials.json'
        self._config = None

    @property
    def config(self) -> dict:
        if self._config is None:
            self._config = self.load()
        return self._config

    def load(self) -> dict:
        """Load configuration from file."""
        if self.config_file.exists():
            with open(self.config_file) as f:
                return json.load(f)
        return self.defaults()

    def save(self, config: dict):
        """Save configuration to file."""
        self.config_dir.mkdir(parents=True, exist_ok=True)
        with open(self.config_file, 'w') as f:
            json.dump(config, f, indent=2)
        self._config = config

    @staticmethod
    def defaults() -> dict:
        """Default configuration."""
        return {
            'threads': os.cpu_count() or 4,
            'output_dir': str(Path.home() / 'mmonitor_results'),
            'server': {
                'url': 'http://localhost:8000',
                'username': '',
                'upload_results': True
            },
            'databases': {
                'emu': '',
                'centrifuger': '',
                'gtdbtk': '',
                'checkm2': '',
                'bakta': ''
            },
            'filtlong': {
                'enabled': True,
                'min_length': 1000,
                'min_mean_q': 7.0,
                'keep_percent': 90.0
            }
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
    return f


def run_snakemake(target: str, config_dict: dict, dry_run: bool = False,
                  cores: int = 4, verbose: bool = False) -> int:
    """
    Run a Snakemake workflow.

    Args:
        target: Snakemake target rule
        config_dict: Configuration dictionary
        dry_run: If True, only show what would be done
        cores: Number of cores to use
        verbose: Enable verbose output

    Returns:
        Exit code (0 for success)
    """
    import subprocess
    import tempfile

    # Write config to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        import yaml
        yaml.dump(config_dict, f)
        config_file = f.name

    try:
        cmd = [
            'snakemake',
            '--snakefile', str(WORKFLOWS_DIR / 'Snakefile'),
            '--configfile', config_file,
            '--cores', str(cores),
            target
        ]

        if dry_run:
            cmd.append('--dry-run')

        if verbose:
            cmd.append('--verbose')

        # Use conda environments
        cmd.extend(['--use-conda', '--conda-frontend', 'mamba'])

        logger.info(f"Running: {' '.join(cmd)}")

        result = subprocess.run(cmd, check=False)
        return result.returncode

    finally:
        # Cleanup temp config
        os.unlink(config_file)


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
                     output_dir, threads, config_file, dry_run,
                     emu_db, min_abundance, use_filtlong):
    """
    Run 16S rRNA taxonomy analysis using EMU.

    \b
    Example:
      mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p project1
      mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p project1 --emu-db /path/to/db
    """
    click.echo(f"Running 16S taxonomy analysis for sample: {sample}")

    # Build config
    cfg = config.config.copy()
    cfg.update({
        'sample_name': sample,
        'project_name': project,
        'subproject_name': subproject,
        'sample_date': date or datetime.now().strftime('%Y-%m-%d'),
        'input_fastq': list(input_files),
        'output_dir': output_dir or cfg.get('output_dir', 'results'),
        'threads': threads or cfg.get('threads', 4),
        'filtlong': {
            'enabled': use_filtlong,
            **cfg.get('filtlong', {})
        },
        'emu': {
            'database': emu_db or cfg.get('databases', {}).get('emu', ''),
            'min_abundance': min_abundance,
            'threads': threads or cfg.get('threads', 4)
        }
    })

    # Validate EMU database
    if not cfg['emu']['database']:
        click.echo("Error: EMU database path required. Use --emu-db or configure in ~/.mmonitor/config.json", err=True)
        click.echo("To download a database: mmonitor database emu download --preset silva138", err=True)
        sys.exit(1)

    # Run Snakemake
    exit_code = run_snakemake(
        'taxonomy_16s',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False)
    )

    if exit_code == 0:
        click.echo(click.style("✓ Analysis complete!", fg='green'))
    else:
        click.echo(click.style("✗ Analysis failed!", fg='red'), err=True)

    sys.exit(exit_code)


@run.command('taxonomy-wgs')
@common_run_options
@click.option('--centrifuger-db', type=click.Path(exists=True),
              help='Path to Centrifuger database prefix')
@click.option('--min-hitlen', type=int, default=22,
              help='Minimum hit length')
@click.option('--use-filtlong/--no-filtlong', default=True,
              help='Enable/disable Filtlong pre-filtering')
@click.pass_context
def run_taxonomy_wgs(ctx, input_files, sample, project, subproject, date,
                     output_dir, threads, config_file, dry_run,
                     centrifuger_db, min_hitlen, use_filtlong):
    """
    Run whole-genome taxonomy analysis using Centrifuger.

    \b
    Example:
      mmonitor run taxonomy-wgs -i reads.fastq -s sample1 -p project1
      mmonitor run taxonomy-wgs -i reads.fastq -s sample1 -p project1 --centrifuger-db /path/to/db
    """
    click.echo(f"Running WGS taxonomy analysis for sample: {sample}")

    # Build config
    cfg = config.config.copy()
    cfg.update({
        'sample_name': sample,
        'project_name': project,
        'subproject_name': subproject,
        'sample_date': date or datetime.now().strftime('%Y-%m-%d'),
        'input_fastq': list(input_files),
        'output_dir': output_dir or cfg.get('output_dir', 'results'),
        'threads': threads or cfg.get('threads', 4),
        'filtlong': {
            'enabled': use_filtlong,
            **cfg.get('filtlong', {})
        },
        'centrifuger': {
            'database': centrifuger_db or cfg.get('databases', {}).get('centrifuger', ''),
            'min_hitlen': min_hitlen,
            'threads': threads or cfg.get('threads', 4)
        }
    })

    # Validate Centrifuger database
    if not cfg['centrifuger']['database']:
        click.echo("Error: Centrifuger database path required. Use --centrifuger-db or configure in ~/.mmonitor/config.json", err=True)
        click.echo("To download a database: mmonitor database centrifuger download --preset bacteria", err=True)
        sys.exit(1)

    # Run Snakemake
    exit_code = run_snakemake(
        'taxonomy_wgs',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False)
    )

    if exit_code == 0:
        click.echo(click.style("✓ Analysis complete!", fg='green'))
    else:
        click.echo(click.style("✗ Analysis failed!", fg='red'), err=True)

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
                 output_dir, threads, config_file, dry_run,
                 mode, medaka_model, meta):
    """
    Run metagenomic assembly with Flye and Medaka polishing.

    \b
    Example:
      mmonitor run assembly -i reads.fastq -s sample1 -p project1
      mmonitor run assembly -i reads.fastq -s sample1 -p project1 --mode nano-hq
    """
    click.echo(f"Running assembly for sample: {sample}")

    # Build config
    cfg = config.config.copy()
    cfg.update({
        'sample_name': sample,
        'project_name': project,
        'subproject_name': subproject,
        'sample_date': date or datetime.now().strftime('%Y-%m-%d'),
        'input_fastq': list(input_files),
        'output_dir': output_dir or cfg.get('output_dir', 'results'),
        'threads': threads or cfg.get('threads', 4),
        'flye': {
            'mode': mode,
            'meta': meta,
            'threads': threads or cfg.get('threads', 4)
        },
        'medaka': {
            'model': medaka_model,
            'threads': threads or cfg.get('threads', 4)
        }
    })

    # Run Snakemake
    exit_code = run_snakemake(
        'assembly',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False)
    )

    if exit_code == 0:
        click.echo(click.style("✓ Assembly complete!", fg='green'))
    else:
        click.echo(click.style("✗ Assembly failed!", fg='red'), err=True)

    sys.exit(exit_code)


@run.command('binning')
@common_run_options
@click.option('--min-contig', type=int, default=2500,
              help='Minimum contig length for binning')
@click.option('--checkm2-db', type=click.Path(exists=True),
              help='Path to CheckM2 database')
@click.pass_context
def run_binning(ctx, input_files, sample, project, subproject, date,
                output_dir, threads, config_file, dry_run,
                min_contig, checkm2_db):
    """
    Run MAG binning with MetaBAT2 and quality assessment with CheckM2.

    Note: Requires assembly output. Run 'mmonitor run assembly' first.

    \b
    Example:
      mmonitor run binning -i reads.fastq -s sample1 -p project1
    """
    click.echo(f"Running binning for sample: {sample}")

    # Build config
    cfg = config.config.copy()
    cfg.update({
        'sample_name': sample,
        'project_name': project,
        'subproject_name': subproject,
        'sample_date': date or datetime.now().strftime('%Y-%m-%d'),
        'input_fastq': list(input_files),
        'output_dir': output_dir or cfg.get('output_dir', 'results'),
        'threads': threads or cfg.get('threads', 4),
        'metabat2': {
            'min_contig': min_contig
        },
        'checkm2': {
            'database': checkm2_db or cfg.get('databases', {}).get('checkm2', ''),
            'threads': threads or cfg.get('threads', 4)
        }
    })

    # Run Snakemake
    exit_code = run_snakemake(
        'binning',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False)
    )

    if exit_code == 0:
        click.echo(click.style("✓ Binning complete!", fg='green'))
    else:
        click.echo(click.style("✗ Binning failed!", fg='red'), err=True)

    sys.exit(exit_code)


@run.command('functional')
@common_run_options
@click.option('--gtdbtk-db', type=click.Path(exists=True),
              help='Path to GTDB-TK database')
@click.option('--bakta-db', type=click.Path(exists=True),
              help='Path to Bakta database')
@click.pass_context
def run_functional(ctx, input_files, sample, project, subproject, date,
                   output_dir, threads, config_file, dry_run,
                   gtdbtk_db, bakta_db):
    """
    Run full functional analysis pipeline.

    Includes: Assembly → Binning → Quality → Taxonomy → Annotation

    \b
    Example:
      mmonitor run functional -i reads.fastq -s sample1 -p project1
    """
    click.echo(f"Running functional analysis for sample: {sample}")

    # Build config
    cfg = config.config.copy()
    cfg.update({
        'sample_name': sample,
        'project_name': project,
        'subproject_name': subproject,
        'sample_date': date or datetime.now().strftime('%Y-%m-%d'),
        'input_fastq': list(input_files),
        'output_dir': output_dir or cfg.get('output_dir', 'results'),
        'threads': threads or cfg.get('threads', 4),
        'gtdbtk': {
            'database': gtdbtk_db or cfg.get('databases', {}).get('gtdbtk', ''),
            'threads': threads or cfg.get('threads', 4)
        },
        'bakta': {
            'database': bakta_db or cfg.get('databases', {}).get('bakta', ''),
            'threads': threads or cfg.get('threads', 4)
        }
    })

    # Run Snakemake
    exit_code = run_snakemake(
        'functional',
        cfg,
        dry_run=dry_run,
        cores=cfg['threads'],
        verbose=ctx.obj.get('verbose', False)
    )

    if exit_code == 0:
        click.echo(click.style("✓ Functional analysis complete!", fg='green'))
    else:
        click.echo(click.style("✗ Functional analysis failed!", fg='red'), err=True)

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
        cfg.setdefault('databases', {})['emu'] = str(result)
        config.save(cfg)
        click.echo(f"Configuration updated with database path")
    else:
        click.echo(click.style("✗ Download failed", fg='red'), err=True)
        sys.exit(1)


@database_emu.command('build')
@click.option('--fasta', '-f', required=True, type=click.Path(exists=True),
              help='Input FASTA file with sequences')
@click.option('--taxonomy', '-t', required=True, type=click.Path(exists=True),
              help='Taxonomy TSV file (columns: tax_id, species, genus, family, etc.)')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output directory')
@click.option('--name', '-n', required=True,
              help='Database name')
@click.option('--threads', type=int, default=4,
              help='Threads for indexing')
def emu_build(fasta, taxonomy, output, name, threads):
    """
    Build custom EMU database.

    \b
    Input files:
      - FASTA file with 16S sequences (headers should contain tax_id)
      - TSV file with taxonomy mapping

    \b
    Example:
      mmonitor database emu build -f sequences.fasta -t taxonomy.tsv -o /data/db -n mydb
    """
    result = db_manager.build_emu_database(
        Path(fasta), Path(taxonomy), Path(output), name, threads
    )

    if result:
        click.echo(click.style(f"✓ EMU database built at: {result}", fg='green'))
    else:
        click.echo(click.style("✗ Build failed", fg='red'), err=True)
        sys.exit(1)


@database_emu.command('restrict')
@click.option('--source', '-s', required=True, type=click.Path(exists=True),
              help='Source EMU database directory')
@click.option('--taxids', '-t', required=True,
              help='Comma-separated list of taxonomy IDs to include')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output directory')
@click.option('--name', '-n', required=True,
              help='Database name')
def emu_restrict(source, taxids, output, name):
    """
    Create restricted EMU database for specific taxa.

    \b
    Example:
      mmonitor database emu restrict -s /data/emu_full -t 1279,1280,1281 -o /data/staph -n staph_db
    """
    taxid_list = [int(t.strip()) for t in taxids.split(',')]

    result = db_manager.restrict_emu_database(
        Path(source), taxid_list, Path(output), name
    )

    if result:
        click.echo(click.style(f"✓ Restricted database created at: {result}", fg='green'))
    else:
        click.echo(click.style("✗ Restriction failed", fg='red'), err=True)
        sys.exit(1)


@database_emu.command('build-from-sources')
@click.option('--name', '-n', required=True, help='Database name (used in manifest)')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output directory')
@click.option('--rrndb-version', default='latest',
              type=click.Choice(['latest', '5.8']),
              help='rrnDB version to use')
@click.option('--refseq-domains', default='bacteria,archaea',
              help='RefSeq domains to include (comma-separated: bacteria,archaea)')
@click.option('--email', default='mmonitor@example.com',
              help='Email for NCBI Entrez queries (required by NCBI)')
@click.option('--threads', '-t', type=int, default=4,
              help='Threads for indexing')
@click.option('--keep-downloads', is_flag=True,
              help='Keep downloaded source files after build')
def emu_build_from_sources(name, output, rrndb_version, refseq_domains,
                           email, threads, keep_downloads):
    """
    Build EMU database from rrnDB + NCBI RefSeq sources.

    This follows the original EMU database construction methodology:
    - Downloads rrnDB 16S sequences
    - Downloads NCBI 16S RefSeq for bacteria/archaea
    - Assigns NCBI taxonomy to all sequences
    - Removes duplicates (identical sequence + same species)
    - Builds minimap2 index
    - Creates FAIR-compliant manifest for reproducibility

    \b
    The resulting database can be shared and rebuilt by others using
    the manifest.json file which contains the exact build command.

    \b
    Example:
      mmonitor database emu build-from-sources -n emu_2024 -o /data/emu_2024
      mmonitor database emu build-from-sources -n emu_2024 -o /data/emu_2024 --rrndb-version 5.8

    \b
    Note: This process downloads several GB of data and may take
    several hours depending on your internet connection.
    """
    from mmonitor.cli.emu_db_builder import EmuDatabaseBuilder

    domains = [d.strip() for d in refseq_domains.split(',')]

    click.echo(f"Building EMU database '{name}'")
    click.echo(f"  rrnDB version: {rrndb_version}")
    click.echo(f"  RefSeq domains: {', '.join(domains)}")
    click.echo(f"  Output: {output}")
    click.echo()
    click.echo("This will download several GB of data and may take hours.")

    if not click.confirm("Continue?"):
        click.echo("Aborted.")
        return

    builder = EmuDatabaseBuilder(
        output_dir=Path(output),
        name=name,
        email=email,
        threads=threads
    )

    try:
        result = builder.build_from_sources(
            rrndb_version=rrndb_version,
            refseq_domains=domains,
            keep_downloads=keep_downloads
        )

        click.echo()
        click.echo(click.style("✓ Database built successfully!", fg='green'))
        click.echo()
        click.echo("Database statistics:")
        click.echo(f"  Total sequences: {builder.manifest.total_sequences}")
        click.echo(f"  Unique species: {builder.manifest.unique_species}")
        click.echo(f"  Duplicates removed: {builder.manifest.duplicates_removed}")
        click.echo()
        click.echo("To use this database:")
        click.echo(f"  mmonitor run taxonomy-16s -i reads.fastq -s sample -p project --emu-db {output}")
        click.echo()
        click.echo("FAIR manifest saved to:")
        click.echo(f"  {output}/manifest.json")
        click.echo()
        click.echo("Others can reproduce this database using the build command in the manifest.")

        # Update config
        cfg = config.config
        cfg.setdefault('databases', {})['emu'] = str(result)
        config.save(cfg)

    except Exception as e:
        click.echo(click.style(f"✗ Build failed: {e}", fg='red'), err=True)
        logger.exception("Build failed")
        sys.exit(1)


@database_emu.command('rebuild')
@click.option('--manifest', '-m', required=True, type=click.Path(exists=True),
              help='Path to manifest.json from a previous build')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output directory for rebuilt database')
@click.option('--threads', '-t', type=int, default=4,
              help='Threads for indexing')
def emu_rebuild(manifest, output, threads):
    """
    Rebuild EMU database from a FAIR manifest.

    This allows exact reproduction of a database build using
    the manifest.json file from a previous build.

    \b
    Example:
      mmonitor database emu rebuild -m /shared/emu_db/manifest.json -o /data/my_emu_db
    """
    from mmonitor.cli.emu_db_builder import EmuDatabaseBuilder, DatabaseBuildManifest

    click.echo(f"Loading manifest: {manifest}")

    manifest_data = DatabaseBuildManifest.load(Path(manifest))

    click.echo(f"Rebuilding database: {manifest_data.name}")
    click.echo(f"  Original build date: {manifest_data.build_date}")
    click.echo(f"  Original command: {manifest_data.build_command}")
    click.echo()

    params = manifest_data.parameters

    builder = EmuDatabaseBuilder(
        output_dir=Path(output),
        name=manifest_data.name,
        threads=threads
    )

    try:
        builder.build_from_sources(
            rrndb_version=params.get('rrndb_version', 'latest'),
            refseq_domains=params.get('refseq_domains', ['bacteria', 'archaea'])
        )

        click.echo(click.style("✓ Database rebuilt successfully!", fg='green'))

    except Exception as e:
        click.echo(click.style(f"✗ Rebuild failed: {e}", fg='red'), err=True)
        sys.exit(1)


@database_emu.command('info')
@click.argument('database_path', type=click.Path(exists=True))
def emu_info(database_path):
    """
    Show information about an EMU database.

    Displays the FAIR manifest including build parameters,
    sources, and checksums for reproducibility.

    \b
    Example:
      mmonitor database emu info /data/emu_db
      mmonitor database emu info /data/emu_db/manifest.json
    """
    from mmonitor.cli.emu_db_builder import DatabaseBuildManifest

    db_path = Path(database_path)

    # Find manifest
    if db_path.is_file() and db_path.name == 'manifest.json':
        manifest_path = db_path
    else:
        manifest_path = db_path / 'manifest.json'

    if not manifest_path.exists():
        click.echo(click.style("No manifest.json found in database directory", fg='yellow'))
        click.echo("This database may have been built with an older version of MMonitor")

        # Show basic file info
        fasta = db_path / 'species_taxid.fasta'
        taxonomy = db_path / 'taxonomy.tsv'

        if fasta.exists():
            size_mb = fasta.stat().st_size / (1024 * 1024)
            click.echo(f"\nFASTA file: {fasta.name} ({size_mb:.1f} MB)")

        if taxonomy.exists():
            import pandas as pd
            df = pd.read_csv(taxonomy, sep='\t')
            click.echo(f"Taxonomy file: {len(df)} species")

        return

    manifest = DatabaseBuildManifest.load(manifest_path)

    click.echo()
    click.echo("=" * 60)
    click.echo(f"EMU Database: {manifest.name}")
    click.echo("=" * 60)
    click.echo()
    click.echo(f"Version:      {manifest.version}")
    click.echo(f"Build date:   {manifest.build_date}")
    click.echo(f"Builder:      {manifest.builder_version}")
    click.echo()
    click.echo("Statistics:")
    click.echo(f"  Sequences:          {manifest.total_sequences:,}")
    click.echo(f"  Unique species:     {manifest.unique_species:,}")
    click.echo(f"  Duplicates removed: {manifest.duplicates_removed:,}")
    click.echo()

    if manifest.sources:
        click.echo("Sources:")
        for source in manifest.sources:
            click.echo(f"  - {source.get('name', 'Unknown')}")
            if 'version' in source:
                click.echo(f"    Version: {source['version']}")
            if 'url' in source:
                click.echo(f"    URL: {source['url']}")
            if 'download_date' in source:
                click.echo(f"    Downloaded: {source['download_date'][:10]}")
        click.echo()

    if manifest.build_command:
        click.echo("Reproducibility:")
        click.echo(f"  Build command: {manifest.build_command}")
        click.echo()

    click.echo("Checksums:")
    if manifest.fasta_checksum:
        click.echo(f"  FASTA:    {manifest.fasta_checksum[:16]}...")
    if manifest.taxonomy_checksum:
        click.echo(f"  Taxonomy: {manifest.taxonomy_checksum[:16]}...")
    if manifest.index_checksum:
        click.echo(f"  Index:    {manifest.index_checksum[:16]}...")

    if manifest.doi:
        click.echo()
        click.echo(f"DOI: {manifest.doi}")

    click.echo()


@database_emu.command('add-sequences')
@click.option('--database', '-d', required=True, type=click.Path(exists=True),
              help='Existing EMU database directory')
@click.option('--fasta', '-f', required=True, type=click.Path(exists=True),
              help='FASTA file with new sequences to add')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output directory for new database')
@click.option('--name', '-n', help='New database name (default: append to original)')
@click.option('--threads', '-t', type=int, default=4,
              help='Threads for indexing')
def emu_add_sequences(database, fasta, output, name, threads):
    """
    Add new sequences to an existing EMU database.

    Creates a new database combining the original sequences
    with additional sequences from the provided FASTA file.

    \b
    Example:
      mmonitor database emu add-sequences -d /data/emu_db -f new_seqs.fasta -o /data/emu_db_v2
    """
    from Bio import SeqIO

    db_path = Path(database)
    new_fasta = Path(fasta)
    out_path = Path(output)

    # Load existing database
    existing_fasta = db_path / 'species_taxid.fasta'
    existing_taxonomy = db_path / 'taxonomy.tsv'

    if not existing_fasta.exists():
        click.echo(click.style("Error: species_taxid.fasta not found in database", fg='red'))
        sys.exit(1)

    out_path.mkdir(parents=True, exist_ok=True)

    # Copy existing sequences
    click.echo("Copying existing sequences...")
    existing_count = 0
    with open(out_path / 'species_taxid.fasta', 'w') as out_f:
        for record in SeqIO.parse(existing_fasta, 'fasta'):
            SeqIO.write(record, out_f, 'fasta')
            existing_count += 1

        # Add new sequences
        click.echo("Adding new sequences...")
        new_count = 0
        for record in SeqIO.parse(new_fasta, 'fasta'):
            # Write in EMU format
            # Assuming header format: >taxid sequence_info
            out_f.write(f">{record.id}\n{str(record.seq)}\n")
            new_count += 1

    click.echo(f"Combined {existing_count} + {new_count} = {existing_count + new_count} sequences")

    # Copy taxonomy
    if existing_taxonomy.exists():
        shutil.copy(existing_taxonomy, out_path / 'taxonomy.tsv')
        click.echo("Copied taxonomy file")

    # Rebuild index
    click.echo("Building minimap2 index...")
    try:
        subprocess.run([
            'minimap2', '-d',
            str(out_path / 'species_taxid.fasta.mmi'),
            str(out_path / 'species_taxid.fasta'),
            '-t', str(threads)
        ], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        click.echo(click.style(f"Index build failed: {e.stderr.decode()}", fg='red'))
        sys.exit(1)

    click.echo(click.style(f"✓ Extended database created at: {out_path}", fg='green'))


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
        cfg.setdefault('databases', {})['centrifuger'] = str(result)
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
        cfg.setdefault('databases', {})['gtdbtk'] = str(result)
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
        cfg.setdefault('databases', {})['checkm2'] = str(result)
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
        cfg.setdefault('databases', {})['bakta'] = str(result)
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
    import json
    click.echo(json.dumps(config.config, indent=2))


@config_group.command('set')
@click.argument('key')
@click.argument('value')
def config_set(key, value):
    """
    Set a configuration value.

    \b
    Example:
      mmonitor config set threads 8
      mmonitor config set databases.emu /path/to/emu/db
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

    # Find server directory
    server_dir = PROJECT_ROOT / 'server'
    if not server_dir.exists():
        # Try parent (in case PROJECT_ROOT is desktop/)
        server_dir = PROJECT_ROOT.parent / 'server'
    if not server_dir.exists():
        click.echo(click.style(f"Error: Server directory not found", fg='red'), err=True)
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
