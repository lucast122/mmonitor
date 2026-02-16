# Database Auto-Download Rules
# These rules automatically download databases when no path is configured.
# Default location: ~/mmonitor_databases/<tool>/
#
# Each rule:
#   1. Checks if a database path is already configured (non-empty)
#   2. If empty, downloads to ~/mmonitor_databases/<tool>/
#   3. Writes the resolved path to a sentinel file for downstream rules

import os
from pathlib import Path

# DB_BASE is defined in the main Snakefile from config["db_base_dir"]


def get_db_path(tool, subkey="database"):
    """Get configured database path, or default download location."""
    configured = config.get(tool, {}).get(subkey, "")
    if configured:
        return configured
    return str(DB_BASE / tool)


# ============ EMU Database ============
rule setup_emu_db:
    """Download EMU default database if not configured"""
    output:
        sentinel=DB_BASE / "emu" / ".db_ready"
    params:
        db_dir=str(DB_BASE / "emu"),
        configured=config.get("emu", {}).get("database", "")
    conda:
        "../envs/taxonomy.yaml"
    log:
        str(DB_BASE / "emu" / "download.log")
    shell:
        """
        mkdir -p {params.db_dir}

        if [ -n "{params.configured}" ] && [ -d "{params.configured}" ]; then
            echo "Using configured EMU database: {params.configured}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        # Check if already downloaded
        if [ -f "{params.db_dir}/species_taxid.fasta" ]; then
            echo "EMU database already exists at {params.db_dir}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        echo "Downloading EMU default database..." | tee {log}
        # EMU uses osf.io for its default database
        cd {params.db_dir}
        emu download-db 2>&1 | tee -a {log} || {{
            echo "emu download-db failed, trying direct download..." | tee -a {log}
            wget -q "https://osf.io/56uf7/download" -O emu_database.tar.gz 2>&1 | tee -a {log}
            tar xzf emu_database.tar.gz 2>&1 | tee -a {log}
            rm -f emu_database.tar.gz
        }}

        touch {output.sentinel}
        """


# ============ CheckM2 Database ============
rule setup_checkm2_db:
    """Download CheckM2 database using built-in download command"""
    output:
        sentinel=DB_BASE / "checkm2" / ".db_ready"
    params:
        db_dir=str(DB_BASE / "checkm2"),
        configured=config.get("checkm2", {}).get("database", "")
    conda:
        "../envs/checkm2.yaml"
    log:
        str(DB_BASE / "checkm2" / "download.log")
    shell:
        """
        mkdir -p {params.db_dir}

        if [ -n "{params.configured}" ] && [ -d "{params.configured}" ]; then
            echo "Using configured CheckM2 database: {params.configured}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        # Check if already downloaded
        if [ -f "{params.db_dir}/uniref100.KO.1.dmnd" ]; then
            echo "CheckM2 database already exists at {params.db_dir}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        echo "Downloading CheckM2 database..." | tee {log}
        checkm2 database --download --path {params.db_dir} 2>&1 | tee -a {log}

        touch {output.sentinel}
        """


# ============ GTDB-TK Database ============
rule setup_gtdbtk_db:
    """Download GTDB-TK database using built-in download command"""
    output:
        sentinel=DB_BASE / "gtdbtk" / ".db_ready"
    params:
        db_dir=str(DB_BASE / "gtdbtk"),
        configured=config.get("gtdbtk", {}).get("database", "")
    conda:
        "../envs/gtdbtk.yaml"
    log:
        str(DB_BASE / "gtdbtk" / "download.log")
    shell:
        """
        mkdir -p {params.db_dir}

        if [ -n "{params.configured}" ] && [ -d "{params.configured}" ]; then
            echo "Using configured GTDB-TK database: {params.configured}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        # Check if already downloaded
        if [ -d "{params.db_dir}/taxonomy" ]; then
            echo "GTDB-TK database already exists at {params.db_dir}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        echo "Downloading GTDB-TK database (this may take a while, ~85GB)..." | tee {log}
        cd {params.db_dir}
        wget -q "https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz" \
            -O gtdbtk_data.tar.gz 2>&1 | tee -a {log}
        tar xzf gtdbtk_data.tar.gz --strip-components=1 2>&1 | tee -a {log}
        rm -f gtdbtk_data.tar.gz

        touch {output.sentinel}
        """


# ============ Bakta Database ============
rule setup_bakta_db:
    """Download Bakta database using built-in download command"""
    output:
        sentinel=DB_BASE / "bakta" / ".db_ready"
    params:
        db_dir=str(DB_BASE / "bakta"),
        db_type=config.get("bakta", {}).get("db_type", "light"),
        configured=config.get("bakta", {}).get("database", "")
    conda:
        "../envs/annotation.yaml"
    log:
        str(DB_BASE / "bakta" / "download.log")
    shell:
        """
        mkdir -p {params.db_dir}

        if [ -n "{params.configured}" ] && [ -d "{params.configured}" ]; then
            echo "Using configured Bakta database: {params.configured}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        # Check if already downloaded
        if [ -f "{params.db_dir}/db/version.json" ] || [ -f "{params.db_dir}/version.json" ]; then
            echo "Bakta database already exists at {params.db_dir}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        echo "Downloading Bakta {params.db_type} database..." | tee {log}
        bakta_db download --output {params.db_dir} --type {params.db_type} 2>&1 | tee -a {log}

        touch {output.sentinel}
        """


# ============ eggNOG Database ============
rule setup_eggnog_db:
    """Download eggNOG-mapper database using built-in download command"""
    output:
        sentinel=DB_BASE / "eggnog" / ".db_ready"
    params:
        db_dir=str(DB_BASE / "eggnog"),
        configured=config.get("eggnog", {}).get("database", "")
    conda:
        "../envs/functional.yaml"
    log:
        str(DB_BASE / "eggnog" / "download.log")
    shell:
        """
        mkdir -p {params.db_dir}

        if [ -n "{params.configured}" ] && [ -d "{params.configured}" ]; then
            echo "Using configured eggNOG database: {params.configured}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        # Check if already downloaded
        if [ -f "{params.db_dir}/eggnog.db" ]; then
            echo "eggNOG database already exists at {params.db_dir}" | tee {log}
            touch {output.sentinel}
            exit 0
        fi

        echo "Downloading eggNOG database..." | tee {log}
        download_eggnog_data.py --data_dir {params.db_dir} -y 2>&1 | tee -a {log}

        touch {output.sentinel}
        """


# ============ Database Setup Targets ============
rule setup_databases_16s:
    """Download all databases needed for 16S taxonomy analysis"""
    input:
        DB_BASE / "emu" / ".db_ready"

rule setup_databases_wgs:
    """Download all databases needed for WGS taxonomy analysis"""
    input:
        # Centrifuger requires manual database build (too large to auto-download)
        []

rule setup_databases_assembly:
    """Download all databases needed for the full assembly pipeline"""
    input:
        DB_BASE / "checkm2" / ".db_ready",
        DB_BASE / "gtdbtk" / ".db_ready",
        DB_BASE / "bakta" / ".db_ready",
        DB_BASE / "eggnog" / ".db_ready"

rule setup_all_databases:
    """Download all databases"""
    input:
        DB_BASE / "emu" / ".db_ready",
        DB_BASE / "checkm2" / ".db_ready",
        DB_BASE / "gtdbtk" / ".db_ready",
        DB_BASE / "bakta" / ".db_ready",
        DB_BASE / "eggnog" / ".db_ready"
