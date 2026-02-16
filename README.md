# ![MMonitor](https://www.mmonitor.org/static/images/mmonitor_banner.png)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub release](https://img.shields.io/github/v/release/lucast122/MMonitor?color=brightgreen)](https://github.com/lucast122/MMonitor/releases)
[![Issues](https://img.shields.io/github/issues/lucast122/MMonitor)](https://github.com/lucast122/MMonitor/issues)
[![Journal](https://img.shields.io/badge/Cell%20Reports%20Methods-2025-green)](https://doi.org/10.1016/j.crmeth.2025.101266)

**MMonitor** is an open-source platform for **real-time metagenome monitoring** from Oxford Nanopore data.
It pairs a **desktop app / CLI** for running pipelines with a **web dashboard** for interactive exploration.

## Publication

**MMonitor for Real-Time Monitoring of Microbial Communities Using Long Reads**
*Cell Reports Methods* (Elsevier, 2025)

**DOI:** [https://doi.org/10.1016/j.crmeth.2025.101266](https://doi.org/10.1016/j.crmeth.2025.101266)

Please cite this paper if you use MMonitor in your research.

---

## Table of Contents
- [Features](#features)
- [Architecture](#architecture)
- [Installation](#installation)
- [Usage](#usage)
  - [GUI](#gui)
  - [CLI](#cli)
  - [Web Dashboard](#web-dashboard)
- [Pipelines](#pipelines)
- [Databases](#databases)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

---

## Features
- **Real-time analysis** as reads are produced (directory watching with incremental updates)
- **GUI & CLI** — run locally on laptops/workstations or headless servers
- **Configurable pipelines** — taxonomy (16S & WGS), assembly, binning, annotation, functional analysis
- **Interactive dashboard** — taxonomy, QC, diversity, functional summaries, MAG browser
- **Cross-platform** — Linux, macOS, and Windows via WSL

---

## Architecture
- **Desktop/CLI** (`desktop/`) — orchestrates Snakemake pipelines and file watching
- **Web server** (`server/`) — Django backend with a React frontend dashboard
- **Pipelines** — Snakemake workflows calling best-of-breed tools (EMU, Centrifuger, Flye, Medaka, MetaBAT2, CheckM2, GTDB-TK, Bakta, eggNOG-mapper, InterProScan)
- **pip-installable** — all pipeline definitions and resources ship inside the Python package

> **Note:** Pre-built app bundles (macOS `.app`, etc.) are no longer provided.
> The Snakemake + conda environment architecture makes standalone binaries impractical.
> Instead, MMonitor is now installed via `pip` as described below.

---

## Installation

### Requirements
- Python 3.10 or 3.11 (3.11 recommended)
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/) (for Snakemake conda environments)

### Linux & macOS

```bash
# 1. Create and activate a conda environment
conda create -n mmonitor python=3.11 -y
conda activate mmonitor

# 2. Clone the repository
git clone https://github.com/lucast122/MMonitor.git
cd MMonitor

# 3. Install MMonitor (with GUI support)
pip install -e "./desktop[gui]"

# Or without GUI (CLI/server only):
pip install -e "./desktop"
```

### Windows (via WSL2)

MMonitor runs on Windows through WSL2 (Windows Subsystem for Linux).

**1. Install WSL** — open PowerShell as Administrator:
```powershell
wsl --install -d Ubuntu
```
Reboot if prompted and complete Ubuntu's first-run setup.

**2. Inside Ubuntu (WSL2):**
```bash
sudo apt update && sudo apt -y install git wget

# Install Miniforge (conda + mamba)
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
# restart shell, then:

conda create -n mmonitor python=3.11 -y
conda activate mmonitor

git clone https://github.com/lucast122/MMonitor.git
cd MMonitor
pip install -e "./desktop[gui]"
```

> **Performance tip:** Store data inside the WSL filesystem (e.g. `~/data`) rather than under `/mnt/c/` for much better I/O performance.

### Verify installation

```bash
mmonitor --version
# mmonitor, version 0.3.0
```

---

## Usage

### GUI

Launch the graphical interface:

```bash
mmonitor gui
```

The GUI lets you:
- Configure server connections and credentials
- Select input directories and set pipeline parameters
- Run analysis pipelines with real-time progress
- Monitor directories for new reads (real-time mode)
- Download and manage databases

> **Note:** On headless servers (no display), use the CLI or web dashboard instead.
> On WSL2, you need WSLg (Windows 11) or an X server for GUI support.

### CLI

MMonitor provides a full command-line interface for all operations:

```bash
# Run 16S taxonomy analysis
mmonitor run taxonomy-16s \
  -i reads.fastq -s sample1 -p myproject \
  --emu-db /path/to/emu_db

# Run WGS taxonomy analysis
mmonitor run taxonomy-wgs \
  -i reads.fastq -s sample1 -p myproject \
  --centrifuger-db /path/to/cfr_db

# Run full assembly pipeline (assembly + binning + annotation + functional)
mmonitor run assembly-full \
  -i reads.fastq -s sample1 -p myproject

# Dry run (show what would be executed)
mmonitor run taxonomy-16s \
  -i reads.fastq -s sample1 -p myproject \
  --emu-db /path/to/emu_db --dry-run

# See all available commands
mmonitor --help
mmonitor run --help
```

### Web Dashboard

Start the local web server and open the dashboard in your browser:

```bash
mmonitor serve
# Dashboard available at http://127.0.0.1:8000/dashboard/

# Custom host/port:
mmonitor serve --host 0.0.0.0 --port 9000

# Don't open browser automatically:
mmonitor serve --no-browser
```

---

## Pipelines

All pipelines are run via Snakemake with automatic conda environment management.

| Command | Pipeline | Tools |
|---------|----------|-------|
| `taxonomy-16s` | 16S rRNA taxonomy | Filtlong, EMU |
| `taxonomy-wgs` | Whole-genome taxonomy | Filtlong, Centrifuger |
| `assembly` | Metagenomic assembly | Flye, Medaka |
| `binning` | MAG binning + QC | MetaBAT2, CheckM2 |
| `functional` | Full functional analysis | All assembly + annotation tools |
| `assembly-full` | Complete WGS pipeline | Flye, Medaka, MetaBAT2, CheckM2, GTDB-TK, Bakta, eggNOG, InterProScan |

---

## Databases

MMonitor can download and manage databases for you:

```bash
# Download EMU database for 16S analysis
mmonitor database emu download --preset silva138

# Download Bakta database for gene annotation
mmonitor database bakta download --type light

# Download CheckM2 database for MAG quality assessment
mmonitor database checkm2 download

# List installed databases
mmonitor database list

# Verify database integrity
mmonitor database verify -t emu
```

Database paths are saved automatically to `~/.mmonitor/config.yaml`.

---

## Configuration

MMonitor stores its configuration in `~/.mmonitor/config.yaml`.

```bash
# Show current configuration
mmonitor config show

# Set a value
mmonitor config set threads 8
mmonitor config set emu.database /path/to/emu_db

# Initialize config with defaults
mmonitor config init
```

You can also export a configuration from the GUI and pass it to the CLI:

```bash
mmonitor run taxonomy-16s -i reads.fastq -s sample1 -p project1 \
  --config-file my_config.yaml
```

---

## Troubleshooting

**`mmonitor gui` fails with "no display"**
You are on a headless server or SSH session without X forwarding. Use `mmonitor serve` for the web dashboard, or connect with `ssh -X` for X11 forwarding.

**Snakemake conda environment errors**
Make sure `conda` or `mamba` is available in your PATH. Snakemake uses conda to create isolated environments for each tool.

**Pipelines fail to find databases**
Run `mmonitor database list` to check installed databases. Set paths with `mmonitor config set emu.database /path/to/db` or pass them directly via CLI flags.

**Import errors after installation**
Make sure you installed with `pip install -e "./desktop[gui]"` (with the `[gui]` extra for GUI support).

---

## Contributing

Contributions are welcome! Please open an issue or pull request on [GitHub](https://github.com/lucast122/MMonitor).

```bash
# Install with dev dependencies
pip install -e "./desktop[gui,dev]"

# Run tests
pytest desktop/tests/
```

---

## Citation

If you use MMonitor in your research, please cite:

> **MMonitor for Real-Time Monitoring of Microbial Communities Using Long Reads**
> *Cell Reports Methods* (Elsevier, 2025)
> DOI: [https://doi.org/10.1016/j.crmeth.2025.101266](https://doi.org/10.1016/j.crmeth.2025.101266)

---

## License

[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0)
