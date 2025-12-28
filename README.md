# ![MMonitor](https://www.mmonitor.org/static/images/mmonitor_banner.png)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub release](https://img.shields.io/github/v/release/lucast122/MMonitor?color=brightgreen)](https://github.com/lucast122/MMonitor/releases)
[![Issues](https://img.shields.io/github/issues/lucast122/MMonitor)](https://github.com/lucast122/MMonitor/issues)

**MMonitor** is an open-source platform for **real-time metagenome monitoring** from Oxford Nanopore data.  
It pairs a **desktop app/CLI** for running pipelines with a **web dashboard** for interactive exploration.




## 📄 Publication

**MMonitor for Real-Time Monitoring of Microbial Communities Using Long Reads**
*Cell Reports Methods* (Elsevier, 2025)

🔗 **DOI:** [https://doi.org/10.1016/j.crmeth.2025.101266](https://doi.org/10.1016/j.crmeth.2025.101266)

**MMonitor** is an open-source software platform for real-time analysis and visualization of metagenomic Oxford Nanopore sequencing data. It supports both GUI- and CLI-based workflows and provides an interactive web dashboard for monitoring taxonomic composition, quality metrics, diversity indices, and metadata correlations over time.

📌 **Please cite this paper if you use MMonitor in your research.**


---

## 🔎 Table of Contents
- [ Features](#-features)
- [ Architecture](#-architecture)
- [ Quick Start](#-quick-start)
  - [macOS (App bundle)](#macos-app-bundle)
  - [Linux (from source)](#linux-from-source)
  - [Windows via WSL (CLI)](#windows-via-wsl-cli)
- [ Pipelines & Tools](#-pipelines--tools)
- [ Databases](#-databases)
- [ Configuration](#-configuration)
- [ Dashboard](#-dashboard)
- [ Troubleshooting](#-troubleshooting)
- [ Contributing](#-contributing)
- [ Citation](#-citation)
- [ License](#-license)

---

##  Features
-  **Real-time analysis** as reads are produced (watch a directory; incremental updates)
-  **GUI & CLI**: run locally on laptops/workstations or headless servers
-  **Configurable pipelines** (taxonomy, assembly/annotation; custom DBs)
-  **Interactive dashboard**: taxonomy, QC, diversity, functional summaries
-  **Automated reporting**: export figures/tables for manuscripts
-  **Cross-platform**: macOS, Linux; Windows via WSL (GUI for Windows coming soon)

---

##  Architecture
- **Desktop/CLI** (`desktop/`): orchestrates pipelines and file watching  
- **Web server** (`server/`): Django + Plotly Dash dashboards  
- **Pipelines**: call best-of-breed tools (e.g., Minimap2, Centrifuger, EMU, Flye, Bakta, CheckM2)  
- **Plug-and-play design**: packaged binaries/environments to avoid user setup where possible

---

##  Quick Start

### macOS (App bundle)
1. Download the latest **macOS release** from **[Releases](https://github.com/lucast122/MMonitor/releases)**.  
2. Unzip and open `MMonitor.app`. (First run may take longer to set up the internal tool environment.)  
3. Point the app to your **reads folder** and **database folder(s)**. Start monitoring.

> If Gatekeeper blocks the app, right-click → **Open** once to trust it.

---

### Linux (from source)
```bash
# 1) Get code
git clone https://github.com/lucast122/MMonitor.git
cd MMonitor

# 2) Create env (conda or micromamba recommended)
conda create -n mmonitor python=3.11 -y
conda activate mmonitor


# 3) Install desktop requirements
pip install -r desktop/requirements.txt

# 4) Run desktop (GUI or CLI)
export PYTHONPATH=$PYTHONPATH:$(pwd)/desktop:$(pwd)/desktop/src
python desktop/src/mmonitor/__main__.py  # starts GUI
# or: python desktop/src/mmonitor/cli.py --help
```
### Windows via WSL (CLI)

GUI for Windows is in development. For now, run the CLI inside WSL2.

Install WSL (Windows 10/11): open PowerShell as Admin:

wsl --install -d Ubuntu


Reboot if prompted and complete Ubuntu’s first-run setup.

Inside Ubuntu (WSL2):
```bash
sudo apt update && sudo apt -y install git wget
# Install micromamba (or use conda)
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C ~/bin/ --strip-components=1 bin/micromamba
echo 'eval "$(/home/$USER/bin/micromamba shell hook -s bash)"' >> ~/.bashrc
source ~/.bashrc

git clone https://github.com/lucast122/MMonitor.git
cd MMonitor
micromamba create -y -n mmonitor python=3.11
micromamba activate mmonitor
pip install -r desktop/requirements.txt
export PYTHONPATH=$PYTHONPATH:$(pwd)/desktop:$(pwd)/desktop/src
```

Run CLI (access Windows files under /mnt/c/Users/<you>/...):
```bash
python desktop/src/mmonitor/cli.py \
  run --reads /mnt/c/Users/<you>/data/fastq \
  --db /mnt/c/Users/<you>/databases/emu_db \
  --out /mnt/c/Users/<you>/mmonitor_results
```


⚠️ Performance tip: for many small files, use directories inside the WSL filesystem (e.g., ~/data) rather than /mnt/c.
