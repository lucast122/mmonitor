# 📊 MMonitor

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub release](https://img.shields.io/github/v/release/lucast122/MMonitor?color=brightgreen)](https://github.com/lucast122/MMonitor/releases)
[![Build](https://img.shields.io/github/actions/workflow/status/lucast122/MMonitor/ci.yml?branch=main)](https://github.com/lucast122/MMonitor/actions)
[![Issues](https://img.shields.io/github/issues/lucast122/MMonitor)](https://github.com/lucast122/MMonitor/issues)

**MMonitor** is an open-source platform for **real-time metagenome monitoring** using Oxford Nanopore sequencing.  
It combines a desktop app for running pipelines with a web dashboard for interactive visualization.

---

## ✨ Features
- ⚡ **Real-time analysis** of Nanopore data as it is generated  
- 🖥️ **GUI & CLI** for flexible usage  
- 🔧 **Configurable pipelines** with custom reference databases  
- 📊 **Interactive dashboard** for taxonomy, QC, diversity, and functional insights  
- 📝 **Automated reporting** with exportable results  
- 🖱️ **Cross-platform** (macOS, Linux; limited Windows support)  

---

## 📦 Installation

### Requirements
- Python 3.10+  
- Git  
- (Recommended) Conda for environment management  
- R packages: `jpeg`, `png`, `RColorBrewer`, `lattice`, `latticeExtra`  
- External tools (e.g. Minimap2, MetaFlye, Medaka, Centrifuger, Emu, Bakta, CheckM2)  

### macOS
1. Download the latest release from this repository.  
2. Unzip and launch `mmonitor`.  
3. First startup can take a while, as the tool runs setup in the background once.  

### Linux
```bash
git clone https://github.com/lucast122/MMonitor.git
cd MMonitor
conda create -n mmonitor python=3.11
conda activate mmonitor
pip install -r desktop/requirements.txt
export PYTHONPATH=$PYTHONPATH:MMonitor/desktop/:MMonitor/desktop/src/
python desktop/src/mmonitor/__main__.py
