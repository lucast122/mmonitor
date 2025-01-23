# 📊 MMonitor Documentation

## Table of Contents

1. [📖 Introduction](#-introduction)
2. [💾 Installation Instructions](#-installation-instructions)
3. [🚀 Quick Start Guide](#-quick-start-guide)
4. [🛠️ Detailed Usage Instructions](#%EF%B8%8F-detailed-usage-instructions)
5. [⚙️ Configuration and Settings](#%EF%B8%8F-configuration-and-settings)
6. [📈 Output and Results](#-output-and-results)
7. [🌐 Web Dashboard](#-web-dashboard)
8. [🐞 Error Handling and Troubleshooting](#-error-handling-and-troubleshooting)
9. [🤝 Contributing](#-contributing)
10. [📄 Licensing and Acknowledgments](#-licensing-and-acknowledgments)

## 📖 Introduction

**MMonitor** is an open-source software platform designed for the real-time analysis and visualization of metagenomic datasets produced by Oxford Nanopore Technologies sequencing. It provides a user-friendly interface for automating taxonomic and functional analysis of metagenomes, offering both a desktop application for running bioinformatics pipelines and a web-based dashboard for interactive result inspection.

### Key Features

- **⚡ Real-time Analysis**: Processes sequencing data as it is generated, enabling immediate insights.
- **🖥️ User-friendly Interface**: Offers both a graphical user interface (GUI) and a command-line interface (CLI) for flexibility.
- **🔧 Customizable Pipelines**: Allows users to configure analysis pipelines and use custom reference databases.
- **📊 Visualization Dashboard**: Provides dynamic insights into taxonomic composition over time, quality scores, diversity indices, and taxonomy-metadata correlations.
- **🖱️ Cross-platform Support**: Available for macOS, Linux, and Windows (limited functionality on Windows).

### Intended Audience

MMonitor is designed for researchers in biology related fields who need to analyze microbial communities using nanopore sequencing (taxonomy and function) in real-time.

---

## 💾 Installation Instructions

### Supported Operating Systems

- **🍎 macOS**
- **🐧 Linux (Unix)**
- **🪟 Windows** (Note: Limited functionality due to pipeline compatibility)

### Prerequisites

- **Hardware Requirements**:
  - 💻 Modern CPU
  - 🧮 16-32 GB of RAM
- **Software Requirements**:
  - 🐍 Python 3.10 or higher
  - 🛠️ Git (for cloning repositories)
  - 🐚 Conda (optional, recommended for environment management)
  - 📦 R and R packages: `jpeg`, `png`, `RColorBrewer`, `lattice`, `latticeExtra`
  - 🧰 External tools (included in builds or installed via requirements): Minimap2, MetaFlye, Medaka, Centrifuger, Emu, Bakta, CheckM2, etc.

### Installation Steps

#### 🍎 macOS

1. **Download the Prebuilt Application**:

   - Visit the [MMonitor website](https://mmonitor.org/) and download the latest version for macOS.

2. **Run the Application**:

   - Unzip the downloaded application.
   - Open the file named `mmonitor`.
   - If prompted with a security warning, you may need to adjust your security settings to allow the app to run.

3. **Verify Installation**:

   - Wait for the application to start.
   - If the application doesn't start, check the console output for error messages.

#### 🐧 Linux (Unix)

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/lucast122/MMonitor.git
   ```

2. **Set Up the Environment** (optional but recommended):

   ```bash
   conda create -n mmonitor python=3.11
   conda activate mmonitor
   ```

3. **Install Required Python Packages**:

   ```bash
   pip install -r MMonitor/desktop/requirements.txt
   ```

4. **Install External Dependencies**:

   ```bash
   conda install -c bioconda minimap2
   ```

5. **Set the PYTHONPATH**:

   ```bash
   export PYTHONPATH=$PYTHONPATH:MMonitor/desktop/
   export PYTHONPATH=$PYTHONPATH:MMonitor/desktop/src/
   ```

6. **Run the Application**:

   ```bash
   python MMonitor/desktop/src/mmonitor/__main__.py
   ```

   - If you encounter an error about missing modules, ensure that the `PYTHONPATH` is correctly set.

#### 🪟 Windows

**Note**: Running analysis pipelines is currently not supported on Windows due to compatibility issues with certain bioinformatics tools.

1. **Clone the Repository**:

   - Open PowerShell and run:

     ```powershell
     git clone https://github.com/lucast122/MMonitor.git
     ```

2. **Set Up the Environment** (optional):

   ```powershell
   conda create -n mmonitor python=3.10
   conda activate mmonitor
   ```

3. **Install Required Python Packages**:

   ```powershell
   pip install -r MMonitor\desktop\requirements.txt
   ```

4. **Run the Application**:

   ```powershell
   python MMonitor\desktop\src\mmonitor\__main__.py
   ```

   - Note: You can use MMonitor on Windows to add metadata and inspect local databases, but running pipelines requires a Unix-based system.

---

## 🚀 Quick Start Guide

### Running a Sample Analysis

1. **🔐 Login or Register**:

   - Visit the [MMonitor website](https://mmonitor.org/) and register an account or log in if you already have one.

2. **🖥️ Launch MMonitor**:

   - Run the MMonitor application as per your operating system instructions.

3. **👤 User Authentication**:

   - In the application, click on **User Authentication**.
   - Enter your username and password that you registered on the website.
   - Save the configuration.

4. **🧪 Process Sequencing Data**:

   - Click on **Process Sequencing Files**.
   - Select the desired analysis:
     - **🔬 Quick Taxonomy 16S Nanopore** for 16S rRNA reads.
     - **🧬 Quick Taxonomy Nanopore** for whole-genome sequencing (WGS) reads.
     - For functional analysis, select **🧩 Functional Analysis Pipeline** and choose the desired steps.

5. **➕ Add Samples**:

   - **Single Sample**:
     - Fill in the required sample information.
     - Click **Add FASTQ from folder** and select the folder containing your sequencing data.
   - **Multiple Samples**:
     - Prepare a CSV file with sample information (an example can be found in the `mmonitor_tutorial` folder).
     - Click **Add multiple samples from CSV** and select your CSV file.

6. **▶️ Run Analysis**:

   - Click **Submit** to start the analysis.
   - Wait for the analysis to complete (a notification will appear).

7. **🔎 View Results**:

   - Visit the [MMonitor Dashboard](https://mmonitor.org/dashboard) and log in.
   - Explore your results using the various visualization tools available.

---

## 🛠️ Detailed Usage Instructions

### User Interface Interaction

- **🖱️ Graphical User Interface (GUI)**:
  - The GUI provides an intuitive way to interact with MMonitor, suitable for users who prefer a visual approach.
- **💻 Command-Line Interface (CLI)**:
  - Advanced users can use the CLI for more control and to run analyses on remote servers.

### Important Commands and Options (CLI)

- **Main Analysis Types**:

  ```bash
  python MMonitor/desktop/src/mmonitor/__main__.py -a [analysis_type] -c [config_file] [options]
  ```

  - `analysis_type`: `taxonomy-wgs`, `taxonomy-16s`, `assembly`, `functional`, `stats`
  - `config_file`: Path to the JSON config file.

- **Options**:

  - `-i`, `--input`: Input files or directories.
  - `-s`, `--sample`: Sample name.
  - `-d`, `--date`: Sample date (YYYY-MM-DD).
  - `-p`, `--project`: Project name.
  - `-u`, `--subproject`: Subproject name.
  - `-m`, `--multicsv`: Path to CSV file for multiple samples.
  - `-b`, `--barcodes`: Use barcode column from CSV for multiplexing.
  - `--emu-db`: Path to custom Emu database.
  - `--centrifuge-db`: Path to Centrifuge database.
  - `-t`, `--threads`: Number of threads to use.
  - `--overwrite`: Overwrite existing records.
  - `-q`, `--qc`: Calculate QC statistics.
  - `-x`, `--update`: Update counts and abundances in the database.

### Customization

- **🗄️ Custom Databases**:
  - Users can specify custom databases for taxonomic profilers like Emu and Centrifuge.
  - Use the **Database Manager** in the GUI or provide paths via CLI options.

- **⚙️ Changing Pipelines**:
  - Currently, pipelines are predefined, but users can adjust parameters and select which steps to run.
  - Advanced users can modify pipeline scripts or contribute to the project for additional customization.

---

## ⚙️ Configuration and Settings

### Configuring the App

- **🔐 User Authentication**:
  - Configure your login credentials in the **User Authentication** window.
  - Credentials are stored in a JSON configuration file (`db_config.json`).

- **🔗 External Tools and Paths**:
  - Paths to external tools are managed automatically in the builds.
  - If running from source, ensure that tools like Minimap2 are installed and available in your system's PATH.

- **📝 Changing Parameters**:
  - Parameters for analyses can be adjusted in the GUI before running a pipeline.
  - In the CLI, parameters are adjusted using command-line options.

### Configuration Files

- **📁 Database Configuration**:
  - The `db_config.json` file stores database connection settings.
  - This file is created and managed by the GUI but can be edited manually if needed.

---

## 📈 Output and Results

### Types of Outputs

- **🦠 Taxonomic Profiles**:
  - Generated using Emu (for 16S rRNA) and Centrifuger (for WGS).
  - Results include species-level taxonomic abundances.

- **📊 Quality Control Reports**:
  - Read lengths, quality scores, GC content, and more.

- **🧬 Functional Analyses**:
  - Metagenome-assembled genomes (MAGs).
  - KEGG metabolic pathway maps.

### Accessing Results

- **💾 Local Storage**:
  - Results are stored in the `pipeline_out` directory within the application resources.
  - File formats include TSV, CSV, and JSON.

- **🌐 Web Dashboard**:
  - Results are uploaded to the MMonitor web server for visualization.
  - Access via the Dashboard after logging in.

---

## 🌐 Web Dashboard

### Interacting with the Dashboard

- **🌐 Access**:
  - Visit the [MMonitor Dashboard](https://mmonitor.org/dashboard) and log in with your credentials.

- **🗺️ Navigation**:
  - Use the menu to navigate between different visualization tools:
    - **🦠 Taxonomy**: View taxonomic profiles.
    - **📊 Diversity**: Explore alpha and beta diversity metrics.
    - **🔍 Quality Control**: Inspect sequencing quality metrics.
    - **🔗 Correlations**: Analyze correlations between taxa and metadata.
    - **🧬 KEGG**: View KEGG metabolic pathway maps.

### Visualizations and Metrics

- **📈 Plots and Graphs**:
  - Stacked bar charts, area plots, horizon graphs, heatmaps, PCoA plots, and more.

- **📊 Metrics**:
  - Alpha Diversity (Shannon, Simpson indices).
  - Beta Diversity (Bray-Curtis distances).

- **🖱️ Interactivity**:
  - Hover over plots for detailed information.
  - Filter and select samples based on criteria.

### Metadata Integration

- **📄 Uploading Metadata**:
  - In the **Correlations** section, upload a CSV file containing metadata.
  - Ensure that sample names in the metadata match those in your analysis.

- **🔗 Correlation Analysis**:
  - Analyze Pearson, Spearman, or Kendall correlations between taxa abundances and metadata variables.

---

## 🐞 Error Handling and Troubleshooting

### Common Issues

- **❌ Dependency Problems**:
  - Missing external tools or libraries.
  - Version incompatibilities.

- **⚠️ Runtime Errors**:
  - Insufficient memory or CPU resources.
  - Incorrect file paths or permissions.

- **🚫 Pipeline Failures**:
  - Errors during assembly, binning, or annotation steps.

### Diagnosing and Fixing Issues

- **📜 Check Logs**:
  - Review console output for error messages.
  - Log files may be generated in the output directories.

- **🔍 Verify Installations**:
  - Ensure all dependencies are correctly installed.
  - For missing R packages, run:

    ```R
    install.packages(c("jpeg", "png", "RColorBrewer", "lattice", "latticeExtra"))
    ```

- **🖥️ System Resources**:
  - Close unnecessary applications to free up memory.
  - Consider running on a system with more RAM or CPU cores.

### Getting Help

- **🐙 GitHub Issues**:
  - Report issues or seek assistance by opening an issue on the [MMonitor Pipeline GitHub repository](https://github.com/lucast122/mmonitor-pipeline/issues).

- **✉️ Contact Support**:
  - Email the developer at [timo-niklas.lucas@uni-tuebingen.de](mailto:timo-niklas.lucas@uni-tuebingen.de).

---

## 🤝 Contributing

### How to Contribute

- **💻 Code Contributions**:
  - Fork the repositories:
    - [MMonitor Pipeline](https://github.com/lucast122/mmonitor-pipeline)
    - [MMonitor Server](https://github.com/lucast122/mmonitor-server)
  - Create a new branch for your feature or bug fix.
  - Submit a pull request for review.

- **🐞 Bug Reports and Feature Requests**:
  - Use GitHub Issues to report bugs or suggest new features.

### Contribution Guidelines

- **📐 Coding Standards**:
  - Follow Python PEP 8 style guidelines.
  - Write clear, concise commit messages.

- **✅ Pull Request Process**:
  - Ensure that all tests pass before submitting.
  - Provide a detailed description of changes.

- **💬 Community Engagement**:
  - Participate in discussions and help answer questions.

---

## 📄 Licensing and Acknowledgments

### License

MMonitor is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

### Acknowledgments

- **👥 Contributors**:
  - **Timo N. Lucas**: Lead developer.
  - **Tobias Laas**: Development contributions.
  - **Simon Konzalla**: Development contributions.

- **💡 Funding**:
  - Funding by CMFI

- **🤝 Collaborations**:
  - MMonitor was developed at the University of Tübingen in collaboration with environmental biotechnologists and microbiome researchers from the Angenent Lab and Ley Lab. Special thanks to:
  - Daniel Huson, Largus Angenent, and Ruth Ley for supervising the project
  - Ulrike Biehain for helping in planning and testing the software, suggesting new features
  - Jeon Byoung Seon and Kurt Gemeinhardt for providing the metagenomic whole genome sequencing data
  - Yihua Liu, Soyoung Ham, for helping to test the application during development and providing feedback


### Conflicts of Interest

- The authors declare no conflicts of interest.

---

**Note**: For detailed information and methods, refer to the publication:
## Automation of Taxonomic and Functional Analysis for the Monitoring of Metagenomes in real-tiem using Nanopore Sequencing and MMonitor
## How to Automate Taxonomic and Functional Analysis of Metagenomes in real-time using Nanopore Sequencing and MMonitor
# Timo N. Lucas, Ulrike Biehain, Anupam Gautam, Tobias Laas, Simon Konzalla, Largus T. Angenent, Ruth E. Ley, Daniel H. Huson

