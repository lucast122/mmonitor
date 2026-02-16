FROM condaforge/mambaforge:latest

LABEL maintainer="MMonitor Team"
LABEL description="MMonitor bioinformatics pipelines (Snakemake + conda)"

# Install Snakemake
RUN mamba install -y -c bioconda -c conda-forge \
    snakemake>=8.0 \
    && mamba clean -a -y

# Copy workflow files
COPY desktop/workflows/ /opt/mmonitor/workflows/

WORKDIR /opt/mmonitor/workflows

# Pre-create conda environments so they don't need to be built at runtime.
# This makes the Docker image larger but startup much faster.
RUN snakemake --use-conda --conda-create-envs-only --cores 1 \
    --configfile config/config.yaml \
    || true

# Default entry point: run Snakemake
ENTRYPOINT ["snakemake", "--snakefile", "/opt/mmonitor/workflows/Snakefile"]
CMD ["--help"]
