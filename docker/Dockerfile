FROM broadinstitute/gatk:4.4.0.0

LABEL maintainer="variant-pipeline"
LABEL description="GATK4 somatic+germline variant calling + annotation"

ENV DEBIAN_FRONTEND=noninteractive

# ── System dependencies ───────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    bwa \
    tabix \
    bcftools \
    && rm -rf /var/lib/apt/lists/*

# ── Python dependencies ───────────────────────────────────────────────────────
# The GATK base image ships Python 3; we add bioinformatics + agent packages
RUN pip install --no-cache-dir \
    anthropic>=0.25.0 \
    cyvcf2>=0.30.22 \
    pandas>=2.0.3 \
    numpy>=1.24.4 \
    matplotlib>=3.7.2 \
    seaborn>=0.12.2 \
    requests>=2.31.0 \
    pyyaml>=6.0.1 \
    tqdm>=4.66.0 \
    plotly>=5.18.0 \
    python-dotenv>=1.0.0

# ── Copy pipeline code ────────────────────────────────────────────────────────
WORKDIR /pipeline
COPY annotation/ annotation/
COPY agent/ agent/
COPY config/ config/

# ── Verify tools ──────────────────────────────────────────────────────────────
RUN gatk --version && \
    samtools --version | head -1 && \
    python -c "import cyvcf2, anthropic, pandas; print('Python deps OK')"

ENTRYPOINT ["/bin/bash"]
