# variant-calling-workshop

A hands-on educational workshop covering GATK best-practice variant calling for both **somatic** (tumor-vs-normal) and **germline** (trio joint genotyping) contexts, with a multi-database annotation layer and an LLM-based interpretation agent (Claude API).

All data is from GATK's public tutorial bundles — no institutional data, no large downloads (~100 MB total).

**[→ View pipeline outputs and figures (DEMO.md)](DEMO.md)**

---

## What this pipeline does

```
GATK Tutorial BAMs (chr17 / chr20 subset, ~100 MB total)
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 1 — GERMLINE (HaplotypeCaller GVCF workflow)    │
│  MarkDuplicates → HC per-sample GVCF (trio)             │
│  → CombineGVCFs → GenotypeGVCFs → Hard filters          │
│  → VEP + gnomAD + ClinVar + CADD annotation             │
└─────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 2 — SOMATIC (Mutect2, tumor-normal mode)        │
│  Precomputed Mutect2 VCF (HCC1143 breast cancer)        │
│  → VEP + gnomAD + ClinVar + CADD + CIViC annotation     │
└─────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 3 — AGENTIC LAYER (Claude API)                  │
│  Annotated TSV → Structured prompt per variant          │
│  → ACMG/AMP tier assignment → Natural language summary  │
└─────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 4 — RESULTS NOTEBOOKS                           │
│  germline_results.ipynb + somatic_results.ipynb         │
│  → VAF, CADD, consequence, de novo, agent output        │
└─────────────────────────────────────────────────────────┘
```

---

## Somatic vs Germline — the conceptual distinction

| Aspect | Somatic | Germline |
|---|---|---|
| **Origin** | Acquired in a single cell lineage after conception | Present in every cell, inherited or de novo |
| **VAF** | Often < 50% (subclonal), can be very low (1–5%) | ~50% (het) or ~100% (hom) in diploid regions |
| **Comparator** | Tumor vs matched normal (same patient) | Joint-called trio (mother, father, son) |
| **Key artifact concern** | FFPE oxidation damage, strand bias, normal contamination | PCR errors, mapping errors near repeats |
| **Caller** | Mutect2 (tumor-normal mode) | HaplotypeCaller (GVCF workflow) |
| **Classification frame** | AMP/ASCO Tier I–IV (actionability) | ACMG/AMP 5-tier (pathogenicity) |

---

## Reference genome

Both modules use **GRCh37 (b37)** coordinates:
- Germline: chr20 subset (~65 Mb), CEPH trio NA12878/NA12891/NA12892
- Somatic: chr17:7M–20M (covers TP53, BRCA1), HCC1143 breast cancer cell line

All annotation APIs are queried against GRCh37:
- VEP: `https://grch37.rest.ensembl.org`
- gnomAD: dataset `gnomad_r2_1` (GRCh37, 125k exomes + 71k genomes)
- CADD: `GRCh37-v1.6` endpoint

---

## GATK best practices — what we follow and what we skip

### What we follow
- **MarkDuplicates**: removes PCR/optical duplicate reads before calling
- **Mutect2 tumor-normal mode**: matched normal subtracts germline variants from somatic callset
- **GetPileupSummaries + CalculateContamination**: quantifies cross-sample contamination
- **LearnReadOrientationModel**: detects FFPE oxidative damage artifacts (C→T orientation bias)
- **FilterMutectCalls**: full Mutect2 filter cascade using contamination + orientation model
- **HaplotypeCaller GVCF mode**: per-sample GVCFs enable joint genotyping without reprocessing
- **CombineGVCFs + GenotypeGVCFs**: joint genotyping across trio
- **Hard filters for germline**: correct alternative to VQSR on small callsets (< 30 samples)

### What we skip, and why
- **BQSR**: minimal benefit on modern sequencers (NovaSeq); requires 16–32 GB RAM on full genomes
- **VQSR**: requires thousands of variants to train; fails on small-interval callsets
- **Somatic BAMs**: HCC1143 tutorial BAMs are not publicly available; pipeline starts from precomputed filtered VCF

---

## Annotation layers

| Layer | Source | What it adds | API |
|---|---|---|---|
| 1 | **VEP** (Ensembl) | Consequence, gene, HGVS, SIFT, PolyPhen, exon/intron | GRCh37 REST |
| 2 | **gnomAD** v2.1 | Global AF + 8 ancestry populations | GraphQL (`gnomad_r2_1`) |
| 3 | **ClinVar** | Pathogenicity, review stars (1–4), condition | NCBI eUtils |
| 4 | **CADD** v1.6 | Deleteriousness phred score (SNVs only) | GET per position |
| 5 | **CIViC** | Cancer drug evidence (somatic only) | GraphQL |

**Checkpointing**: each layer saves a JSON checkpoint in `.annotation_cache/`. If the pipeline is interrupted, it resumes from the last completed layer — avoiding expensive API re-calls.

**API rate limits**: gnomAD enforces strict rate limiting (~1 req/2s). For 1111 germline variants this takes ~37 min. Production deployments should use local tabix lookups against downloaded gnomAD VCFs.

---

## Pipeline structure

```
variant-calling-workshop/
├── README.md
├── environment.yml                    ← conda environment (variant-pipeline)
├── .env                               ← API keys (gitignored)
├── .gitignore
├── config/
│   └── pipeline_config.yml            ← all paths and parameters
├── data/
│   ├── download_tutorial_data.sh      ← GATK GCS tutorial data download
│   ├── bams/germline/                 ← mother.bam, father.bam, son.bam
│   ├── ref/                           ← b37 chr20 reference FASTA
│   └── resources/                     ← dbSNP, gnomAD AF VCF, PON, intervals
├── germline/
│   ├── 01_preprocessing.sh            ← MarkDuplicates
│   ├── 02_haplotypecaller.sh          ← HC per-sample GVCF (trio)
│   ├── 03_genotyping.sh               ← CombineGVCFs + GenotypeGVCFs
│   ├── 04_hardfilter_annotate.sh      ← Hard filters + 5-layer annotation
│   └── notebooks/
│       └── germline_results.ipynb     ← depth, CADD, consequences, rare variants, agent output
├── somatic/
│   ├── 01_preprocessing.sh            ← MarkDuplicates (requires HCC1143 BAMs)
│   ├── 02_mutect2_calling.sh          ← Mutect2 + pileups + contamination
│   ├── 03_filtering.sh                ← LearnReadOrientationModel + FilterMutectCalls
│   ├── 04_annotation.sh               ← 5-layer annotation starting from precomputed VCF
│   └── notebooks/
│       └── somatic_results.ipynb      ← VAF, consequence, CIViC evidence, agent output
├── annotation/
│   ├── annotate_variants.py           ← orchestrates all 5 layers + checkpointing
│   ├── vep_annotator.py               ← VEP REST (GRCh37 region endpoint)
│   ├── gnomad_annotator.py            ← gnomAD v2.1 GraphQL (gnomad_r2_1 dataset)
│   ├── clinvar_annotator.py           ← ClinVar eUtils (esearch + esummary)
│   ├── cadd_annotator.py              ← CADD GET /api/v1.0/GRCh37-v1.6/{chrom}:{pos}
│   ├── civic_annotator.py             ← CIViC GraphQL (somatic drug evidence)
│   └── hgvs_formatter.py              ← HGVS notation normalizer
├── agent/
│   ├── variant_interpreter.py         ← Claude API agent (--model, --top, --all flags)
│   └── prompts.py                     ← ACMG germline + AMP somatic prompt templates
├── results/
│   ├── germline/annotated/            ← annotated.tsv, annotated.vcf.gz, interpreted.json
│   └── somatic/annotated/             ← annotated.tsv, annotated.vcf.gz, interpreted.json
└── logs/
    ├── germline/                      ← per-step logs
    └── somatic/                       ← per-step logs
```

---

## Agentic interpretation layer

`agent/variant_interpreter.py` sends each candidate variant to Claude with a structured prompt and parses a tiered clinical summary.

**Two-stage classification:**
1. **Rule-based pre-classification** — variants with BA1 (AF > 5%), clear MODIFIER impact, or no annotation are classified without an LLM call (saves API cost)
2. **LLM interpretation** — the top N unclassified variants (by priority score) are sent to Claude

**Germline tiers** (ACMG/AMP):
`PATHOGENIC` | `LIKELY_PATHOGENIC` | `VUS` | `LIKELY_BENIGN` | `BENIGN`

**Somatic tiers** (AMP/ASCO/CAP):
`TIER_I` | `TIER_II` | `VUS_SOMATIC` | `PASSENGER` | `UNCLASSIFIED`

**CLI flags:**
```bash
python agent/variant_interpreter.py \
  --tsv results/germline/annotated/annotated.tsv \
  --type germline \
  --top 50 \                          # number of variants to send to LLM
  --model claude-haiku-4-5-20251001 \ # default (cheapest); override with --model
  --output results/germline/annotated/interpreted.json
```

**Cost guide** (approximate, 50 variants):

| Model | Cost | Quality |
|---|---|---|
| `claude-haiku-4-5-20251001` | ~$0.01 | Good for structured ACMG tasks |
| `claude-sonnet-4-5-20250929` | ~$0.06 | Better reasoning |
| `claude-opus-4-5-20251101` | ~$0.30 | Best |

**API key**: store in `.env` at project root (gitignored):
```
ANTHROPIC_API_KEY=sk-ant-...
```

---

## Quick start

```bash
# 1. Create conda environment
conda env create -f environment.yml
conda activate variant-pipeline

# 2. Download tutorial data (~100 MB)
bash data/download_tutorial_data.sh

# 3. Run germline pipeline (chr20, CEPH trio)
bash germline/01_preprocessing.sh
bash germline/02_haplotypecaller.sh
bash germline/03_genotyping.sh
bash germline/04_hardfilter_annotate.sh

# 4. Run somatic annotation (starts from precomputed Mutect2 VCF)
bash somatic/04_annotation.sh

# 5. Run agentic interpretation
python agent/variant_interpreter.py \
    --tsv results/germline/annotated/annotated.tsv \
    --type germline --top 50 \
    --output results/germline/annotated/interpreted.json

python agent/variant_interpreter.py \
    --tsv results/somatic/annotated/annotated.tsv \
    --type somatic --top 26 \
    --output results/somatic/annotated/interpreted.json

# 6. Explore results
jupyter lab germline/notebooks/germline_results.ipynb
jupyter lab somatic/notebooks/somatic_results.ipynb
```

---

## Limitations and production equivalents

| This workshop | Production equivalent |
|---|---|
| BQSR skipped | Required for older platforms; 16–32 GB RAM |
| VQSR → hard filters | Required for WGS > 30 samples; needs HapMap/Omni |
| chr17/chr20 subsets only | Full genome: ~500 GB disk, 64 GB RAM |
| VEP via REST API | Local cache (~50 GB) for throughput |
| gnomAD via GraphQL API (slow) | Local tabix on downloaded VCF (instant) |
| CADD via web API (SNVs only) | Local CADD score files for indels too |
| No CNV/SV calling | GATK CNV, Manta, DRAGEN SV |
| No tumor purity estimation | PURPLE or similar |
| No phasing | Compound het detection requires phasing |
| Single LLM call per variant | Batch + cache for production throughput |

---

## Data sources

- **GATK tutorial bundle**: `gs://gatk-tutorials/workshop_2002/` — chr17 somatic + chr20 germline trio
- **Reference genome**: GRCh37 (b37) chr20 subset
- **gnomAD**: v2.1.1 (GRCh37) via GraphQL API
- **ClinVar**: via NCBI eUtils API
- **CADD**: v1.6 GRCh37 via web API
- **VEP**: Ensembl GRCh37 REST API
- **CIViC**: via GraphQL API

---

## License

MIT — educational use only, not for clinical decision-making.
