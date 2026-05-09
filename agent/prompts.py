"""
prompts.py — Structured prompt templates for variant interpretation

These templates are designed to give Claude enough structured context to
produce clinically grounded, reproducible variant interpretations.

Design principles:
  1. Structured input: all annotation fields in labelled sections
     (Claude performs better with labelled fields than raw text)
  2. Role priming: sets domain expectations before the question
  3. Output format specified: ensures consistent, parseable responses
  4. Temperature 0.2: low temperature for clinical reproducibility
  5. Tier assignment separated from free-text: enables downstream filtering

SOMATIC TIER DEFINITIONS:
  DRIVER     — variant has strong evidence as a cancer-causing driver
               (known oncogene hotspot OR loss-of-function in known TSG,
                COSMIC count > 50, functional evidence in literature)
  LIKELY_DRIVER — variant in a driver gene but not at a canonical hotspot,
                  or COSMIC count 10-50, or CIViC evidence B/C
  PASSENGER  — variant in a cancer gene but no evidence of driver function
  VUS_SOMATIC— variant of unknown somatic significance
  ARTIFACT   — technical artifact that passed filters (low VAF + artifact features)

GERMLINE TIER DEFINITIONS (ACMG/AMP 2015 framework):
  PATHOGENIC        — meets ≥ 2 PVS / 1 PVS + 1 PS / 2 PS / ... criteria
  LIKELY_PATHOGENIC — meets criteria below Pathogenic threshold
  VUS               — uncertain significance (conflicting or insufficient evidence)
  LIKELY_BENIGN     — meets BP criteria
  BENIGN            — meets BA1 (common in gnomAD) or ≥ 2 BP criteria

ACMG CRITERIA ENCODED IN PROMPTS:
  PVS1 — null variant in LOF gene
  PS1  — same amino acid change as known Pathogenic
  PM2  — absent from gnomAD (or extremely rare)
  PP3  — multiple in silico predictors support pathogenicity (CADD ≥ 25)
  BA1  — allele frequency > 5% in gnomAD
  BP4  — multiple in silico predictors support benign (CADD < 15)
  BP7  — synonymous, no splice effect predicted
"""

SYSTEM_PROMPT = """You are a clinical bioinformatician and molecular geneticist specializing
in variant interpretation for both somatic oncology and germline rare disease contexts.
You produce concise, accurate, evidence-based variant summaries following ACMG/AMP 2015
guidelines (germline) and ESMO/ASCO oncogenicity criteria (somatic).

Be factual. Do not speculate beyond the evidence provided. If evidence is insufficient,
state that clearly. Use standard clinical genetics terminology."""


SOMATIC_PROMPT_TEMPLATE = """
Interpret the following somatic variant from a tumor-normal sequencing assay.

## VARIANT COORDINATES
- Location: {chrom}:{pos} {ref}>{alt}
- HGVS genomic: {hgvsg}
- HGVS coding:  {hgvsc}
- HGVS protein: {hgvsp}

## GENE & TRANSCRIPT
- Gene symbol:  {gene_symbol} ({ensembl_gene_id})
- Gene biotype: {gene_biotype}
- Transcript:   {transcript_id} (canonical: {is_canonical}, MANE: {mane_select})

## CONSEQUENCE
- SO consequence: {consequence}
- IMPACT level:   {impact}
- Exon/Intron:    {exon} / {intron}
- Amino acid change: {amino_acids} (codons: {codons})

## IN SILICO PREDICTIONS
- SIFT:      {sift_score} ({sift_pred})
- PolyPhen:  {polyphen_score} ({polyphen_pred})
- CADD phred: {cadd_phred}
- Conservation: phyloP={phylop}, GERP={gerp}

## POPULATION FREQUENCY (gnomAD v4)
- Global AF:    {gnomad_af}
- Popmax AF:    {gnomad_popmax_af} ({gnomad_popmax_pop})
- AFR: {gnomad_afr_af} | AMR: {gnomad_amr_af} | EAS: {gnomad_eas_af}
- NFE: {gnomad_nfe_af} | SAS: {gnomad_sas_af} | ASJ: {gnomad_asj_af}

## CLINICAL DATABASES
- ClinVar significance: {clinvar_clnsig} (★{clinvar_review_stars})
- ClinVar disease:      {clinvar_clndn}
- Known variants (dbSNP/ClinVar ID): {existing_variation}

## CANCER DATABASES
- CIViC evidence level: {civic_evidence_level} ({civic_evidence_type})
- CIViC drug:           {civic_drug}
- CIViC cancer type:    {civic_cancer_type}

## TUMOR DATA
- Variant allele frequency (VAF): {vaf}
- Sequencing depth at site:       {depth}

## TASK
1. Provide a 3–5 sentence clinical interpretation of this somatic variant.
   Address: gene function, oncogenic mechanism (if applicable), VAF interpretation,
   population frequency context, and therapeutic implications (if any evidence exists).

2. Assign a tier: DRIVER | LIKELY_DRIVER | PASSENGER | VUS_SOMATIC | ARTIFACT

3. List the key evidence points supporting your tier assignment (bullet points).

Format your response as:
INTERPRETATION: <3-5 sentences>
TIER: <one of the options above>
EVIDENCE:
- <point 1>
- <point 2>
...
"""


GERMLINE_PROMPT_TEMPLATE = """
Interpret the following germline variant identified in a clinical sequencing context.

## VARIANT COORDINATES
- Location: {chrom}:{pos} {ref}>{alt}
- HGVS genomic: {hgvsg}
- HGVS coding:  {hgvsc}
- HGVS protein: {hgvsp}

## GENE & TRANSCRIPT
- Gene symbol:  {gene_symbol} ({ensembl_gene_id})
- Gene biotype: {gene_biotype}
- Transcript:   {transcript_id} (canonical: {is_canonical}, MANE: {mane_select})

## CONSEQUENCE
- SO consequence: {consequence}
- IMPACT level:   {impact}
- Exon/Intron:    {exon} / {intron}
- Amino acid change: {amino_acids} (codons: {codons})

## IN SILICO PREDICTIONS (ACMG PP3/BP4)
- SIFT:      {sift_score} ({sift_pred})
- PolyPhen:  {polyphen_score} ({polyphen_pred})
- CADD phred: {cadd_phred}   [≥25 supports PP3; <15 supports BP4]
- Conservation: phyloP={phylop}, GERP={gerp}

## POPULATION FREQUENCY — ACMG BA1 / PM2
- Global AF:    {gnomad_af}       [>5% → BA1 Benign; absent → PM2 Pathogenic]
- Popmax AF:    {gnomad_popmax_af} ({gnomad_popmax_pop})
- faf95:        {gnomad_popmax_af} [filtering AF at 95% CI — used for BA1/PM2]
- AFR: {gnomad_afr_af} | AMR: {gnomad_amr_af} | EAS: {gnomad_eas_af}
- NFE: {gnomad_nfe_af} | SAS: {gnomad_sas_af} | FIN: {gnomad_fin_af}

## CLINICAL SIGNIFICANCE — ACMG PS1/PP5
- ClinVar significance: {clinvar_clnsig} (★{clinvar_review_stars})
- ClinVar disease:      {clinvar_clndn}
- ClinVar last evaluated: {clinvar_last_evaluated}
- Submitter count:      {clinvar_submission_count}

## REGULATORY CONTEXT
- Regulatory feature:   {regulatory_id} ({regulatory_type})
- Known variant IDs:    {existing_variation}

## GENOTYPE (trio)
- Genotype: {genotype}   [0/1=het, 1/1=hom-alt, 0/0=ref]
- Depth:    {depth}

## TASK
1. Provide a 3–5 sentence clinical interpretation following ACMG/AMP 2015 guidelines.
   Address: gene disease association (if known), consequence type,
   population frequency, ClinVar evidence, and in silico support.

2. Assign an ACMG classification: PATHOGENIC | LIKELY_PATHOGENIC | VUS | LIKELY_BENIGN | BENIGN

3. List the ACMG criteria that apply (e.g., PVS1, PS1, PM2, PP3, BA1, BP4, BP7).

Format your response as:
INTERPRETATION: <3-5 sentences>
CLASSIFICATION: <ACMG class>
ACMG_CRITERIA: <comma-separated list of applicable criteria>
EVIDENCE:
- <point 1>
- <point 2>
...
"""


def build_somatic_prompt(variant: dict) -> str:
    """Fill the somatic prompt template with variant annotation fields."""
    from collections import defaultdict
    safe = defaultdict(lambda: ".", {k: v if v else "." for k, v in variant.items()})
    return SOMATIC_PROMPT_TEMPLATE.format_map(safe)


def build_germline_prompt(variant: dict) -> str:
    """Fill the germline prompt template with variant annotation fields."""
    from collections import defaultdict
    safe = defaultdict(lambda: ".", {k: v if v else "." for k, v in variant.items()})
    return GERMLINE_PROMPT_TEMPLATE.format_map(safe)
