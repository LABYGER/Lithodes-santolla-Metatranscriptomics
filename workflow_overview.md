# Overview of Analytical Workflow
**Project Title:** Environmental Variability Drives Functional Plasticity in the Gill-Associated Microbiome of Lithodes santolla

This document describes the methodological workflow used for bioinformatic analyses and figure generation presented in this study.

---

## Tools and Versions

| Tool             | Version   | Purpose                           |
|------------------|-----------|-----------------------------------|
| **FastQC**       | `0.12.1`  | Raw reads quality check            |
| **Trimmomatic**  | `0.39`    | Adapter and quality trimming       |
| **SortMeRNA**    | `4.3.6`   | rRNA depletion                     |
| **Trinity**      | `2.15.1`  | De novo transcriptome assembly     |
| **Bowtie2**      | `2.5.4`   | Read mapping                       |
| **CD-HIT**       | `4.8.1`   | Redundancy removal                 |
| **BUSCO**        | `5.8.2`   | Assembly completeness assessment   |
| **TransDecoder** | `5.7.1`   | ORF prediction                     |
| **DIAMOND**      | `2021.11` | Protein homology search            |
| **MEGAN6 CE**    | `6.25.10` | Taxonomic binning                  |
| **eggNOG-mapper**| `2.1.12`  | Functional annotation              |
| **KAAS (KEGG)**  | Web       | KEGG pathway annotation            |
| **Salmon**       | `1.10.1`  | Transcript quantification           |
| **DESeq2 (R)**   | `1.40.0`  | Differential expression analysis    |
| **GOSeq (R)**    | `1.54.0`  | GO enrichment analysis              |

---

## Data Pre-processing

### Quality Control:
- Tool: **FastQC v0.11.9**
- Command: `fastqc *.fastq.gz`

### Adapter Trimming:
- Tool: **Trimmomatic v0.39**
- Command example: `ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:150`

### Ribosomal RNA Removal:
- Tool: **SortMeRNA v4.3.6**
- Databases: SILVA, Rfam

---

##  Transcriptome Assembly and Assessment

### De novo Assembly:
- Tool: **Trinity v2.15.1**
- Command: `Trinity --seqType fq --max_memory 100G --left ... --right ...`

### Assembly Quality:
- Mapping: **Bowtie2 v2.5.4**
- Completeness: **BUSCO v5.8.2** (arthropoda_odb10, eukaryota_odb10, bacteria_odb12, archaea_odb12)
- Redundancy Reduction: **CD-HIT v4.8.1** (95% identity)

---

##  Annotation and Taxonomic Assignment

### ORF Prediction:
- Tool: **TransDecoder v5.7.1**
- Parameters: minimum ORF 50 amino acids

### Sequence Similarity Search:
- Tool: **DIAMOND v2.1.8**
- Database: NCBI nr (downloaded 2024-11-02)
- Mode: more-sensitive, e-value < 1E-5

### Taxonomic Assignment:
- Tool: **MEGAN6 CE v6.25.10**, LCA method

### Functional Annotation:
- Tool: **eggNOG-mapper v2**, **KAAS (KEGG Automatic Annotation Server)**

---

##  Differential Expression Analysis

### Abundance Estimation:
- Tool: **Salmon within Trinity**
- Command: `align_and_estimate_abundance.pl`

### Normalization:
- Method: TMM (Trimmed Mean of M-values)

### Statistical Testing:
- Tool: **DESeq2 v1.40.2** in R (v4.4.1)
- Criteria: FDR < 0.05, |log2FC| > 2

---

##  Functional Enrichment

### GO Enrichment:
- Tool: **GOseq v1.50.0** in Trinity
- Inputs: DEG lists with GO annotations

### KEGG Pathways:
- Tool: **clusterProfiler v4.8.3** in R
- Visualization: ggplot2, enrichplot, pathview

---

## Visualization

- Complete list of scripts in `scripts/`
