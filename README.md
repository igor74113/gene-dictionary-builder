# Gene Dictionary Builder

This repository builds versioned gene dictionary JSON artifacts from authoritative public genomics sources.

A gene dictionary is a deterministic mapping from a gene identifier to a stable integer index.  
These dictionaries are intended for use in:
- sparse vector representations (e.g. COO / CSR)
- reproducible comparisons across samples and pipelines
- downstream annotation, weighting, and explanation layers

---

## What This Repository Produces

Running the builder generates two complementary JSON files, each with a distinct role.

### 1) gene_dictionary_pc_v1.json  
Canonical gene index (GENCODE protein-coding)

- Source: GENCODE human GTF
- Scope: protein-coding genes only
- Gene ID system: Ensembl Gene IDs (ENSG…)
- Ordering: deterministic (sorted Ensembl IDs)
- Purpose: defines the coordinate system for gene-level feature vectors

This dictionary answers the question:

“What genes exist, and which column index does each gene occupy?”

This file is treated as frozen by default.  
Once created, it will not be regenerated on subsequent runs unless explicitly forced.

---

### 2) gene_dictionary_clinvar_v1.json  
ClinVar-derived evidence / clinical-attention gene set

- Source: ClinVar tab-delimited exports  
  (prefers gene_summary.txt.gz, falls back to variant_summary.txt.gz)
- Gene ID system: NCBI GeneID (Entrez)
- Scope: genes that appear in ClinVar submissions
- Ordering: deterministic (sorted numeric GeneID)
- Purpose: evidence overlay, prioritization, and interpretability

This dictionary answers the question:

“Which genes have been observed or discussed in clinical variant submissions?”

This file is expected to change over time and is regenerated on every run.

---

## How the Dictionaries Are Intended to Be Used

| Use case | Protein-coding dictionary | ClinVar dictionary |
|--------|---------------------------|--------------------|
| Define sparse vector dimensions | Yes | No |
| Ensure vectors are comparable | Yes | No |
| Gene-level feature indexing | Yes | No |
| Weighting / prioritization | Optional | Yes |
| Annotation & explanation | Yes (decode indices) | Yes (clinical context) |

Only the GENCODE-based dictionary should be used to define vector indices.  
The ClinVar-derived dictionary is an overlay and must not be used as a coordinate system.

---

## Update Behavior

### Default behavior
When the script is run:

- If gene_dictionary_pc_v1.json already exists, it is not rebuilt
- gene_dictionary_clinvar_v1.json is always rebuilt

This prevents accidental changes to the gene index space while allowing regular updates of clinical evidence.

### Fresh directory
If neither file exists, both dictionaries are generated.

### Optional override
If enabled via CLI flag (if present in the script):

```bash
python3 build_gene_dictionaries.py --rebuild-gencode
