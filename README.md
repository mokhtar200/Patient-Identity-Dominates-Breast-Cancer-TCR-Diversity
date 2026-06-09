# TCR beta-Chain Repertoire Analysis — Breast Cancer (PRJNA301507)

Power-law modelling and diversity analysis of T-cell receptor beta (TRB) repertoires from breast cancer patients. Compares Treg, CD4+ Memory, and CD45RA+ (Other) CD4+ T-cell subsets across blood, tumor tissue, and lymph node compartments.

## Data

- **Source**: Wang et al. (2017) *Cancer Immunol Res* — NCBI BioProject [PRJNA301507](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA301507)
- **Samples**: 202 TCRβ clonotype tables (CSV), 157 passed QC
- **Patients**: 33 breast cancer patients
- **Subsets**: Treg (71), CD4_Memory (38), Other/CD45RA+ (48)
- **Compartments**: Blood (94), Tumor (48), Lymph node (15)
- **Data files**: `data/SRR*.csv` with columns: AASeq, cloneCount, cloneFraction, Vregion, Jregion, Chain

## Directory Structure

```
Patient-Identity-Dominates-Breast-Cancer-TCR-Diversity/
├── README.md
├── scripts/
│   ├── 01_build_metadata.py       
│   ├── 02_main_pipeline.py        
│   └── 03_additional_analyses.py  
├── results/                        
├── figures/                       

```

## Analyses Performed

- Data loading & QC — Standardisation (CDR3 length 5–25aa, TRBV/J genes, ≥100 clonotypes, ≥1000 reads)
- Alpha diversity — Shannon, invSimpson, Gini, Pielou, Clonality, Hill numbers, Chao1, top-10% fraction
- Beta diversity — Bray-Curtis, PCoA, PERMANOVA, Morisita-Horn overlap
- Power-law modelling — Pooled & patient-level alpha exponents, Vuong model selection, x-min sensitivity
- Clonal-sharing network — Jaccard inter-patient network, Louvain communities, Kruskal-Wallis test
- Mixed-effects model — ICC estimation with cluster-bootstrap CI
- Centrality vs alpha — Spearman correlation
- V/J gene segment usage frequencies
- CDR3 length distribution & amino acid composition
- Differential clonal abundance (log2 FC, FDR-corrected)
- Compartment-specific diversity & Treg dominance
- PERMDISP & ANOSIM beta-diversity tests
- Clonal dominance profiling (top1/5/10/20 fractions)
- Public vs private clonotype sharing
- Random Forest classification (Treg vs CD4_Memory)
- Variance partitioning (mixed models)
- CDR3 sequence similarity network (Levenshtein distance 1)
- Within-patient subset overlap (Jaccard)
- Jensen-Shannon divergence between subsets

## Requirements

- Python 3.11+
- numpy, pandas, scipy, matplotlib, seaborn
- powerlaw, python-igraph, leidenalg
- statsmodels, scikit-learn

Install: `pip install -r requirements.txt`

## Usage

```bash
# 1. Build clean metadata
python scripts/01_build_metadata.py

# 2. Run primary analysis
python scripts/02_main_pipeline.py

# 3. Run additional analyses
python scripts/03_additional_analyses.py
```

## Key Results

| Metric | Treg | CD4_Memory | Other |
|--------|------|------------|-------|
| Pooled α | 2.214 | 2.293 | 2.340 |
| Patient α (median) | 1.55 | 1.41 | 1.64 |
| ICC(Shannon) | — | 0.831 [0.640–0.909] | — |
| PERMANOVA | — | R²=0.013, p=0.76 | — |
| Network density | — | 1.0 (33 nodes, 5 communities) | — |
| Public CDR3s | — | 35,362 / 743,702 (4.8%) | — |
| JSD | vs Mem: 0.74 | vs Other: 0.74 | vs Treg: 0.63 |
