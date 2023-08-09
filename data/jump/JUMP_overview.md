# JUMP Dataset Overview

The JUMP dataset includes a large production dataset with three perturbation modalities: chemical compounds, overexpression of genes using ORFs, and gene knockout using CRISPR guides (cpg0016), as well as three pilot datasets (cpg0000, cpg0001, cpg0002). The data is generated at Broad and partner sites.

For our purposes we only focus on the compounds from the primary JUMP dataset and ignore the pilot datasets. The primary dataset (cpg0016) contains multiple subsets:
Compound dataset.
Open Reading Frame (ORF) overexpression dataset.
CRISPR knockout dataset.

The genetic perturbation data was generated at specific sites (CRISPR at S13, ORFs at S4), while compound replicates were generated at multiple sites and there are overlaps in compounds or genes tested in different datasets. The compound dataset contains images and profiles from the Cell Painting assay for over 116,750 unique compounds, over-expression of 12,602 genes, and knockout of 7,975 genes using CRISPR-Cas9. The data is collected from human osteosarcoma cells (U2OS) and is estimated to be 115 TB in size, capturing 1.6 billion cells and their single-cell profiles. We are only interested in the compound dataset, however.

The compound data was generated across 11 data sources. Each source, apart from source_7, exchanged their nominated compounds with either other sources where they were assayed using different instruments and microscopes. JUMP-Target-2-Compound, a positive control plate of 306 diverse compounds, was run with every batch of data generation. 

Control plates include negative control plates, JUMP-Target-2-Compound plates, and JUMP-Target-1-Compound plates (positive controls). Controls are integrated to identify and correct for experimental artifacts:
* Within-plate controls include negative controls, positive controls, and untreated wells.
* Negative control plates encompass untreated cells, non-targeting CRISPR guides, and DMSO-treated cells.
* Positive control wells consist of distinct compound signatures and protein overexpression controls.

Download: https://github.com/jump-cellpainting/datasets/blob/main/sample_notebook.ipynb


### Overall Plate Layout

Unique JUMP identifiers was assigned to all compound, ORF, and CRISPR perturbations, starting with "JCP2022_."

1) Compound Plates:
* 384-well plates designed with controls in outer columns and treatment wells in inner columns to minimize well-position effects.
* Four replicates of eight compound positive controls in outermost columns, followed by DMSO in next-inner columns.
* Remaining wells contain one replicate of compound treatments. 
* 1536-well plates follow the same layout as 384-well plates.

2) ORF Plates:
* Pre-designed with negative control and untreated wells spread across the plate.
* Positive controls added to specific wells.
* Remaining wells contain ORF treatments, with one replicate per plate map and five replicate plates produced.

3) CRISPR Plates:
* Outer columns have positive and negative controls.
* Eight compound positive controls in outermost columns, followed by negative controls and DMSO in other columns.

4) JUMP-Target-1-Compound Plates:
* Contain 306 compounds and DMSO, mostly in singlicates.
* 64 DMSO wells spread across the plate.
* Some compounds have two replicates.

5) JUMP-Target-2-Compound Plates:
* Same content as JUMP-Target-1-Compound Plates, different layout. 




