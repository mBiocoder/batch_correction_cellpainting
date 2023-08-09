# JUMP Dataset Overview

The JUMP dataset includes a large production dataset with three perturbation modalities: chemical compounds, overexpression of genes using ORFs, and gene knockout using CRISPR guides (cpg0016), as well as three pilot datasets (cpg0000, cpg0001, cpg0002). The data is generated at Broad and partner sites.

For our purpose, **we only use the primary JUMP dataset (cpg0016)**. 
  - generated in 12 centers/sources (source_7 and source_13 are the same center)
      - genetic perturbation data was generated at specific sites (CRISPR at S13, ORFs at S4)
  - 116k chemical and > 15k genetic perturbations (over-expression of 12,602 genes, and knockout of 7,975 genes using CRISPR-Cas9)
  - Human U2OS cells
  - Notebook demonstrating how the data can be downloaded: https://github.com/jump-cellpainting/datasets/blob/main/sample_notebook.ipynb


### Overall Plate Layout

Unique JUMP identifiers were assigned to all compounds, ORF, and CRISPR perturbations, starting with "JCP2022_."

1) **Compound Plates**:
    * 384-well plates designed with controls in outer columns and treatment wells in inner columns to minimize well-position effects.
    * **Four replicates of eight compound positive controls in the outermost columns, followed by DMSO in the next-inner columns.**
    * Remaining wells contain one replicate of compound treatments. 
    * 1536-well plates follow the same layout as 384-well plates.

2) ORF Plates:
    * Pre-designed with negative control and untreated wells spread across the plate.
    * Positive controls added to specific wells.
    * Remaining wells contain ORF treatments, with one replicate per plate map and five replicate plates produced.

3) CRISPR Plates:
    * Outer columns have positive and negative controls.
    * Eight compound positive controls in the outermost columns, followed by negative controls and DMSO in other columns.

**In each batch of the cpg0016, control plates were run**: 

4) JUMP-Target-1-Compound Plates:
   - **Run with every batch of data generation**
        - 306 diverse compounds
            - 46 are included as positive controls
            - Negative control: DMSO (64 DMSO wells spread across the plate.)

5) **JUMP-Target-2-Compound Plates**:
    * Same set of compounds, but different Broad sample IDs and different layouts. JUMP-Target-2-Compound is used in the production of the JUMP dataset.

Plates run by **some** data acquisition partners periodically (e.g., one per batch, or at the beginning and end of each batch): 

6) Negative control plates
   * An entire plate of untreated cells
   * An entire plate of cells treated with non-targeting guides (for CRISPR)
   * An entire plate where all wells were treated with DMSO 

### Metadata

Schema: https://github.com/jump-cellpainting/datasets/blob/main/metadata/README.md

