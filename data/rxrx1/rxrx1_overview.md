# Recursion rxrx1 Dataset Overview

Download: https://www.rxrx.ai/rxrx1

The file `embeddings.csv` contains 125,510 embeddings of 128 features each.

| Level                 | Count |
|-----------------------|-------|
| experiment = batch    | 51    |
| plates per experiment | 4     |
| wells per plate       | 308   |
| images per well       | 2     |

Cell lines/ cell types used are
* HEPG2 is a human hepatocellular carcinoma cell line that is commonly used to study liver metabolism and drug toxicity.
* HUVEC (human umbilical vein endothelial cells) is a cell line that is derived from the endothelium of the human umbilical vein. These cells are often used in vascular research and drug discovery.
* RPE (retinal pigment epithelial cells) is a cell line that is derived from the retinal pigment epithelium, which is a layer of cells in the eye. RPE cells are commonly used in ophthalmic research to study diseases such as age-related macular degeneration.
* U2OS is a human osteosarcoma cell line that is frequently used in cancer research, particularly in studies related to DNA damage response and repair.

One cell type is used per experiment/batch. In total, 4 cell types were used: HUVEC (24 times), RPE (11 times), HepG2 (11 times) and U2OS (5 times). 
This makes a total of 125,664 images, out of those, 154 were excluded during quality control, resulting in the final
number of 125,510 images/embeddings.

Generally, there are
* 1 untreated well (negative control)
* 30 wells treated with control siRNAs (positive control, 30 different siRNAs that are the same across each plate and
  each experiment)
* 277 wells treated with different siRNAs (treatment).

In each experiment, the same 1,108 different non-control siRNAs are used (277*4 = 1,108), meaning one non-control siRNA is introduced into one well per experiment. 
The locations of the non-control siRNAs are randomized, so that the well positions of replicates differ. It's not
explicitly stated on the website, but the positive control well positions also vary between plates.

**Note**: The negative control is always in well B02. Some plates have additional negative controls (e.g. HUVEC-24,
plate 1). In these cases, there are less positive control and/or treatment wells.

## Metadata

| Column     | Description                                                                             |
|------------|-----------------------------------------------------------------------------------------|
| site_id    | Unique ID of each image (composed of cell_type, experiment, plate, well and site)       |
| well_id    | Unique well ID (composed of cell_type, experiment, plate and well)                      |
| cell_type  | Name of the cell type tested (one of HEPG2, HUVEC, RPE, U2OS)                           |
| dataset    | `train` or `test`, depending on which split a site belongs to                           |
| experiment | Unique ID of the experiment (composed of cell_type and a number)                        |
| plate      | Plate number within experiment (one of 1, 2, 3, 4)                                      |
| well       | Well location on the plate                                                              |
| site       | Location in the well where image was taken (1 or 2)                                     |
| well_type  | `treatment`, `negative_control` or `positive_control`                                   |
| sirna      | ThermoFisher ID of the siRNA introduced to the well. `EMPTY` for negative control wells |
| sirna_id   | siRNAs mapped to integers                                                               |
