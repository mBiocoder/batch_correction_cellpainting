# Recursion RxRx19b Dataset Overview

RxRx19b is the second component of the RxRx19 dataset series released by Recursion sharing data from a high-dimensional human cellular assay for COVID-19 associated disease. 

In severe cases, COVID-19 is known to ignite a chronic state of inflammation, culminating in acute respiratory distress syndrome (ARDS). This high-dimensional cellular model of severe, inflammatory COVID-19 was treated with a library of approved drugs in service of finding compounds that can relieve the damaging, hyper-inflammatory state.

RxRx19b modeles the cytokine storm in endothelial cells associated with late-stage COVID-19 using cocktails of circulating proteins. Cocktails were prepared to mimic concentrations of circulating soluble factors in healthy and COVID-19-infected patients.

Vascular endothelial cells (HUVEC) were treated with each cocktail and resulting morphological changes were observed.
High-dimensional phenotype resulting from the morphological changes was screened with 1,856 FDA-approved drugs and tool compounds.
Thus, this dataset represents inflammatory effects and potential treatments in the context of COVID-19 ARDS.


Download: https://www.rxrx.ai/rxrx19b

The file `embeddings.csv` contains 70,384 embeddings of 128 features each.

| Level                 | Count             |
|-----------------------|-------------------|
| cell type             | HUVEC             |
| experiment number     | only HUVEC-1      |
| plates                | 53                |
| well locations        | up to 1380        |
| site                  | 1                 |
| channel               | 3                 |

**Note:** The wells in the dataset may have different sizes or cell numbers associated with them. A higher count for a particular well indicates that it has more observations (cells) compared to wells with lower counts.
It's important to note that the count alone doesn't provide information about the actual size or number of cells in each well. It simply reflects the number of occurrences of each unique well value in this dataset.

RxRx19b consists of 70,384 fluorescence microscopy images and their deep learning embeddings. Each image dimension is 2048x2048x6. The evaluated perturbation consists of 1,856 small molecules at 4-6 concentrations in three COVID-19-associated cytokine storm conditions.

Specifically, the conditions are
* Severe cytokine storm: This condition involves exposing the cells to a cocktail of circulating proteins that mimic the cytokine storm observed in severe COVID-19 patients.
* Healthy: This condition involves exposing the cells to a cocktail of circulating proteins that mimic the levels observed in healthy individuals.
* No cytokines: This condition involves exposing the cells to a control condition with no added cytokines.

In total we have these many disease_conditions:
* storm-severe:    67364
* healthy:         2756
* blank:            264

In a single plate (except plate 53) we have:
* storm-severe:    1276
* healthy:         52
* blank:           0

Plate 53 has the same condition of "blank" 264 times, 52 times "healthy" and 1012 times "storm-severe".

## Metadata

| Column             | Description                                                                     |
|------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| site_id            | Unique identifier of a given site                                               |
| well_id            | Unique identifier of a given well                                               |
| cell_type          | Cell type tested, namely HUVEC                                                  |
| experiment         | Experiment identifier                                                           |
| plate              | Plate number within the experiment                                              |
| well               | Location on the plate                                                           |
| site               | Location in the well of the image taken                                         |
| disease_condition  | The disease condition tested in the well (healthy, healthy cytokine cocktail; storm-severe, severe cytokine storm cocktail; or blank, no cytokines)                                                |
| treatment          | Compound tested in the well (if any)                                            |
| treatment_conc     | Compound concentration tested (in uM)                                           |
| SMILES             | Formula of tested compound (as CXSMILES/ChemAxon Extended SMILES)               |
