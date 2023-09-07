import scipy
from scipy.stats import median_abs_deviation

scipy.stats.median_absolute_deviation = median_abs_deviation
import anndata as ad
import pandas as pd
import numpy as np
import pycytominer as pc

metadata = pd.read_csv("../../data/jump/target2_metadata.txt", header=None)[0].tolist()
feature_names = pd.read_csv("../../data/jump/target2_feature_names.txt", header=None)[
    0
].tolist()
profiles_normalized = pd.read_csv(
    "../../data/jump/target2_profiles_normalized.csv", compression="gzip", index_col=0
)

# run feature selection
print("Start feature selection")
pc.feature_select(
    profiles_normalized,
    features=feature_names,
    operation=[
        "variance_threshold",
        "correlation_threshold",
        "drop_outliers",
        "noise_removal",
    ],
    noise_removal_perturb_groups="Metadata_JCP2022",
    output_file="../../data/jump/target2_profiles_featureselected.csv",
    compression_options="gzip",
)
print("Finished feature selection")
