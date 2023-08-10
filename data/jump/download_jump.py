import pandas as pd
import scanpy as sc
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument(
    "-out", type="str", help="path to the output h5ad file", required=True
)
parser.add_argument(
    "-meta_dir",
    type="str",
    default="batch_correction_cellpainting/data/jump/",
    help="path to the directory with the metadata files",
    required=False,
)
args = parser.parse_args()


# read metadata files
plates = pd.read_csv(f"{args.meta_dir}plate.csv.gz")
wells = pd.read_csv(f"{args.meta_dir}well.csv.gz")
compound = pd.read_csv(f"{args.meta_dir}compound.csv.gz")
compound_meta = pd.read_table(f"{args.meta_dir}JUMP-Target-2_compound_metadata.tsv")
micro_filter = pd.read_csv(f"{args.meta_dir}microscope_filter.csv")
micro_config = pd.read_csv(f"{args.meta_dir}microscope_config.csv")

# helper formatter functions
profile_formatter = (
    "s3://cellpainting-gallery/cpg0016-jump/"
    "{Metadata_Source}/workspace/profiles/"
    "{Metadata_Batch}/{Metadata_Plate}/{Metadata_Plate}.parquet"
)


# Download TARGET2 plates
print("Start downloading TARGET2")
sample = plates.query("Metadata_PlateType=='TARGET2'")

dframes = []
columns = None
for _, row in sample.iterrows():
    s3_path = profile_formatter.format(**row.to_dict())
    dframes.append(
        pd.read_parquet(s3_path, storage_options={"anon": True}, columns=columns)
    )
dframes1 = pd.concat(dframes)


# Download COMPOUND plates (controls only)
print("Start downloading COMPOUND")
sample = plates.query("Metadata_PlateType=='COMPOUND'")

# the ids stem from the DMSO (JCP2022_033924) and POSCON8 plates
controls = wells[
    (wells.Metadata_JCP2022 == "JCP2022_085227")
    | (wells.Metadata_JCP2022 == "JCP2022_037716")
    | (wells.Metadata_JCP2022 == "JCP2022_025848")
    | (wells.Metadata_JCP2022 == "JCP2022_046054")
    | (wells.Metadata_JCP2022 == "JCP2022_035095")
    | (wells.Metadata_JCP2022 == "JCP2022_064022")
    | (wells.Metadata_JCP2022 == "JCP2022_050797")
    | (wells.Metadata_JCP2022 == "JCP2022_012818")
    | (wells.Metadata_JCP2022 == "JCP2022_033924")
]
controls = controls.drop("Metadata_JCP2022", axis=1)

dframes = []
columns = None
for _, row in sample.iterrows():
    s3_path = profile_formatter.format(**row.to_dict())
    temp = pd.read_parquet(s3_path, storage_options={"anon": True}, columns=columns)
    temp = temp.merge(
        controls, on=["Metadata_Source", "Metadata_Plate", "Metadata_Well"]
    )
    dframes.append(temp)
dframes2 = pd.concat(dframes)


# Combine Results
dframes1.insert(loc=0, column="Metadata_PlateType", value="TARGET2")
dframes2.insert(loc=0, column="Metadata_PlateType", value="COMPOUND")
dframes = pd.concat([dframes1, dframes2])


# Add metadata
metadata = compound.merge(wells, on="Metadata_JCP2022")
ann_dframe = metadata.merge(
    dframes, on=["Metadata_Source", "Metadata_Plate", "Metadata_Well"]
)


# Create and save AnnData object
adata = sc.AnnData(ann_dframe.iloc[:, 7:].to_numpy())
adata.var_names = ann_dframe.columns[7:]
adata.obs = ann_dframe.iloc[:, 0:7]


# Label the wells
adata.obs["Metadata_WellType"] = "treatment"
adata.obs.loc[
    adata.obs.Metadata_JCP2022 == "JCP2022_033924", "Metadata_WellType"
] = "DMSO"
adata.obs.loc[
    (adata.obs.Metadata_JCP2022 == "JCP2022_085227")
    | (adata.obs.Metadata_JCP2022 == "JCP2022_037716")
    | (adata.obs.Metadata_JCP2022 == "JCP2022_025848")
    | (adata.obs.Metadata_JCP2022 == "JCP2022_046054")
    | (adata.obs.Metadata_JCP2022 == "JCP2022_035095")
    | (adata.obs.Metadata_JCP2022 == "JCP2022_064022")
    | (adata.obs.Metadata_JCP2022 == "JCP2022_050797")
    | (adata.obs.Metadata_JCP2022 == "JCP2022_012818"),
    "Metadata_WellType",
] = "poscon"


# Add microscope and compound metadata
micro = micro_config.merge(micro_filter, on="Metadata_Filter_Configuration")
micro.Metadata_Source = "source_" + micro.Metadata_Source.astype(str)
adata.obs = adata.obs.merge(micro, on="Metadata_Source")

adata.obs = pd.merge(
    adata.obs,
    compound_meta,
    left_on="Metadata_InChIKey",
    right_on="InChIKey",
    how="left",
)


# Save anndata object
adata.write(args.out)
