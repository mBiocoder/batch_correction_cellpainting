import scarches
import pandas as pd
import urllib.request
from bs4 import BeautifulSoup
from typing import Literal
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import seaborn as sns

alphabet = list("abcdefghijklmnopqrstuvwxyz".upper())


def get_scarches_version():
    """Return the version of the SCARCHES package."""

    return scarches.__version__


def get_sirna_target_gene(
    sirna: str, return_as: Literal["gene_symbol", "entrez"] = "gene_symbol"
) -> str:
    """Function to get the target gene(s) of a siRNA given its ThermoFisher ID

    Given the ThermoFisher ID of a siRNA, its target gene(s) are fetched from the
    ThermoFisher website.
    If a siRNA has multiple target genes, the returned string consists of all targets
    separated by an underscore "_".
    If no targets can be found for a siRNA, pd.NA is returned.

    Args:
        sirna:      ThermoFisher ID of a siRNA
        return_as:  determines whether the gene symbol (gene_symbol) or ENTREZ
                    gene ID (entrez) of the target(s) are returned

    Returns:
        str:    target gene(s) of the given siRNA, retrieved from ThermoFisher website

    """

    page = urllib.request.urlopen(
        "https://www.thermofisher.com/order/genome-database/AssayDetails?rev=V302&"
        f"assayid={sirna}&assaytype=sirna&CID=&ICID="
    )
    soup = BeautifulSoup(page, "html.parser")
    results = soup.find_all("dl")

    if len(results) == 0:
        # retrun NA if no target gene was found
        return pd.NA

    genes = []
    if return_as == "gene_symbol":
        for r in results:
            gene = str(r).split("<dt>Gene Symbol</dt>\n<dd>")[1].split("</dd>")[0]
            genes.append(gene)
    else:
        for r in results:
            gene = (
                str(r)
                .split("<dt>Entrez Gene ID</dt>\n<dd><")[1]
                .split(">")[1]
                .split("<")[0]
            )
            genes.append(gene)
    return "_".join(genes)


def letter_to_int(letter: str) -> int:
    """Get the index of a well row identified by a letter of the alphabet.

    This function is used by :func:`~bcc.utils.plot_well_type_positions`.

    Given the well row as letter, return the index such that A = 1, B = 2 and so on.
    For any input that is not a single uppercase letter of the alphabet, an Error is raised.

    Args:
        letter: uppercase letter identifying the well row

    Returns:
        int:    index of the well row

    """

    if not letter in alphabet:
        raise ValueError("The input must be one uppercase letter of the alphabet")
    return alphabet.index(letter) + 1


def plot_well_type_positions(
    metadata_df: pd.DataFrame,
    well_column: str,
    type_column: str,
    plate_column: str,
    nrows: int = 8,
    ncols: int = 12,
    plots_per_row: int = 3,
) -> Figure:
    """Plot an overview of what is in the wells of one or multiple plates

    This function visualizes the metadata of biological experiments on
    microtiter plates. One or multiple plates are shown and the wells are
    colored according to the type of sample they contain.

    Args:
        metadata_df:    dataframe with the metadata, containing at least the
                        columns `well_column`, `type_column`, `plate_column`
        well_column:    name of the column in `metadata_df` containing the
                        well number in the format of a single uppercase letter
                        (row) followed by a number (column), example: A01
        type_column:    name of the column in `metadata_df` containing the type
                        of sample the corresponding well contains. The function
                        is made for a categorical column, so ideally this column
                        contains strings
        plate_column:   name of the column in `metadata_df` containing the plate
                        identifiers. This column can be numerical or string
        nrows:          optional, number of rows of the plates used.
                        Default: 8 (96 well plate)
        ncols:          optional, number of columns of the plates used.
                        Default: 12 (96 well plate)
        plots_per_row:  optional, number of plots to put in one row when plotting
                        multiple plates. Default: 3

    Returns:
        Figure: figure containing the visualization

    """
    if well_column not in metadata_df.columns:
        raise ValueError(
            "The given well column name does not exist in the given dataframe!"
        )
    if type_column not in metadata_df.columns:
        raise ValueError(
            "The given type column name does not exist in the given dataframe!"
        )
    if plate_column not in metadata_df.columns:
        raise ValueError(
            "The given plate column name does not exist in the given dataframe!"
        )

    # number of plates to visualize
    p = len(pd.unique(metadata_df[plate_column]))

    # number of rows and columns needed in the figure
    fig_columns = min(p, plots_per_row)
    fig_rows = (
        p // plots_per_row
        if p // plots_per_row == p / plots_per_row
        else p // plots_per_row + 1
    )

    # prepare dataframe for plotting
    df = metadata_df[[plate_column, type_column, well_column]].copy()
    df = df.drop_duplicates()
    df[type_column] = df[type_column].astype(str)
    df["well_col_num"] = df[well_column].str[1:].astype(int)
    df["well_row_letter"] = df[well_column].str[:1]
    df["well_row_num"] = df.well_row_letter.apply(lambda letter: letter_to_int(letter))

    plate_names = pd.unique(df[plate_column])

    # create colormap
    color_labels = df[type_column].unique()
    rgb_values = (
        sns.color_palette("colorblind", len(color_labels))
        if len(color_labels) <= 10
        else sns.color_palette("hls", len(color_labels))
    )
    color_map = dict(zip(color_labels, rgb_values))

    # create plot for each plate
    fig = plt.figure(figsize=(7 * fig_columns, 5 * fig_rows))
    axes = []
    for plate in range(1, p + 1):
        ax = fig.add_subplot(fig_rows, fig_columns, plate)
        ax.set_xticks(list(range(1, ncols + 1)))
        ax.set_yticks(list(range(1, nrows + 1)))
        ax.set_yticklabels(alphabet[:nrows])
        one_plate_df = df[df[plate_column] == plate_names[plate - 1]]
        grouped = one_plate_df.groupby(type_column)
        for name, group in grouped:
            group.plot(
                x="well_col_num",
                y="well_row_num",
                kind="scatter",
                ax=ax,
                xlabel=f"plate {plate_names[plate-1]}",
                xlim=[0, ncols + 1],
                ylim=[0, nrows + 1],
                ylabel="",
                s=100,
                label=name,
                color=color_map[name],
            )
        ax.grid(False)
        ax.invert_yaxis()
        axes.append(ax)

    _ = fig.suptitle(
        f"Distribution of {type_column} over the different plates and wells"
    )

    # add legend
    handles, labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        mask = [x not in labels for x in l]
        mask = [i for i, x in enumerate(mask) if x]
        handles += [h[i] for i in mask]
        labels += [l[i] for i in mask]
        ax.get_legend().remove()
    leg = axes[0].legend(handles, labels, title=type_column)
    leg.set_bbox_to_anchor((-0.3, 0.8, 0.2, 0.2))

    return fig
