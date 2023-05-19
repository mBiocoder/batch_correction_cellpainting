import scarches
import pandas as pd
import urllib.request
from bs4 import BeautifulSoup
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from typing import Optional
import sys
import os
sys.path.append('../bcc/rxrx1-utils/rxrx/')
import io as rio


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


def visualize(data: str, experiment: str, plate: int, well: str, site: int, channels: Optional[list] = [1, 2, 3, 4, 5, 6], ax=None) -> np.ndarray:
    """Visualize the respective rxrx1 image given data, experiment, plate, well, and site information using the loa_site() function provided by Recursion.

    Args:
        data: "train" or "test" set.
        experiment: Experiment identifier.
        plate: Plate number.
        well: Well identifier.
        site: Site number.
        channels: List of channel numbers to visualize (default: [1, 2, 3, 4, 5, 6]).
        ax: Optional matplotlib Axes object to plot the image (default: None).

    Returns:
        The RGB image as a NumPy array.
    """

    img_tensor = rio.load_site(data, experiment, plate, well, site, channels)
    img = rio.convert_tensor_to_rgb(img_tensor)
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.imshow(img)
    ax.set_axis_off()
    return img



def visualize_cells_with_sirna(data: str, experiment: str, well: str, site: int, sirna: str, max_plate_number: int, plate: Optional[int] = None, channels: list = [1, 2, 3, 4, 5, 6]):
    """
    Visualizes cells with a specific siRNA treatment across different plates.
    
    Args:
        data: "train" or "test" set.
        experiment: Name of the experiment.
        well: Well identifier.
        site: Site identifier.
        sirna: siRNA treatment.
        max_plate_number: Maximum number of plates for this dataset.
        plate: Plate number. If None, visualize all plates.
        channels: List of channel numbers to visualize.

    Returns:
        None
    """
    well_size = 10
    title_fs = 25

    fig = plt.figure(figsize=(5 * well_size, 3 * well_size))

    if plate is None:
        for plate in range(1,  max_plate_number + 1):
            ax = fig.add_subplot(2, 4, plate)
            visualize(data, experiment, plate, well, site=site, channels=channels, ax=ax)
            ax.set_title(f'{experiment} Plate {plate}, siRNA {sirna}, well {well}, site {site}', fontdict={'fontsize': title_fs})

    else:
        y = rio.load_site_as_rgb(data, experiment, plate, well, site)
        plt.figure(figsize=(3, 3))
        plt.axis('off')
        plt.title(f'{experiment} Plate {plate}, siRNA {sirna}, well {well}, site {site}', fontdict={'fontsize': 7})
        _ = plt.imshow(y)

