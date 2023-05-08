import scarches
import pandas as pd
import urllib.request
from bs4 import BeautifulSoup
from typing import Literal


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
