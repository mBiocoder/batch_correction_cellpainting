import unittest

import bcc
import pandas as pd


class utils_TestClass(unittest.TestCase):
    def setUp(self):
        self.scarches_version = bcc.utils.get_scarches_version()

    def test_scarches_over_0_1_0(self):
        major, minor, revision = [int(x) for x in self.scarches_version.split(".")]

        assert all([major >= 0, minor >= 1, revision >= 0])

    def test_sirna_target_webscraping(self):
        # this siRNA has no target annotated (as of 08.05.2023)
        x = bcc.utils.get_sirna_target_gene("s21506")
        assert pd.isnull(x)
        x = bcc.utils.get_sirna_target_gene("s21506", "entrez")
        assert pd.isnull(x)
        x = bcc.utils.get_sirna_target_gene("nonsense")
        assert pd.isnull(x)
        x = bcc.utils.get_sirna_target_gene("nonsense", "entrez")
        assert pd.isnull(x)

        # test gene name for one or multiple targets
        x = bcc.utils.get_sirna_target_gene("s21480")
        assert x == "SUB1"
        x = bcc.utils.get_sirna_target_gene("s21422")
        assert x == "DHRS4L2_DHRS4"

        # test entrez ID for one or multiple targets
        x = bcc.utils.get_sirna_target_gene("s21480", "entrez")
        assert x == "10923"
        x = bcc.utils.get_sirna_target_gene("s21422", "entrez")
        assert x == "317749_10901"
