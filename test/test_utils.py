import unittest

import bcc

class utils_TestClass(unittest.TestCase):
    def setUp(self):
        
        self.scarches_version = bcc.utils.get_scarches_version()
        
    def test_scarches_over_0_1_0(self):
        
        major, minor, revision = [int(x) for x in self.scarches_version.split(".")]
        

        assert all([major >= 0, minor >= 1, revision >= 0])