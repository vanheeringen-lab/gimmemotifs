import unittest
from gimmemotifs.db import MotifDb


class TestMotifDatabases(unittest.TestCase):
    """ A test class for motif database functionality """

    def setUp(self):
        pass

    def test_list(self):
        """Test motif db list"""

        l = sorted(MotifDb.list_databases())
        expected = [
            "cis-bp",
            "encode",
            "factorbook",
            "hocomoco",
            "homer",
            "image",
            "jaspar",
            "rsat",
            "swissregulon",
        ]
        self.assertEqual(l, expected)

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
