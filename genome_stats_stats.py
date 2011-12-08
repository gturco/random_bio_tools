import unittest
import collections

import sys
import itertools
import logging

logging.basicConfig(level=logging.INFO)
sys.path.append("../scripts")
from flatfeature import Bed
from paper import gc_stats, stats

class TestPaper(unittest.TestCase):
    def setUp(self):
        self.fasta = 'data/rice_v6.features.fasta'
    def test_gc_stats(self):
        gc_dict = gc_stats(self.fasta,'>')
        #print gc_dict['Os09g39960']
        self.assertEqual(float(0), gc_dict['Os01g01100']['gc'])
        self.assertEqual(48.431685273790535, gc_dict['Os09g39960']['gc'])

    def test_gc_stats_query_only(self):
        gc_dict = gc_stats('data/rice_v6_setaria64/rice_v6_setaria64.cns_real.fasta','>q')
        self.assertEqual(38.17204301075269,gc_dict['q__6|24305193|24305378|scaffold_4|33902373|33902192|1e-51']['gc'])
    def test_gc_stats_subject_only(self):
        gc_dict = gc_stats('data/rice_v6_setaria64/rice_v6_setaria64.cns_real.fasta', '>s')
        self.assertEqual(37.362637362637365,gc_dict['s__6|24305193|24305378|scaffold_4|33902373|33902192|1e-51']['gc'])

    def test_stats(self):
        # did only ran rices ....
        stats('data/rice_v6_setaria64/rice_v6.features.fasta','>')
        stats('data/rice_v6_setaria64/rice_v6_setaria64.cns_real.fasta', '>q')
        stats('data/rice_t_sorghum_v1/rice_t_sorghum_v1.cns_real.fasta', '>q')
        stats('data/rice_t_sorghum_v1/rice_t.features.fasta', '>')
if __name__ == '__main__':
    unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
