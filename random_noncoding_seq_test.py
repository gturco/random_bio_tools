import unittest
import collections

from flatfeature import Bed
from random_noncoding_seq import gene_size_dict,merged_feats,union,find,recursive_merge,non_coding,recursive_merge_both
from collections import deque

class TestPaper(unittest.TestCase):
    def setUp(self):
        self.tab_file = '/Users/gturco/Desktop/rice_set_size.tsv'
        self.size_dict = gene_size_dict(self.tab_file)
        self.bed_file = '/Users/gturco/data/paper3/rice_b_set.nolocaldups.with_new_cns_mask.bed'
    def test_dict(self):
        size_dict = gene_size_dict(self.tab_file)
        print size_dict["Os01g02200"]
        self.assertEqual(len(size_dict["Os01g02200"]), 13)
   # def test_merge_overlapping_feats(self):
   #     x = [(25651121, 25652114), (25646452, 25653625), (25652132, 25652179),(25674078, 25674126), (25665254, 25671356), (25663959,25663977), (25653294, 25653313), (25654427, 25658494),(25650556, 25651103), (25659144, 25667383), (25647439,25647454), (25659633, 25659704), (25672987, 25673116),(25647262, 25647279), (25673843, 25677732), (25650142,25650525), ( 5653712, 25653905), (25656436,25656492), (25675392, 25675467), (25652288, 25652586),(25659742, 25659762), (25678175, 25681280)]
   #     merged = merge_overlapping_feats(x,Bed(self.bed_file))
   #     print merged
#  #      self.assertEqual(38.17204301075269,gc_dict['q__6|24305193|24305378|scaffold_4|33902373|33902192|1e-51']['gc'])
    def test_find(self):
        x = (1,5)
        y = (3,20)
        self.assertEqual(True, find(x,y))

        x = (6,9)
        y = (30,40)
        self.assertEqual(False, find(x,y))

    def test_union(self):
        x = union((1,5),(3,20))
        self.assertEqual(x,deque([(1,20)]))

        x = union((6,9),(30,40))
        self.assertEqual(x,deque([(6,9), (30,40)]))

        x = union((1,20), (6,9))
        self.assertEqual(x,deque([(1,20)]))

        x = union((1,20), (30,40))
        self.assertEqual(x,deque([(1,20), (30,40)]))


    def test_merged_feats(self):
        feats = deque([(1,5),(3,20),(4,14),(6,9),(30,40),(60,80),(70,79),(100,110)])
        prv_count = 1
        new_feats = feats
        print "feats {0}".format(new_feats)
        while prv_count != len(new_feats) and len(new_feats) > 1:
        #    print "new_feats_old:{0}".format(new_feats)
            prv_count = len(new_feats)
            new_feats = merged_feats(new_feats)
            print len(feats)
        #print "new_feats :{0}".format(list(new_feats))
        expected = deque([(1,20),(30,40),(60,80),(100,110)])
        self.assertEqual(expected, new_feats)

    def test_recursive_merge(self):
        x = [(25646452, 25653625), (25650142, 25650525), (25650556,
           25651103),(25651121, 25652114), (25652132, 25652179),
           (25652288,25652586)]
        x.sort()
        expected = [(25646452,25653625)]
        results = recursive_merge_both(deque(x))
        self.assertEqual(results,expected)
        #print (recursive_merge(y))
    
    def test_non_coding(self):
        feats =[(30,33),(38,40),(55,60)]
        search_start = 0
        search_end = 68
        cns_size = 5
        non_coding_seq = non_coding(feats,cns_size,search_start,search_end)
        expected = [[41, 42, 43, 44, 45, 46, 47, 48, 49], [0, 1, 2, 3, 4, 5, 6,
            7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24], [61], [62], [63]]
        self.assertEqual(expected,non_coding_seq)
    def test_non_coding_geneends(self):
        feats =[(30,33),(38,40),(55,68)]
        search_start = 20
        search_end = 68
        cns_size = 5
        non_coding_seq = non_coding(feats,cns_size,search_start,search_end)
        expected = [[41, 42, 43, 44, 45, 46, 47, 48, 49], [20, 21, 22, 23, 24]]
        self.assertEqual(expected,non_coding_seq)
    def test_non_coding_large(self):
        feats = [(23044091,23044218),(23045020,23045059),(23042605,23042660)]
        feats.sort()
        search_start = 2304000
        search_end = 23045091
        cns_size = 30
        non_coding_seq = non_coding(feats,cns_size,search_start,search_end)
        #print non_coding_seq

    def test_unon(self):
        #x= [(5785295, 5787370), (5787723, 5788769), (5788850, 5788891),
        #        (5789359, 5797146), (5790143, 5792771), (5793155, 5793305)]
        x = [(1,2),(3,4),(5,7),(12,20),(13,24),(26,30)]
        y = recursive_merge_both(deque(x))
        print y
if __name__ == '__main__':
    unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
