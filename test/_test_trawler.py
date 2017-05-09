from __future__ import print_function
from motiftools.tools import *
t = Trawler()
print(t._run_program("/mnt/data/gimmemotifs/chen/CTCF/CTCF_prediction.fa", params={"background":"/mnt/data/gimmemotifs/chen/CTCF/CTCF_bg_genomic_matched.fa"}))

