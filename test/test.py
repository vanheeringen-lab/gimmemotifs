import unittest
from tempfile import NamedTemporaryFile
import os
from gimmemotifs.motif import *
from gimmemotifs.fasta import Fasta
from time import sleep


data_dir = "data/pwmscan"

with open(os.path.join(data_dir, "TATA.pwm")) as f:
    motif = read_motifs(f, fmt="pwm")[0]
print(motif)
prom = Fasta(os.path.join(data_dir, "promoters.fa"))
result = motif.pwm_scan(prom, nreport=1)
print(result)
