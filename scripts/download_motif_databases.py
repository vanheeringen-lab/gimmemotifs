import os
import sys
import time
import urllib.request

from gimmemotifs.motif import read_motifs

HOMER_URL = "http://homer.ucsd.edu/homer/custom.motifs"
HOMER_NAME = "HOMER.pwm"
JASPAR_URL = "http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE{}_non-redundant_pfms_jaspar.txt"
JASPAR_NAME = "JASPAR2018{}.pwm"
JASPAR_GROUPS = ["", "vertebrates", "plants", "insects", "nematodes", "fungi", "urochordates"]

HOCOMOCO_URL = "http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/{0}/mono/HOCOMOCOv11_core_pcms_{0}_mono.txt"
HOCOMOCO_NAME = "HOCOMOCOv11_{}.pwm"

date = time.strftime("%Y-%m-%d")

outdir = "data/motif_databases"
if not os.path.exists(".git") or not os.path.exists(outdir):
    print("This script is for GimmeMotifs development only!")
    sys.exit(1)


### Homer ###
with open(os.path.join(outdir, HOMER_NAME), "w") as f:
    print("# Homer motif database", file=f)
    print("# Retrieved from: {}".format(HOMER_URL), file=f)
    print("# Version: {}".format(date), file=f)
    with urllib.request.urlopen(HOMER_URL) as response:
        for line in response:
            line = line.decode().strip()
            if line.startswith(">"):
                line = "_".join(line.split("\t")[:2])
            print(line, file=f)


### JASPAR ###
for group in JASPAR_GROUPS:
    if group != "":
        group = "_" + group
    outfile = os.path.join(outdir, JASPAR_NAME.format(group))
    url = JASPAR_URL.format(group)
    with open(outfile, "w") as f:
        with urllib.request.urlopen(url) as response:
            for line in response:
                line = line.decode().strip()
                if line.startswith(">"):
                    line = "_".join(line.split("\t")[:2])
                print(line, file=f)
        
    motifs = read_motifs(outfile, fmt="jaspar")
    with open(outfile, "w") as f:
        print("# JASPAR2018{} motif database".format(group), file=f)
        print("# Retrieved from: {}".format(url), file=f)
        print("# Version: {}".format(date), file=f)
        for motif in motifs:
            print(motif.to_pwm(), file=f)
 
### HOCOMOCO ###
for group in ["HUMAN", "MOUSE"]:
    outfile = os.path.join(outdir, HOCOMOCO_NAME.format(group))
    url = HOCOMOCO_URL.format(group)
    with open(outfile, "w") as f:
        print("# HOCOMOCOv11_{} motif database".format(group), file=f)
        print("# Retrieved from: {}".format(url), file=f)
        print("# Version: {}".format(date), file=f)
        with urllib.request.urlopen(url) as response:
            for line in response:
                line = line.decode().strip()
                if line.startswith(">"):
                    line = "_".join(line.split("\t")[:2])
                print(line, file=f)
