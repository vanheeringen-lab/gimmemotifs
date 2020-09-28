# import mygene
import pandas as pd
import pybedtools
from genomepy import Genome
import sys
from gimmemotifs.fasta import Fasta

# mg = mygene.MyGeneInfo()

# xli = ["gata3"]

# out = mg.querymany(xli, scopes="symbol", fields="genomic_pos", species="all")

# for hit in out:
#     print(hit)
# #     if "genomic_pos" in hit:
# #         print("{}:{}-{}\t{}".format(
# #             hit["genomic_pos"]["chr"],
# #             hit["genomic_pos"]["start"],
# #             hit["genomic_pos"]["end"],
# #             hit["query"],
# #             ))

# sys.exit()
from functools import singledispatch


@singledispatch
def scan(obj):
    # default implementation
    raise NotImplementedError(f"Not implemented for {type(obj)}")


@scan.register(pd.DataFrame)
def _scan_dataframe(df, columns=["chrom", "start", "end"], genome="hg38"):
    if not set(columns).issubset(df.columns):
        raise ValueError(f"Expected columns {columns}")

    if len(columns) == 3:
        # Assume this is chromosome start, end
        g = Genome(genome)
        seqs = list(
            (
                df[columns[0]]
                + ":"
                + df[columns[1]].astype(str)
                + "-"
                + df[columns[2]].astype(str)
            ).values
        )
        return g.track2fasta(seqs)
    elif len(columns) == 1:
        # Assume this is some kind of gene_id
        return df[columns[0]].values


# @scan.register(pybedtools.BedTool)
# @profile
def _scan_bedtool(bed, genome="hg38"):
    g = Genome(genome)
    intervals = [g[f.chrom][f.start : f.stop] for f in bed]
    return intervals


# @profile
def _scan_bedtool2(bed, genome="hg38"):
    g = Genome(genome)
    return Fasta(bed.sequence(fi=g.filename).seqfn).seqs


import requests

rest_url = "https://rest.ensembl.org/info/species"
r = requests.get(rest_url, headers={"Content-Type": "application/json"})

if not r.ok:
    r.raise_for_status()

json = r.json()

print(json)


# #'9:120165822-120173708, 'ensemblgene': 'MGP_SPRETEiJ_G0033934', 'start': 120165822,
# df = pd.DataFrame({"chrom":["chr1", "chr2"], "start":[100, 1000], "end":[200, 200]})
# b = pybedtools.BedTool("5k.bed")
# # for f in b:
# #     print(f)
# #     break
# seqs = _scan_bedtool(b, genome="Spur_3.1")
# seqs = _scan_bedtool2(b, genome="Spur_3.1")

# g = Genome("hg19")
# print(g["chr1"][1000000:1000100])

# FASTA file
# BED file
# region file
# Gene file
#   - promoter (all species)
#   - closest accessible region (human)
#   - sum / mean / max of regions within distance of promoter
#   - sum / mean / max of regions within weighted distance of promoter
