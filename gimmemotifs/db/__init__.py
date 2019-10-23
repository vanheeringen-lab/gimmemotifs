"""Motif databases."""
import glob
import os
from urllib.request import urlopen, urlretrieve
import re
import time
import tarfile
import shutil
from tempfile import mkdtemp
import zipfile

import pandas as pd

from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import get_jaspar_motif_info

DEFAULT_OUT = "data/motif_databases/"


class MotifDb(object):

    """MotifDb base class.

    Use to get a list of available databases:
    >>> MotifDb.list_databases()
    []
    """

    _dbs = {}
    name = None

    date = time.strftime("%Y-%m-%d")

    @classmethod
    def create(cls, name):
        """Create a motif database object based on the db name.

        Parameters
        ----------
        name : str
            Name of the provider (eg. JASPAR, HOMER, ...)

        Returns
        -------
        db : MotifDb instance
            MotifDb instance.
        """
        try:
            return cls._dbs[name.lower()]()
        except KeyError:
            raise Exception("Unknown motif database")

    @classmethod
    def register_db(cls, dbname):
        """Register method to keep list of dbs."""

        def decorator(subclass):
            """Register as decorator function."""
            cls._dbs[dbname] = subclass
            subclass.name = dbname
            return subclass

        return decorator

    @classmethod
    def list_databases(self):
        """List available databases."""
        return self._dbs.keys()

    def __hash__(self):
        return hash(str(self.__class__))

    def create_annotation(self, name, anno):
        base = os.path.splitext(name)[0]
        fname = base + ".motif2factors.txt"
        with open(fname, "w") as f:
            print("Motif\tFactor\tEvidence\tCurated", file=f)
            for motif, factors in anno.items():
                for factor, status, curated in factors:
                    print(
                        "{}\t{}\t{}\t{}".format(motif, factor, status, curated), file=f
                    )


register_db = MotifDb.register_db


@register_db("jaspar")
class JasparMotifDb(MotifDb):
    """
    JASPAR motif database
    """

    URL = "http://jaspar.genereg.net/download/CORE/JASPAR{0}_CORE{1}_non-redundant_pfms_transfac.txt"
    UNVALIDATED_URL = "http://jaspar.genereg.net/download/collections/JASPAR{}_UNVALIDATED{}_pfms_transfac.txt"
    NAME = "JASPAR{}{}.pfm"
    GROUPS = [
        "",
        "vertebrates",
        "plants",
        "insects",
        "nematodes",
        "fungi",
        "urochordates",
    ]

    def download(self, version="2020", outdir=DEFAULT_OUT):
        # JASPAR
        for group in self.GROUPS:
            if group != "":
                group = "_" + group
            outfile = os.path.join(outdir, self.NAME.format(version, group))

            for i, base_url in enumerate([self.URL, self.UNVALIDATED_URL]):
                url = base_url.format(version, group)
                if i == 0:
                    mode = "w"
                else:
                    mode = "a"
                with open(outfile, mode) as f:
                    with urlopen(url) as response:
                        for line in response:
                            line = line.decode().strip()
                            if line.startswith(">"):
                                line = "_".join(line.split("\t")[:2])
                            print(line, file=f)

            motifs = read_motifs(outfile, fmt="transfac")

            anno = self.annotate_factors(motifs)
            with open(outfile, "w") as f:
                print("# JASPAR{}{} motif database".format(version, group), file=f)
                print("# Retrieved from:", file=f)
                for base_url in [self.URL, self.UNVALIDATED_URL]:
                    print("#     * {}".format(base_url.format(version, group)), file=f)
                print(
                    "# Note: this file also contains the unvalidated motifs from JASPAR.",
                    file=f,
                )
                print(
                    "#       These have not been confirmed by orthogonal evidence and have ",
                    file=f,
                )
                print("#       a motif id that starts with UN.", file=f)
                print("# Date: {}".format(self.date), file=f)
                for motif in motifs:
                    print(motif.to_pwm(), file=f)

            # if group == "_vertebrates":
            self.create_annotation(
                os.path.join(outdir, self.NAME.format(version, group)), anno
            )

    def annotate_factors(self, motifs):
        anno = {}
        for motif in motifs:
            mtype = "Unknown"
            if hasattr(motif, "metadata") and "data_type" in motif.metadata:
                mtype = motif.metadata["data_type"]
            else:
                info = get_jaspar_motif_info(motif.id.split("_")[0])
            try:
                mtype = info["type"]
            except Exception:
                pass
            if mtype == "universal protein binding microarray (PBM)":
                mtype = "PBM"
            factors = re.sub(r"\([^)]+\)", "", motif.id.split("_")[1]).split("::")
            if motif.id.startswith("MA"):
                direct = "Y"
            else:
                direct = "N"
            anno[motif.id] = [[f, mtype, direct] for f in factors]
        return anno


@register_db("homer")
class HomerMotifDb(MotifDb):
    """
    HOMER motif database
    """

    NAME = "HOMER.pfm"
    URL = "http://homer.ucsd.edu/homer/custom.motifs"

    def download(self, outdir=DEFAULT_OUT):
        # Homer
        pfm_out = os.path.join(outdir, self.NAME)
        with open(pfm_out, "w") as f:
            print("# Homer motif database (v4.10)", file=f)
            print("# Retrieved from: {}".format(self.URL), file=f)
            print("# Date: {}".format(self.date), file=f)
            with urlopen(self.URL) as response:
                for line in response:
                    line = line.decode().strip()
                    if line.startswith(">"):
                        line = "_".join(line.split("\t")[:2])
                    print(line, file=f)

        motifs = read_motifs(pfm_out)
        anno = self.annotate_factors(motifs)
        self.create_annotation(os.path.join(outdir, self.NAME), anno)

    def annotate_factors(self, motifs):
        anno = {}
        p = re.compile(r"\w+_([\w.-]+)\(([\w,-]+\))")
        for motif in motifs:
            name, source, _ = motif.id.split("/")
            try:
                m = p.search(name)
                name_factor = m.group(1).lower()
                source_factor = source.split("-")[1]
                for tag in ["gfp", "v5", "biotin", "myc"]:
                    source_factor = source_factor.replace("." + tag, "")

                if name_factor.replace("-", "") == source_factor.lower():
                    anno[motif.id] = [[source_factor, "ChIP-seq", "Y"]]
                else:
                    pass
            except Exception:
                pass
            # anno[motif.id] = factors
        return anno


@register_db("hocomoco")
class HocomocoMotifDb(MotifDb):
    """
    HOCOMOCO motif database
    """

    # BASE_URL = "http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/{0}/mono/"
    # ANNO_URL = BASE_URL + "HOCOMOCOv11_core_annotation_{0}_mono.tsv"
    # URL = BASE_URL + "HOCOMOCOv11_core_pcms_{0}_mono.txt"
    # NAME = "HOCOMOCOv11_{}.pfm"

    BASE_URL = "http://hocomoco10.autosome.ru/final_bundle/{0}/mono/"
    ANNO_URL = BASE_URL + "HOCOMOCOv10_annotation_{0}_mono.tsv"
    URL = BASE_URL + "HOCOMOCOv10_pcms_{0}_mono.txt"
    NAME = "HOCOMOCOv10_{}.pfm"

    def download(self, outdir=DEFAULT_OUT):
        for group in ["HUMAN", "MOUSE"]:
            outfile = os.path.join(outdir, self.NAME.format(group))
            url = self.URL.format(group)
            with open(outfile, "w") as f:
                print("# HOCOMOCOv10_{} motif database".format(group), file=f)
                print("# Retrieved from: {}".format(url), file=f)
                print("# Date: {}".format(self.date), file=f)
                with urlopen(url) as response:
                    for line in response:
                        line = line.decode().strip()
                        if line.startswith(">"):
                            line = "_".join(line.split("\t")[:2])
                        print(line, file=f)
            motifs = read_motifs(outfile)
            anno = self.annotate_factors(motifs, self.ANNO_URL.format(group))
            self.create_annotation(os.path.join(outdir, self.NAME.format(group)), anno)

    def annotate_factors(self, motifs, url):
        anno = {}

        with urlopen(url) as response:
            for line in response:
                vals = line.decode().strip().split("\t")
                anno[vals[0]] = vals[1]

        anno = {motif.id: [[anno[motif.id], "ChIP-seq", "Y"]] for motif in motifs}
        return anno


@register_db("encode")
class EncodeMotifDb(MotifDb):
    """
    ENCODE motif database
    Kheradpour and Kellis, 2013, doi:10.1093/nar/gkt1249
    """

    URL = "http://compbio.mit.edu/encode-motifs/motifs.txt"
    NAME = "ENCODE.pfm"

    def download(self, outdir=DEFAULT_OUT):
        outfile = os.path.join(outdir, self.NAME)
        with open(outfile, "w") as f:
            print("# ENCODE motif database", file=f)
            print("# Retrieved from: {}".format(self.URL), file=f)
            print("# Date: Dec. 2013", file=f)
            with urlopen(self.URL) as response:
                for line in response:
                    line = line.decode().strip()
                    if line.startswith(">"):
                        line = line.replace("\t", " ")
                    print(line, file=f)
        motifs = read_motifs(outfile)
        anno = self.annotate_factors(motifs)
        self.create_annotation(os.path.join(outdir, self.NAME), anno)

    def annotate_factors(self, motifs):
        anno = {}
        for motif in motifs:
            if "disc" in motif.id:
                source = motif.id.split(" ")[-1]
                factor = source.split("_")[0]
                anno[motif.id] = [[factor, "ChIP-seq", "N"]]
            elif "jolma" in motif.id:
                vals = motif.id.split(" ")
                factor = vals[-1].split("_")[0]
                anno[motif.id] = [[factor, "HT-SELEX", "N"]]
            elif "transfac" in motif.id:
                vals = motif.id.split(" ")
                factor = vals[-2].split("_")[0]
                anno[motif.id] = [[f, "TRANSFAC", "Y"] for f in factor.split("::")]
            elif "jaspar" in motif.id:
                vals = motif.id.split(" ")
                jaspar_factor = vals[-1].split("_")[0]
                factor = vals[-2].split("_")[0]
                if (
                    len(jaspar_factor) > 1
                    and jaspar_factor[1] == jaspar_factor[1].lower()
                ):
                    factor = factor.capitalize()
                anno[motif.id] = [[f, "JASPAR", "Y"] for f in factor.split("::")]
            elif "bulyk" in motif.id:
                vals = motif.id.split(" ")
                factor = vals[-2].split("_")[0].capitalize()
                bulyk_factor = vals[-1].split("_")[0].capitalize()
                if factor != bulyk_factor and factor.startswith("Znf"):
                    factor = bulyk_factor
                anno[motif.id] = [[factor, "PBM", "N"]]
            else:
                raise ValueError("Don't recognize motif {}".format(motif.id))
        return anno


@register_db("factorbook")
class FactorbookMotifDb(MotifDb):
    """
    Factorbook
    """

    ANNO_URL = (
        "https://genome.cshlp.org/content/suppl/2012/08/22/22.9.1798.DC1/TableS1.xls"
    )
    NAME = "factorbook.pfm"

    def download(self, outdir=DEFAULT_OUT):
        # Factorbook is only supplied in non-redundant form as a supplemental pdf
        # For now, use the non-redundant version included with GimmeMotifs
        infile = "data/motif_databases/factorbook.pfm"
        outfile = os.path.join(outdir, self.NAME)
        motifs = read_motifs(infile)
        with open(outfile, "w") as f:
            for motif in motifs:
                print(motif.to_pwm(), file=f)
        anno = self.annotate_factors(motifs)
        self.create_annotation(os.path.join(outdir, self.NAME), anno)

    def annotate_factors(self, motifs):
        anno = {}
        df = pd.read_excel(
            "https://genome.cshlp.org/content/suppl/2012/08/22/22.9.1798.DC1/TableS1.xls"
        )
        t = {}
        for factor, motif_names in (
            df[["HGNC ID", "canonical motif"]].dropna().drop_duplicates().values
        ):
            for m in motif_names.split(";"):
                t[m] = t.get(m, []) + [factor]

        for motif in motifs:
            name = motif.id.split(".")[0]
            if name in t:
                for factor in t[name]:
                    anno[motif.id] = anno.get(motif.id, []) + [
                        [factor, "ChIP-seq", "N"]
                    ]
        return anno


@register_db("swissregulon")
class SwissregulonMotifDb(MotifDb):
    """
    SwissRegulon
    """

    URL = "http://swissregulon.unibas.ch/data/hg19_f5/hg19_weight_matrices_v2"
    ANNO_URL = "http://swissregulon.unibas.ch/data/hg19_f5/hg19_mat_TF_associations.txt"
    # URL = "http://swissregulon.unibas.ch/data/hg19/weight_matrices"
    # ANNO_URL = "http://swissregulon.unibas.ch/data/hg19/mat_TF_associations.hg"
    NAME = "SwissRegulon.pfm"

    def download(self, outdir=DEFAULT_OUT):
        outfile = os.path.join(outdir, self.NAME)
        with open(outfile, "w") as f:
            with urlopen(self.URL) as response:
                for line in response:
                    line = line.decode().strip()
                    print(line, file=f)

        motifs = read_motifs(outfile, fmt="transfac")
        with open(outfile, "w") as f:
            print("# SwissRegulon motif database (hg19:FANTOM5)", file=f)
            print("# Retrieved from: {}".format(self.URL), file=f)
            print("# Date: {}".format(self.date), file=f)
            for motif in motifs:
                if len(motif) > 0:
                    print(motif.to_pwm(), file=f)

        motifs = read_motifs(outfile)
        anno = self.annotate_factors(motifs)
        self.create_annotation(os.path.join(outdir, self.NAME), anno)

    def annotate_factors(self, motifs):
        anno = {}

        with urlopen(self.ANNO_URL) as response:
            for line in response:
                line = line.decode().strip()
                # print(line)
                motif, *factors = line.split("\t")
                factors = [f.split(":")[2] for f in factors[1:]]
                for factor in factors:
                    anno[motif] = anno.get(motif, []) + [[factor, "Unknown", "N"]]
        return anno


@register_db("image")
class ImageMotifDb(MotifDb):
    """
    IMAGE
    """

    URL = "http://bioinformatik.sdu.dk/solexa/webshare/IMAGE/IMAGE_v1.1.tar.gz"
    NAME = "IMAGE.pfm"

    def download(self, outdir=DEFAULT_OUT):
        tmpdir = mkdtemp()
        file_tmp = urlretrieve(self.URL, filename=None)[0]
        tar = tarfile.open(file_tmp)
        fname = "IMAGE/utils/Collection.motif"
        members = [tar.getmember(fname)]
        tar.extractall(tmpdir, members=members)
        outfile = os.path.join(outdir, self.NAME)

        motifs = read_motifs(os.path.join(tmpdir, fname))
        with open(outfile, "w") as f:
            print("# IMAGE motif database (v1.1)", file=f)
            print("# Retrieved from: {}".format(self.URL), file=f)
            print("# Date: {}".format(self.date), file=f)
            for motif in motifs:
                print(motif.to_pwm(), file=f)
        shutil.rmtree(tmpdir)

        motifs = read_motifs(outfile)
        anno = self.annotate_factors(motifs)
        self.create_annotation(os.path.join(outdir, self.NAME), anno)

    def annotate_factors(self, motifs):
        anno = {}
        tmpdir = mkdtemp()
        file_tmp = urlretrieve(self.URL, filename=None)[0]
        tar = tarfile.open(file_tmp)
        fname = "IMAGE/utils/Genename_Motif.txt"
        members = [tar.getmember(fname)]
        tar.extractall(tmpdir, members=members)
        with open(os.path.join(tmpdir, fname)) as f:
            for line in f:
                vals = line.strip().split("\t")
                if len(vals) == 3:
                    factor, motif, status = vals
                    anno[motif] = anno.get(motif, []) + [[factor, status, "N"]]
        shutil.rmtree(tmpdir)
        return anno


@register_db("cis-bp")
class CisbpMotifDb(MotifDb):
    """
    CIS-BP 2.00
    """

    VERSION = "2.00"
    BASE = "http://cisbp.ccbr.utoronto.ca/data/{}/DataFiles/Bulk_downloads/EntireDataset/".format(  # noqa: E501
        VERSION
    )
    ANNO_URL = BASE + "/TF_Information_all_motifs.txt.zip"
    URL = BASE + "/PWMs.zip"
    NAME = "CIS-BP.pfm"

    def download(self, outdir=DEFAULT_OUT):
        tmpdir = mkdtemp()
        file_tmp = urlretrieve(self.URL, filename=None)[0]

        with zipfile.ZipFile(file_tmp, "r") as zip_ref:
            zip_ref.extractall(tmpdir)

        motifs = []
        for fname in glob.glob(os.path.join(tmpdir, "pwms/*")):
            m_id = os.path.splitext(os.path.basename(fname))[0]
            for m in read_motifs(fname, fmt="transfac"):
                if len(m) > 0:
                    m.id = m_id
                    motifs.append(m)
        outfile = os.path.join(outdir, self.NAME)
        with open(outfile, "w") as f:
            print("# CIS-BP motif database (v{})".format(self.VERSION), file=f)
            print("# Retrieved from: {}".format(self.URL), file=f)
            print("# Date: {}".format(self.date), file=f)
            for motif in motifs:
                print(motif.to_pwm(), file=f)

        shutil.rmtree(tmpdir)

        motifs = read_motifs(outfile)
        anno = self.annotate_factors(motifs)
        self.create_annotation(os.path.join(outdir, self.NAME), anno)

    def annotate_factors(self, motifs):
        anno = {}
        df = pd.read_table(self.ANNO_URL)
        df = df.loc[
            df["TF_Species"].isin(["Homo_sapiens"]) & (df["TF_Status"] != "N"),
            ["Motif_ID", "TF_Name", "MSource_Type", "TF_Status"],
        ]
        df["TF_Status"] = df["TF_Status"].str.replace("D", "Y").str.replace("I", "N")
        df["MSource_Type"] = df["MSource_Type"].str.replace("HocoMoco", "ChIP-seq")
        df = df.drop_duplicates()
        df = df.set_index("Motif_ID")
        df.columns = ["Factor", "Evidence", "Curated"]
        df = df.loc[df.index.intersection([m.id for m in motifs])].dropna()

        for m_id, row in df.iterrows():
            anno[m_id] = anno.get(m_id, []) + [row]

        return anno


@register_db("rsat")
class RsatMotifDb(MotifDb):
    """
    RSAT clustered motifs
    """

    URL = "http://pedagogix-tagc.univ-mrs.fr/rsat/data/published_data/Castro_2016_matrix-clustering/Application_4/{}/cor0.8_Ncor0.65/All_{}_motifs_cluster_root_motifs.tf"  # noqa: E501
    NAME = "RSAT_{}.pfm"

    def download(self, outdir=DEFAULT_OUT):
        for tax in ["insects", "plants", "vertebrates"]:
            tax_ = tax
            if not tax.endswith("es"):
                tax_ = tax[:-1]
            url = self.URL.format(tax.capitalize(), tax_)
            print(url)
            name = self.NAME.format(tax)

            file_tmp = urlretrieve(url, filename=None)[0]
            motifs = read_motifs(file_tmp, fmt="transfac")
            outfile = os.path.join(outdir, name)
            with open(outfile, "w") as f:
                print("# RSAT non-redundant {} motif database".format(tax), file=f)
                print("# Retrieved from: {}".format(url), file=f)
                print("# Date: {}".format(self.date), file=f)
                for motif in motifs:
                    print(motif.to_pwm(), file=f)

            anno = self.annotate_factors(motifs)
            self.create_annotation(os.path.join(outdir, self.NAME.format(tax)), anno)

    def annotate_factors(self, motifs):
        anno = {}
        return anno
