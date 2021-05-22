"""
Make a new motif2factors file based on orthology. This gets complicated
because not all factors in the motif2factors file are based on gene names;
some factors are gene ids, and some are aliases or other symbols.
"""
import os
import re
import json
import shutil
import sys
import urllib
import sqlite3
import pathlib
import tempfile
import subprocess
import multiprocessing.dummy
from typing import List
from textwrap import wrap
from functools import lru_cache
from urllib.error import HTTPError

import pyfaidx
import genomepy
import numpy as np
import pandas as pd
from loguru import logger
from genomepy.utils import get_genomes_dir

from gimmemotifs.motif import read_motifs


FASTA_LINEWIDTH = 80
BLACKLIST_TFS = [
    "dobox4",  # does not exist
    "dobox5",  # does not exist
    "foxa",  # not sure which FOXA(1,2,3)
    "ep300",
]
RENAME_TFS = {"NR1A4": "NR4A1", "SREBP1a": "SREBF1"}


def motif2factor_from_orthologs(
    database: str = "gimme.vertebrate.v5.0",
    database_references: List[str] = None,
    extra_orthologs_references: List[str] = None,
    new_reference: List[str] = None,
    tmpdir: str = None,
    outdir: str = ".",
    strategy: str = "lenient",
    threads: int = 24,
):
    """
    Make a motifs2factors file based on gene orthology.

    This function first downloads the genomes of the new reference, the old
    reference (human and mouse), and a wide range of different vertebrate
    genomes (a range of genomes is necessary to get a better ortholog inference)

    Then the peptide sequence of each gene is extracted from the genome +
    annotation, and orthologs between these are derived by orthofinder.

    Finally, based on these orthologs, a new motifs2factors file is created from
    the old one.
    """
    if not database_references:
        database_references = ["GRCh38.p13", "GRCm38.p6"]

    if not extra_orthologs_references:
        extra_orthologs_references = [
            "danRer11",  # zebrafish
            "UCB_Xtro_10.0",  # xenopus
            "GRCg6a",  # chicken
            "BraLan2",  # lancet fish (out-group ish)
            "oryLat2",  # medaka
            "ARS-UCD1.2",  # cow
            "phaCin_unsw_v4.1",  # koala
            "rCheMyd1.pri",  # turtle
        ]

    _check_install()

    tmpdir = tempfile.mkdtemp() if tmpdir is None else tmpdir

    logger.info(f"Making a new reference for: {' & '.join(new_reference)}.")
    logger.info(
        f"You say the {database} database is based on: {' & '.join(database_references)}."
    )
    logger.info(
        f"For better orthology inference we are also using these assemblies: {' & '.join(extra_orthologs_references)}."
    )
    logger.info(f"Using strategy: {strategy} for orthology/name inference.")
    logger.info(f"tmpdir: {tmpdir}.")
    logger.info(f"outdir: {outdir}.")

    all_genomes = set(database_references + extra_orthologs_references + new_reference)

    # download all required genomes
    logger.info("Downloading all assemblies.")
    _download_genomes_with_annot(all_genomes, tmpdir)

    # convert each genome + annotation into the primary genes (longest protein per gene)
    logger.info("Taking the longest protein per gene per assembly.")
    for genome in all_genomes:
        annot2primpep(genome, tmpdir)

    # run orthofinder on our primary genes
    logger.info("Running orthofinder to find orthologs.")
    orthofinder_result = _orthofinder(f"{tmpdir}/prim_genes", threads)

    # now parse the output of orthofinder
    logger.info("Storing everything in a database.")
    orthogroup_db = f"{tmpdir}/orthologs.sqlite"
    load_orthogroups_in_db(orthogroup_db, all_genomes, orthofinder_result)

    # get all motifs and related factors from motif database
    motifs = read_motifs(database)
    motifsandfactors = {
        motif.id: [val for sublist in motif.factors.values() for val in sublist]
        for motif in motifs
    }

    # process the references
    logger.info("Now writing your motif2factors files.")
    for genome in new_reference:
        make_motif2factors(
            f"{outdir}/{genome}.{database}",
            new_reference=genome,
            database_references=database_references,
            motifsandfactors=motifsandfactors,
            database=orthogroup_db,
            strategy=strategy,
            motifs=motifs,
        )

    # cleanup
    shutil.rmtree(tmpdir)


def _check_install():
    dependencies = ["gffread", "orthofinder"]
    if any(shutil.which(dependency) is None for dependency in dependencies):
        logger.warning(
            f"""Running gimme motif2factors requires {" & ".join(dependencies)} to be installed.

You can easily install these dependencies with conda in the environment where gimmemotifs is installed:
conda install {" ".join(dependencies)}"""
        )
        sys.exit(1)


def _orthofinder(peptide_folder, threads):
    """
    Run orthofinder on the peptide folder
    """
    # run orthofinder on the primary transcripts
    result = subprocess.run(
        ["orthofinder", "-f", peptide_folder, "-t", threads], capture_output=True
    )

    logger.debug(f"""stdout of orthofinder:\n {result.stdout.decode("utf-8")}""")
    logger.debug(f"""stderr of orthofinder:\n {result.stderr.decode("utf-8")}""")

    orthofinder_result = re.search(
        "Results:\n    (.*)", result.stdout.decode("utf-8")
    ).group(1)
    return orthofinder_result


def _download_genomes_with_annot(genomes, genomes_dir):
    # make sure each genome has a gene annotation
    no_annotations = [genome for genome in genomes if not _has_annotation(genome)]
    assert (
        len(no_annotations) == 0
    ), f"genome(s): {','.join(no_annotations)} seem not to have an annotation for it."

    # download the genomes
    # add check to see if not already in genomes dir?
    default_genomes_dir = get_genomes_dir(check_exist=False)
    logger.debug(f"using default genome dir: {default_genomes_dir}")
    for genome in genomes:
        # check if already in default genomes dir, if so, skip downloading and directly copy
        if all(
            os.path.exists(f"{default_genomes_dir}/{genome}/{genome}.{extension}")
            for extension in ["fa", "annotation.gtf"]
        ):
            logger.info(f"{genome} was already downloaded, using that version.")
            # if except, probably a continuation from previous run
            try:
                os.mkdir(f"{genomes_dir}/{genome}")
                os.symlink(
                    f"{default_genomes_dir}/{genome}/{genome}.fa",
                    f"{genomes_dir}/{genome}/{genome}.fa",
                )
                os.symlink(
                    f"{default_genomes_dir}/{genome}/{genome}.annotation.gtf",
                    f"{genomes_dir}/{genome}/{genome}.annotation.gtf",
                )
            except Exception:
                pass
        else:
            logger.info(f"Downloading {genome} through genomepy.")
            genomepy.install_genome(genome, annotation=True, genomes_dir=genomes_dir)
            result = subprocess.run(
                ["gunzip", f"{genomes_dir}/{genome}/{genome}.annotation.gtf.gz"],
                capture_output=True,
            )
            logger.debug(f"""stdout of genomepy:\n {result.stdout.decode("utf-8")}""")
            logger.debug(f"""stderr of genomepy:\n {result.stderr.decode("utf-8")}""")


def annot2primpep(genome, outdir):
    """
    Go from a gene annotation (gtf) to a primary proteins file. Uses gffread to
    do the heavy lifting.

    The gene identifier of the fasta consists of both the gene name and the
    gene id, separated by a pipe symbol, like so:

    >GATA4|ENSG00000136574
    MYQSLAMAANHGPPPGAYEAGGPGAFMHGAGAASSPVYVPTPR
    GTQQGSPGWSQAGADGAAYTPPPVSPRFSFPGTTGSLAAAAAA
    """
    # setup our result folder
    pathlib.Path(f"{outdir}/prim_genes").mkdir(parents=True, exist_ok=True)

    # for each genome, make a .pep.fa. This pep.fa contains the LONGEST protein for each gene
    # use gffread to convert our annotation.gtf into all possible peptides
    result = subprocess.run(
        [
            "gffread",
            "-y",
            f"{outdir}/{genome}/{genome}.pep.fa",
            "-g",
            f"{outdir}/{genome}/{genome}.fa",
            f"{outdir}/{genome}/{genome}.annotation.gtf",
            "-S",
            "--table",
            "gene_name,gene_id",
        ],
        capture_output=True,
    )
    logger.debug(f"""stdout of gffread:\n {result.stdout.decode("utf-8")}""")
    logger.debug(f"""stderr of gffread:\n {result.stderr.decode("utf-8")}""")

    # read the resulting .pep.fa
    proteins = pyfaidx.Fasta(f"{outdir}/{genome}/{genome}.pep.fa", read_long_names=True)

    # and only keep the longest protein per gene
    records = dict()
    for record in proteins:
        # get the gene name and gene id of the protein
        prot_name = "|".join(record.long_name.split("\t")[1:])

        # get the protein sequence
        protein = str(record)

        # skip proteins with stop codon
        # (probably mitochondrial protein since they have a different codon table)
        if "*" in protein:
            logger.debug(
                f"skipping {prot_name} since it contains a * symbol (probably mitochondrial read)."
            )
            continue

        # skip when we already have a longer edition of the same gene
        if prot_name in records and len(records[prot_name]) > len(protein):
            continue

        # add to our records
        records[prot_name] = protein

    # and finally save our result
    with open(f"{outdir}/prim_genes/{genome}.pep.fa", "w") as f:
        for record, sequence in records.items():
            f.write(f">{record}\n")
            f.write("\n".join(wrap(sequence, width=FASTA_LINEWIDTH)))
            f.write("\n")


def load_orthogroups_in_db(db, genomes, orthofinder_result):
    """
    Save the results of orthofinder (tsv file) in a relational database
    (SQLite). This makes it possible to easily switch between genes and
    orthogroups.

    We create three tables; genes, orthogroups, and assemblies.
    """
    orthogroups = pd.read_table(
        f"{orthofinder_result}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
    )
    unassigned = pd.read_table(
        f"{orthofinder_result}/Orthogroups/Orthogroups_UnassignedGenes.tsv"
    )

    # remove old db if it exists (TODO: shouldnt be necessary?)
    if os.path.exists(db):
        os.remove(db)
    conn = sqlite3.connect(db)
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS orthogroups
        (
        id integer PRIMARY KEY,
        orthogroup text UNIQUE NOT NULL
        )
        """
    )

    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS assemblies
        (
        id integer PRIMARY KEY,
        assembly text
        )
        """
    )

    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS genes
        (
        id integer PRIMARY KEY,
        gene_name text,
        gene_name_lower text,
        gene_id text,
        gene_id_lower text,
        assembly NOT NULL,
        orthogroup,
        FOREIGN KEY (orthogroup) REFERENCES orthogroups (id),
        FOREIGN KEY (assembly) REFERENCES assemblies (assembly)
        )
        """
    )

    # fill up our database
    # all orthogroups
    for orthogroup in orthogroups["HOG"]:  # <-- all Hierarchical OrthoGroups
        conn.execute(f"INSERT INTO orthogroups VALUES(NULL, '{orthogroup}')")
    conn.execute("INSERT INTO orthogroups VALUES(NULL, 'UNASSIGNED')")

    # all assemblies
    for assembly in genomes:
        conn.execute(f"INSERT INTO assemblies VALUES(NULL, '{assembly}')")

    # all genes per species
    for assembly in genomes:
        for orthogroup, genes in zip(
            orthogroups["HOG"], orthogroups[f"{assembly}.pep"]
        ):
            if isinstance(genes, float) and np.isnan(genes):
                continue

            for gene in genes.split(", "):
                # for each gene we store its gene_name (if present,
                # otherwise .) and gene_id. We load them both in the database.
                gene_name, gene_id = gene.split("|")
                if gene_name != ".":
                    conn.execute(
                        f"INSERT INTO genes VALUES(NULL, '{gene_name}', '{gene_name.lower()}', '{gene_id}', '{gene_id.lower()}', '{assembly}', '{orthogroup}')"
                    )
                else:
                    conn.execute(
                        f"INSERT INTO genes VALUES(NULL, NULL, NULL, '{gene_id}', '{gene_id.lower()}', '{assembly}', '{orthogroup}')"
                    )

    # also add unasigned genes
    for assembly in genomes:
        for gene in unassigned[f"{assembly}.pep"]:
            if isinstance(gene, float) and np.isnan(gene):
                continue

            # for each gene we store its gene_name (if present,
            # otherwise .) and gene_id. We load them both in the database.
            gene_name, gene_id = gene.split("|")
            if gene_name != ".":
                conn.execute(
                    f"INSERT INTO genes VALUES(NULL, '{gene_name}', '{gene_name.lower()}', '{gene_id}', '{gene_id.lower()}', '{assembly}', NULL)"
                )
            else:
                conn.execute(
                    f"INSERT INTO genes VALUES(NULL, NULL, NULL, '{gene_id}', '{gene_id.lower()}', '{assembly}', NULL)"
                )

    conn.commit()
    cur = conn.cursor()

    # now make some indices to make it fast!
    cur.execute("CREATE INDEX idx_name_lower ON genes (gene_name_lower, assembly)")
    cur.execute("CREATE INDEX idx_id_lower ON genes (gene_id_lower, assembly)")
    cur.execute("CREATE INDEX idx_orthogroup_assembly ON genes (orthogroup, assembly)")
    conn.commit()

    return db


def make_motif2factors(
    prefix,
    new_reference,
    database_references,
    motifsandfactors,
    database,
    strategy,
    motifs,
):
    """
    Make a motifs2factors file based on an existing ortholog database.

    We loop over each motif, find which factors belong to it in the old
    reference, then find all the orthogroups those belong to, and finally assign
    all the genes of the new reference species that belong to those orthogroups.
    """
    factor2orthogroups.cache_clear()
    conn = sqlite3.connect(database)
    cur = conn.cursor()

    with open(f"{prefix}.motif2factors.txt", "w") as f:
        f.write("Motif\tFactor\tEvidence\tCurated\n")
        for motif, factors in motifsandfactors.items():
            # remember if we found any orthologs for all the factors belonging
            # to our motif
            motif_set = False

            # get all the orthogroups belonging to all our factors belonging to
            # our motif. This is not as trivial as it sounds, since gene naming
            # is a mess
            orthologousgroups = {
                item
                for factor in factors
                for item in factor2orthogroups(
                    factor, tuple(database_references), database, strategy
                )
            }
            if len(orthologousgroups) > 0:
                for orthologousgroup in orthologousgroups:
                    res = cur.execute(
                        f"SELECT * FROM genes WHERE orthogroup='{orthologousgroup}' AND assembly='{new_reference}'"
                    ).fetchall()

                    if len(res):
                        for _, gene_name, _, gene_id, *_ in res:
                            motif_set = True
                            if gene_name is not None:
                                f.write(f"{motif}\t{gene_name}\tOrthologs\tN\n")
                            else:
                                f.write(f"{motif}\t{gene_id}\tOrthologs\tN\n")
            if not motif_set:
                f.write(f"{motif}\tNO ORTHOLOGS FOUND\tOrthologs\tN\n")

    with open(f"{prefix}.pfm", "w") as f:
        for motif in motifs:
            print(motif.to_pwm(), file=f)


@lru_cache(maxsize=99999)
def factor2orthogroups(factor, references, database, strategy):
    """
    Returns all the orthogroups a factor belongs to. This is NOT a trivial task.

    Some of the factors shouldn't be considered (blacklisted), some have an
    ambiguous name that should be renamed, and in general the factor naming is
    inconsistent.

    We first check whether or not it fits in our orthogroup database (easy!), if
    not, we start querying mygeneinfo to get more gene aliases to search for in
    our orthogroups. If that didn't work we query mygeneinfo again, but with a
    **very** lenient query, and hope that works.
    """
    if factor.lower() in BLACKLIST_TFS:
        return []

    if factor in RENAME_TFS:
        factor = RENAME_TFS[factor]

    # see if the factor is in the original database
    if strategy in ["strict", "medium", "lenient"]:
        orthogroups = _factor2orthogroups_sql(factor, references, database)

    # if not, check if we can find anything by quering mygene.info
    if strategy in ["medium", "lenient"]:
        if len(orthogroups) == 0:
            gene_symbols = _unknownfactor2symbols(factor, fields=["name", "symbol"])
            for gene_symbol in gene_symbols:
                orthogroups += _factor2orthogroups_sql(
                    gene_symbol, references, database
                )

    # if still none found, we try a very broad myinfo query (higher chance of false positive)
    if strategy in ["lenient"]:
        if len(orthogroups) == 0:
            gene_symbols = _unknownfactor2symbols(
                factor,
                fields=[
                    "alias",
                    "other_names",
                    "accession",
                    "accession.protein",
                    "refseq",
                    "refseq.protein",
                    "ensembl",
                    "ensembl.gene",
                ],
            )
            for gene_symbol in gene_symbols:
                orthogroups += _factor2orthogroups_sql(
                    gene_symbol, references, database
                )

    orthogroups = list(set(orthogroups))

    if len(orthogroups) == 0:
        logger.info(f"No orthologs found for TF {factor}.")

    return orthogroups


def _factor2orthogroups_sql(factor, references, database):
    """
    Query our sql database for all orthogroups a factor belongs to
    """
    conn = sqlite3.connect(database)
    cur = conn.cursor()

    orthogroups = list()
    res = cur.execute(
        "SELECT orthogroup FROM genes "
        "WHERE ("
        + " OR ".join(f"assembly='{assembly}'" for assembly in references)
        + ")"
        + f"    AND (gene_name_lower='{factor.lower()}' OR gene_id_lower='{factor.lower()}')"
    ).fetchall()

    for orthogroup in res:
        if orthogroup is not None:
            orthogroups.append(*orthogroup)
    return orthogroups


def _unknownfactor2symbols(factor, fields):
    """
    Query mygeneinfo for different aliases of our gene
    """

    def mygeneinfo(field):
        request = f"http://mygene.info/v3/query?q={field}:{factor}&fields=entrezgene,name,symbol,taxid,other_names,ensembl.gene"
        try:
            with urllib.request.urlopen(request) as url:
                data = json.loads(url.read().decode())
                return data["hits"]
        except HTTPError:
            return list()

    # multithreaded to do multiple queries at once
    p = multiprocessing.dummy.Pool(len(fields))
    hits = p.map(mygeneinfo, fields)
    hits = [item for sublist in hits for item in sublist]
    symbols = set()

    # only keep aliases, gene names and HGNC symbols
    for hit in hits:
        for field in ["alias", "name", "symbol", "ensembl"]:
            if field in hit:
                if "'" not in hit[field]:
                    if field == "ensembl" and "gene" in hit[field]:
                        symbols.add(hit[field]["gene"])
                    else:
                        symbols.add(hit[field])
                elif isinstance(hit[field], list):
                    symbols.update(
                        {each["gene"] for each in hit[field] if "gene" in each}
                    )
    return list(symbols)


def _has_annotation(genome):
    """
    Return True if a genome has an annotation for it, else False
    """
    providers = [
        genomepy.ProviderBase.create(provider)
        for provider in ["ensembl", "ucsc", "ncbi"]
    ]
    for provider in providers:
        if genome in provider.genomes:
            link = provider.get_annotation_download_link(genome)
            return link is not None
    return False
