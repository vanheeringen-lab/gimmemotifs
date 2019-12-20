from gimmemotifs.fasta import *
import tempfile
import os
import pytest


@pytest.fixture()
def fasta_file():
    return "test/data/fasta/test.fa"


@pytest.fixture()
def id_fa_files():
    return ["test/data/fasta/test2.fa", "test/data/fasta/test3.fa"]


@pytest.fixture()
def fasta_obj(fasta_file):
    return Fasta(fasta_file)


def test1_index(fasta_obj):
    """ Fasta as a dictionary """
    assert fasta_obj["seq1"] == "AAAA"
    assert fasta_obj["seq2"] == "ACGT"
    assert fasta_obj["seq3"] == "CCCCGGGG"


def test2_items(fasta_obj):
    """ Fasta.items() """
    assert len(list(fasta_obj.items())) == 3


def test3_writefasta(fasta_file, fasta_obj):
    """ Write fasta-formatted file"""
    temp = tempfile.NamedTemporaryFile()
    tempname = temp.name
    fasta_obj.writefasta(tempname)
    with open(fasta_file) as f:
        with open(tempname) as f_ref:
            assert f.read().strip() == f_ref.read().strip()


def test4_gt_in_id(id_fa_files):
    for id_fa_file in id_fa_files:
        f = Fasta(id_fa_file)
        seqs = [seq for seq in f.seqs]
        print(seqs)
        for i in range(3):
            assert seqs[i] == "ACTG" * (i + 1)
