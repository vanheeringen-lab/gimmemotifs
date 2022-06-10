import pytest

from gimmemotifs.preprocessing import (
    _load_peaks,
    _load_reads,
    combine_peaks,
    coverage_table,
    peakfile2bedfile,
    read_peak_file_to_df,
)


def test_peakfile2bedfile():
    peak_bed = peakfile2bedfile("test/data/filetype/bed.example2.bed")
    with open(peak_bed) as f:
        assert (
            next(f)
            == "Chr01\t18774479\t18774679\tp300-LELS-stage105-1_peak_1\t13.21507\t+\n"
        )
        assert (
            next(f)
            == "Chr01\t109122001\t109122201\tp300-LELS-stage105-1_peak_2\t24.28237\t+\n"
        )


def test__load_peaks():
    peak_bed = peakfile2bedfile("test/data/filetype/bed.example2.bed")
    regions = _load_peaks(peak_bed)
    assert len(regions) == 4
    assert regions[0] == "Chr01:18774479-18774679"


def test__load_reads():
    peak_bed = peakfile2bedfile("test/data/filetype/bed.example2.bed")
    regions = _load_peaks(peak_bed)
    datafile = "test/data/filetype/bed.example2.bed"
    df = _load_reads(peak_bed, datafile, regions)
    assert df.shape == (4, 1)
    assert df.columns == ["bed.example2"]
    assert df.index[0] == "Chr01:18774479-18774679"
    assert df.at["Chr01:18774479-18774679", "bed.example2"] == 1.0


def test_coverage_table():
    fname = "test/data/filetype/bed.example2.bed"
    df = coverage_table(fname, [fname], log_transform=False, ncpus=1)
    assert df.shape == (4, 1)
    assert df.columns == ["bed.example2"]
    assert df.index[0] == "Chr01:18774479-18774679"
    assert df.at["Chr01:18774479-18774679", "bed.example2"] == 1.0

    df = coverage_table(fname, [fname], log_transform=True, ncpus=1)
    assert df.shape == (4, 1)
    assert df.columns == ["bed.example2"]
    assert df.index[0] == "Chr01:18774479-18774679"
    assert round(df.at["Chr01:18774479-18774679", "bed.example2"], 2) == 0.69

    df = coverage_table(
        fname, [fname], log_transform=False, normalization="quantile", ncpus=1
    )
    assert df.shape == (4, 1)
    assert df.columns == ["bed.example2"]
    assert df.index[0] == "Chr01:18774479-18774679"
    assert df.at["Chr01:18774479-18774679", "bed.example2"] == 1

    # fname2 = "test/data/filetype/bed.example3.txt"
    df = coverage_table(fname, [fname, fname], top=2, ncpus=1)
    assert df.shape == (2, 2)

    df2 = coverage_table(fname, [fname, fname], top=2, ncpus=2)
    assert df2.shape == (2, 2)
    assert df.equals(df2)


def test_read_peak_file_to_df():
    narrowpeak = "test/data/filetype/narrowpeak.example1.narrowpeak"
    df = read_peak_file_to_df(narrowpeak)
    cols = ["chrom", "start", "end", "name", "value"]
    assert all(col in df.columns for col in cols)  # list(df.columns).__eq__(cols)  #
    assert df.shape == (10, len(cols) + 2)

    with pytest.raises(ValueError):
        # not a summits.bed (end-start>1)
        read_peak_file_to_df("test/data/filetype/bed.example2.bed")


def test_combine_peaks():
    # TODO: better test data, better tests
    peaks = ["test/data/filetype/narrowpeak.example1.narrowpeak"]
    genome = "test/data/scan/genome/scan_test.fa"
    df = combine_peaks(peaks, genome)
    for _, row in df.iterrows():
        assert int(row["end"]) - int(row["start"]) == 1
