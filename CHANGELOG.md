# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added

- The `gimme motifs` command supports new *de novo* motif prediction tools: DREME, ProSampler, YAMDA, DiNAMO and RPMCMC.
- A set of non-redundant motifs is selected using recursive feature elimination for `gimme motifs`.
- Motif scan results for non-redundant motifs are now included in the output.
- Plot motif logos using different styles of visualization (information content, frequency, energy or [Ensembl](http://www.ensembl.info/2018/10/15/new-ensembl-motif-features/)).
- New scoring scheme that uses a z-score based on genomic sequences with a similar GC%. 
- CIS-BP motif database version 2.0 ([Lambert et al. 2019](https://www.nature.com/articles/s41588-019-0411-1)).
- JASPAR 2020 motif database.

### Changed

- Only three motif prediction tools are selected by default for `gimme motifs`: BioProspector, Homer and MEME.
- The `xl` motif setting (motif width 6-20) is selected by default for `gimme motifs`. 
- Output names for files and reports of `gimme motifs` are now more consistent.
- Command line tools `gimme roc` has been removed. This functionality has now been merged with `gimme motifs`. The `gimme motifs` command now scans for both known and *de novo* motifs.
- Replaced weblogo/seqlogo with [logomaker](https://logomaker.readthedocs.io/). 
- GC% background now uses a genomic index of GC% frequencies.
- GC% z-score is now the default score for all `gimme` tools.

### Fixed

- FASTA files with `>` symbols in the header are now correctly parsed.
- Fixed GC% background. The GC% background is now correct for all input formats. 
- `factorial` import for scipy >= 1.3.0.
- All output columns in HTML reports can now be sorted correctly.

### Removed

- Deprecated modules and scripts.

## [0.13.1] - 2018-12-04

### Added

- Improved docstrings of several modules.
- Added new API examples.

### Fixed

- The MEME motif tools should now be recognized after install.
- Output of MEME 5.0.2 is now parsed correctly.
- If the inputfile of `gimme motifs` is not recognized, a clear error message is printed.
- Duplicate factors are removed from the motif factors list.

### Changed

- MEME is no longer included with GimmeMotifs. When installing via conda meme will be included. If GimmeMotifs is installed via pip, then MEME needs to be installed separately. 
- Changed "user" background to "custom" background.
- Updated Posmo to run with a wider variety of settings.

## [0.13.0] - 2018-11-19

### Added

- Multiple other motif databases (JASPAR, HOMER, HOCOMOCO, CIS-BP, ENCODE,
  Factorbook, SwissRegulon, IMAGE).
- Helper script to combine peaks (summit files from MACS2)
- Helper script to create coverage table (similar to bedtools multicov)
- Option to report z-score normalized motif scores.
- Added precision-recall AUC to stats and `gimme roc`.
- `gimme motifs` now supports narrowPeak input.
- Updated documentation with an explanation of the score that `gimme maelstrom` reports.

### Changed

- The `maelstrom` tools now use z-score normalized motif scores.
- Improved efficiency of motif scanning (>10X speed improvement).
- Removed dependency on R for rank aggregation.
- Dropped support for Python 2.
- Use versioneer for versioning.
- Removed the default genome in config file.
- Config file is now independent from GimmeMotifs version and will be created by
  default on first use.
- Simplified setup.py script.
- Updated parameters for ChIPMunk motif finder.

### Fixed 

- Fixed the seqcor similarity metric to use a non-random sequence and to also take
  the reverse complement of motif 2 into account.
- Improved the speed of `gimme roc`.
- Fixed memory leak of `gimme roc`.
- Fixed `scale` for newer `pandas`/`sklearn` combo
- FIxed bug related to backgroundgradient with new pandas

## [0.12.0] - 2018-07-10

**Please note:** the way GimmeMotifs uses genome FASTA files has changed in a
major way. It is no longer necessary to index genomes. GimmeMotifs now uses the
`faidx` index. You can use [genomepy](http://github.com/simonvh/genomepy) to
manage your genomes.

### Added

- Release checklist (for developer use).
- Helpful message when the genome cannot be found.
- Added a `-N/--nthreads` option to command-line tools to specify the number of
  threads independent from the config file.
- Added a `npcus` argument to many functions.
- Support for the `feather` data format as input for `gimme maelstrom`.
- Progress bar and resume option for maelstrom.
- Motif database is copied to maelstrom output directory.
- Added `MaelstromResult` class to work with maelstrom output.

### Changed

- Moved rank aggregation to separate function.
- Parallel implementation of rank aggregation.
- Updated documentation on how to add genomes in the new version.
- Switched to using `faidx` for FASTA indexing and `genomepy` for genome
  management.

### Removed

- The `gimme index` script (replaced by `genomepy`).

### Fixed

- Upload with correct description to PyPi.
- Fixed warning of `ks_pvalue` when `p == 0`.
- Fixed issue with nested multiprocessing pools.
- Fix numpy version because of DeprecationWarning in sklearn.
- Updated xgboost dependency, where the API had changed.
