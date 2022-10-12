# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## Unreleased


## [0.17.2] - 2022-10-12

### Changed

- made xgboost an optional dependency (to save space on bioconda)
- an existing config will now update available tools when accessed (e4b3275)
- applied the bioconda patch to compile_externals.py (11b0c2c)
- `coverage_table` and `combine_peaks` have their positional arguments under positional arguments (20819ee)
- `coverage_table` should be slightly faster now (20819ee)

### Fixed

- biofluff dependency back in requirements
- pinned conda and mamba versions in `.travis.yaml`
  - temp fix until conda>=4.12 can install mamba properly
- documentation is working again!
- gimmemotifs now supports pandas >=1.30

### Removed

- pyarrow dependency


## [0.17.1] - 2022-06-02

### Added

- `requirements.yaml` contains all conda dependencies.
  - packages available from one channel have been pinned (for solving speed)
  - packages have minimum versions where known (for solving speed)

### Changed

- alphabetized tools everywhere (how could you live like that!?)
- updated `setup.py`
- updated installation instructions

### Fixed

- Yamda is now recognized in the config
- most tools work with the editable installation again 
- all tests work for unix
  - there were still some flakey values, where randomness is involved.
- background.py updated to work with the specified minimum `genomepy` version
- all `sphinx-build docs build` warnings
- motifs require to have unique ids when clustering, thanks @akmorrow13!
- motif2factors removes apostrophes so it wont crash :)
- removed a print

### Removed

- a bunch of redundant requirement files.
- OSX tests. Possibly temporary.
  - The tests haven't working for ages, so I have no idea where to begin.
  - and Travis asks 5x credits for OSX machines...


## [0.17.0] - 2021-12-22

### Added

* Added `--genomes_dir` argument to `gimme motif2factors`.
* Added `--version` flag.
* Function `sample()` for fast sequence sampling from a `Motif()` instance.
* Added JASPAR 2022 motif databases.
* Updated Homer motif database.
* Operators:
  * `+` - take the combination of two motifs (average), based on pfm, which means that motifs with higher counts will be weighed more heavily.
  * `&` - take the combination of two motifs (average), based on the ppm, which means that both motifs will be weighed equally.
  * `<<` - "shift" motif left (adding a non-informative position to the right side)
  * `>>` - "shift" motif right (adding a non-informative position to the left side)
  * `~` - reverse complement
  * `*` - multiply the pfm by a value
* Progress bar for scanning.
* `list_installed_libraries()` to list available motif libraries.

### Changed

* `Motif()` class completely restructured:
  * Split into multiple files with coherent function.
  * Uses `numpy.array` internally.
  * All functions that mention `pwm` renamed to `ppm` (position-probability matrix), as the definition of a PWM is usually a log-odds matrix, not a probability matrix.
    * `to_pwm()` is deprecated, use `to_ppm()` instead.
  * Changed functions `pwm_min_score()` and `pwm_max_score()` to properties `max_score` and `min_score`.
  * All internal data is correctly updated when `Motif()` is changed, for instance by trimming (#218).


### Fixed

* `gimme motif2factors` can now unzip genome fastas.
* `gimme motif2factors` will sanitize genome names.
* Fixed bugs related to partial rerun of `gimme motif2factors`.
* Fixed unhandled `OSError` during installation on Mac.
* Fixed bug related to `RFE()` (#226).
* Positional probability matrix now sum to 1 over all positions (#209).
* Fixed issue with pandas >= 1.3.
* Fixed issue with `non_reducing_slice` import from pandas.
* Fix threshold calculation if more than 20,000 sequences are supplied.
* Fix issue with config file getting corrupted.
* Fix FPR threshold calculation.

### Removed


## [0.16.1] - 2021-06-28

Bugfix release.

### Added

* Added warning when the number of sequences used for de novo motif prediction is low.

### Fixed

* Fixed bug with `gimme motif2factors`.
* Fixed "Motif does not occur in motif database when running maelstrom" (#192).
* Fixed bugs related to runs where no (significant) motifs is found.

## [0.16.0] - 2021-05-28

Many bugfixes, thanks to @kirbyziegler, @irzhegalova, @wangmhan, @ClarissaFeuersteinAkgoz and @fgualdr for reporting and proposing solutions!
Thanks to @Maarten-vd-Sande for the speed improvements.

### Added

* `gimme motif2factors` command to annotate a motif database with TFs from different species
  based on orthogroups.
* Informative error message with link to fix when cache is corrupted (running on a cluster).
* Print an informative error message if the input file is not in the correct format.

### Changed

* Speed improvements to motif scanning, which is now up to 2X faster!
* Size of input regions is now automatically adjusted (#123, #128, #129)
* Quantile normalization in `coverage_table` now uses multiple CPUs.

### Fixed

* Fixes issue where % of motif occurence would be incorrectly reported in `gimme maelstrom` output (#162).
* Fix issues with running Trawler (#181)
* Fix issues with running YAMDA (#180)
* Fix issues with parsing XXmotif output (#178)
* Fix issue where command line argument (such as single strand) are ignored (#177)
* Fix pyarrow dependency (#176)
* The correct % of regions with motif is now reported (#162)
* Fix issue with running `gimme motifs` with the HOMER database (#135)
* Fix issue with the `--size` parameter in `gimme motifs`, which now works as expected (#128)

## [0.15.3] - 2021-02-01

### Fixed

* `_non_reducing_slice` vs `non_reducing_slice` for pandas>=1.2 (#168)
* When using original region size, skip regions smaller than 10bp and warn if no
  regions are left. 
* Fixed creating statistics report crashed with `KeyError: 'Factor'` (#170)
* Fixed bug with creating GC bins for a genome with unusual GC% (like Plasmodium).
* Fixed bug that occurs when upgrading pyarrow with an existing GimmeMotifs
  cache.


## [0.15.2] - 2020-11-26

### Changed

- Refactoring to make `coverage_table` and `combine_peaks` available via API.

### Fixed

- Fix issue with -s parameter of `gimme motifs` (#146)
- Fix issues (hopefully) with scanning large input files.

## [0.15.1] - 2020-10-07

### Added

- `Motif.plot_logo()` accepts an `ax` argument.

### Fixed

- Support for pandas>=1.1
- `coverage_table` doesn't add a newline at the end of the file.

## [0.15.0] - 2020-09-29

### Added

- Added additional columns to `gimme maelstrom` output for better intepretation (correlation of motif to signal and % of regions with motif).
- Added support for multi-species input in `genome@chrom:start-end` format.
- `gimme maelstrom` warns if data is not row-centered and will center by default.
- `gimme maelstrom` selects a set of non-redundant (or less redundant) motifs by default.
- Added SVR regressor for `gimme maelstrom`.
- Added quantile normalization to `coverage_table`.

### Removed

- Removed the lightning classifiers and regressors as the package is no longer actively maintained.

### Changed

- Visually improved HTML output.
- Score of `maelstrom` is now an aggregate z-score based on combining z-scores from individual methods using Stouffer's method. The z-scores of individual methods are generated using the inverse normal transform.
- Reorganized some classes and functions.

### Fixed

- Fixed minor issues with sorting columns in HTML output.
- `gimme motifs` doesn't crash when no motifs are found.
- Fixed error with Ensembl chromosome names in `combine_peaks`.

## [0.14.4] - 2020-04-02

### Fixed

- Fixed "TypeError: an integer is required (got type str)" when creating GC index.
- Fixed `combine_peaks` with Ensembl chromosome names (thanks @JGAsmits). 
- Fixed bug with pandas>=1.0.

## [0.14.3] - 2020-02-19

### Fixed

- Fixed 'AttributeError: can't delete attribute' in `gimme maelstrom` and `gimme motifs` (#108, #109).

## [0.14.2] - 2020-01-31

Bugfix release

### Added

- The `combine_peaks` script now supports `.narrowPeak` files.

### Removed

- Removed seqlogo dependency.
- Removed obsolete code.

### Changed

- Refactored the `tools` section.

### Fixed

- Configuration issue with `size` instead of `width` (#103).
- Updated `tqdm` requirement (#98).


## [0.14.1] - 2019-12-19

Bugfix release

### Fixed

- Fix function for locating a pwm/pfm motif database.
- Added configparser dependency

## [0.14.0] - 2019-12-05

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
