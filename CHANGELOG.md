# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Unreleased]

## [0.12.0] - 2018-07-10

### Added

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

- Fixed warning of `ks_pvalue` when `p == 0`.
- Fixed issue with nested multiprocessing pools.
- Fix numpy version because of DeprecationWarning in sklearn.
- Updated xgboost dependency, where the API had changed.
