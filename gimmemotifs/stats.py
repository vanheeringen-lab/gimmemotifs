"""Calculate motif enrichment statistics."""
from multiprocessing import Pool
import logging

import numpy as np
from scipy.stats import rankdata
import pandas as pd

from gimmemotifs import rocmetrics
from gimmemotifs.scanner import scan_to_best_match, Scanner
from gimmemotifs.motif import read_motifs, Motif
from gimmemotifs.config import MotifConfig
from gimmemotifs.utils import pfmfile_location

logger = logging.getLogger("gimme.stats")


def calc_stats_iterator(
    fg_file=None,
    bg_file=None,
    fg_table=None,
    bg_table=None,
    motifs=None,
    stats=None,
    genome=None,
    zscore=True,
    gc=True,
    ncpus=None,
):
    """Calculate motif enrichment metrics.

    Parameters
    ----------
    fg_file : str, optional
        Filename of a FASTA, BED or region file with positive sequences.

    bg_file : str, optional
        Filename of a FASTA, BED or region file with negative sequences.

    fg_table : str, optional
        Filename of a table with motif scan results of positive sequences.

    bg_table : str, optional
        Filename of a table with motif scan results of negative sequences.

    motifs : str, list or Motif instance, optional
        A file with motifs in pfm format, a list of Motif instances or a
        single Motif instance. If motifs is `None`, the default motif
        database is used.

    genome : str, optional
        Genome or index directory in case of BED/regions.

    stats : list, optional
        Names of metrics to calculate. See gimmemotifs.rocmetrics.__all__
        for available metrics.

    ncpus : int, optional
        Number of cores to use.

    Returns
    -------
    result : dict
        Dictionary with results where keys are motif ids and the values are
        dictionary with metric name and value pairs.
    """
    if not stats:
        stats = rocmetrics.__all__

    if fg_table is None:
        if fg_file is None:
            raise ValueError("Need either fg_table or fg_file argument")
    elif fg_file is not None:
        raise ValueError("Need either fg_table or fg_file argument, not both")

    if bg_table is None:
        if bg_file is None:
            raise ValueError("Need either bg_table or bg_file argument")
    elif bg_file is not None:
        raise ValueError("Need either bg_table or bg_file argument, not both")

    if fg_table is not None or bg_table is not None:
        remove_stats = []
        for s in stats:
            func = getattr(rocmetrics, s)
            if func.input_type == "pos":
                remove_stats.append(s)
        if len(remove_stats) != 0:
            logger.warn(
                "Cannot calculate stats that require position from table of motif scores."
            )
            logger.warn(f"Skipping the following statistics: {', '.join(remove_stats)}")
            stats = [s for s in stats if s not in remove_stats]

    if isinstance(motifs, Motif):
        all_motifs = [motifs]
    else:
        if type([]) == type(motifs):
            all_motifs = motifs
        else:
            motifs = pfmfile_location(motifs)
            all_motifs = read_motifs(motifs, fmt="pwm")
    if fg_table is not None or bg_table is not None:
        filtered_motifs = pd.read_csv(
            fg_table, sep="\t", index_col=0, nrows=1, comment="#"
        ).columns
        filtered_motifs = filtered_motifs.intersection(
            pd.read_csv(bg_table, sep="\t", index_col=0, nrows=1, comment="#").columns
        )
        all_motifs = [m for m in all_motifs if m.id in filtered_motifs]

    if ncpus is None:
        ncpus = int(MotifConfig().get_default_params()["ncpus"])

    if fg_file is not None or bg_file is not None:
        if zscore or gc:
            # Precalculate mean and stddev for z-score calculation
            s = Scanner(ncpus=ncpus)
            s.set_motifs(all_motifs)
            s.set_genome(genome)
            s.set_meanstd(gc=gc)

    chunksize = 240
    for i in range(0, len(all_motifs), chunksize):
        result = {}
        logger.debug(
            "chunk %s of %s", (i / chunksize) + 1, len(all_motifs) // chunksize + 1
        )
        motifs = all_motifs[i : i + chunksize]

        if fg_table is None:
            fg_total = scan_to_best_match(
                fg_file, motifs, ncpus=ncpus, genome=genome, zscore=zscore, gc=gc
            )
        else:
            fg_total = pd.read_csv(
                fg_table, sep="\t", usecols=[m.id for m in motifs], comment="#"
            ).to_dict(orient="list")
            for m in fg_total:
                fg_total[m] = [(x, None) for x in fg_total[m]]

        if bg_table is None:
            bg_total = scan_to_best_match(
                bg_file, motifs, ncpus=ncpus, genome=genome, zscore=zscore, gc=gc
            )
        else:
            bg_total = pd.read_csv(
                bg_table, sep="\t", usecols=[m.id for m in motifs], comment="#"
            ).to_dict(orient="list")
            for m in bg_total:
                bg_total[m] = [(x, None) for x in bg_total[m]]

        logger.debug("calculating statistics")

        if ncpus == 1:
            it = _single_stats(motifs, stats, fg_total, bg_total)
        else:
            it = _mp_stats(motifs, stats, fg_total, bg_total, ncpus)

        for motif_id, s, ret in it:
            if motif_id not in result:
                result[motif_id] = {}
            result[motif_id][s] = ret
        yield result


def calc_stats(
    fg_file=None,
    bg_file=None,
    fg_table=None,
    bg_table=None,
    motifs=None,
    stats=None,
    genome=None,
    zscore=True,
    gc=True,
    ncpus=None,
):
    """Calculate motif enrichment metrics.

    Parameters
    ----------
    fg_file : str
        Filename of a FASTA, BED or region file with positive sequences.

    bg_file : str
        Filename of a FASTA, BED or region file with negative sequences.

    fg_table : str
        Filename of a table with motif scan results of positive sequences.

    bg_table : str
        Filename of a table with motif scan results of negative sequences.

    motifs : str, list or Motif instance
        A file with motifs in pwm format, a list of Motif instances or a
        single Motif instance.

    genome : str, optional
        Genome or index directory in case of BED/regions.

    stats : list, optional
        Names of metrics to calculate. See gimmemotifs.rocmetrics.__all__
        for available metrics.

    ncpus : int, optional
        Number of cores to use.

    Returns
    -------
    result : dict
        Dictionary with results where keys are motif ids and the values are
        dictionary with metric name and value pairs.
    """
    result = {}
    for batch_result in calc_stats_iterator(
        fg_file=fg_file,
        bg_file=bg_file,
        fg_table=fg_table,
        bg_table=bg_table,
        motifs=motifs,
        genome=genome,
        stats=stats,
        ncpus=ncpus,
        zscore=zscore,
        gc=gc,
    ):
        for motif_id in batch_result:
            if motif_id not in result:
                result[motif_id] = {}
            for s, ret in batch_result[motif_id].items():
                result[motif_id][s] = ret
    return result


def _single_stats(motifs, stats, fg_total, bg_total):
    # Initialize multiprocessing pool

    for motif in motifs:
        motif_id = motif.id
        fg_vals = fg_total[motif_id]
        bg_vals = bg_total[motif_id]
        for s in stats:
            func = getattr(rocmetrics, s)
            if func.input_type == "score":
                fg = [x[0] for x in fg_vals]
                bg = [x[0] for x in bg_vals]
            elif func.input_type == "pos":
                fg = [x[1] for x in fg_vals]
                bg = [x[1] for x in bg_vals]
            else:
                raise ValueError("Unknown input_type for stats")

            ret = func(fg, bg)
            yield str(motif), s, ret


def _mp_stats(motifs, stats, fg_total, bg_total, ncpus):
    # Initialize multiprocessing pool
    pool = Pool(ncpus, maxtasksperchild=1000)

    jobs = []
    for motif in motifs:
        motif_id = motif.id
        fg_vals = fg_total[motif_id]
        bg_vals = bg_total[motif_id]
        for s in stats:
            func = getattr(rocmetrics, s)
            if func.input_type == "score":
                fg = [x[0] for x in fg_vals]
                bg = [x[0] for x in bg_vals]
            elif func.input_type == "pos":
                fg = [x[1] for x in fg_vals]
                bg = [x[1] for x in bg_vals]
            else:
                raise ValueError("Unknown input_type for stats")

            j = pool.apply_async(func, (fg, bg))
            jobs.append([str(motif), s, j])
    pool.close()
    pool.join()

    for motif_id, s, job in jobs:
        ret = job.get()
        yield motif_id, s, ret


def star(stat, categories):
    stars = 0
    for c in sorted(categories):
        if stat >= c:
            stars += 1
        else:
            return stars
    return stars


def add_star(stats):
    all_stats = {
        "mncp": [2, 5, 8],
        "roc_auc": [0.6, 0.75, 0.9],
        "max_enrichment": [10, 20, 30],
        "enr_at_fpr": [4, 8, 12],
        "fraction_fpr": [0.4, 0.6, 0.8],
        "ks_significance": [4, 7, 10],
        "numcluster": [3, 6, 9],
    }

    for motif, s2 in stats.items():
        for bg, s in s2.items():
            stats[motif][bg]["stars"] = int(
                np.mean([star(s[x], all_stats[x]) for x in all_stats.keys() if x in s])
                + 0.5
            )
    return stats


def rank_motifs(stats, metrics=("roc_auc", "recall_at_fdr")):
    """Determine mean rank of motifs based on metrics."""
    rank = {}
    combined_metrics = []
    motif_ids = stats.keys()
    background = list(stats.values())[0].keys()
    for metric in metrics:
        mean_metric_stats = [
            np.mean([stats[m][bg][metric] for bg in background]) for m in motif_ids
        ]
        ranked_metric_stats = rankdata(mean_metric_stats)
        combined_metrics.append(ranked_metric_stats)

    for motif, val in zip(motif_ids, np.mean(combined_metrics, 0)):
        rank[motif] = val

    return rank


def write_stats(stats, fname, header=None):
    """write motif statistics to text file."""
    # Write stats output to file

    for bg in list(stats.values())[0].keys():
        f = open(fname.format(bg), "w")
        if header:
            f.write(header)

        stat_keys = sorted(list(list(stats.values())[0].values())[0].keys())
        f.write("{}\t{}\n".format("Motif", "\t".join(stat_keys)))

        for motif in stats:
            m_stats = stats.get(str(motif), {}).get(bg)
            if m_stats:
                f.write(
                    "{}\t{}\n".format(
                        "_".join(motif.split("_")[:-1]),
                        "\t".join([str(m_stats[k]) for k in stat_keys]),
                    )
                )
            else:
                logger.warn(
                    "No stats for motif {0}, skipping this motif!".format(motif.id)
                )
            # motifs.remove(motif)
        f.close()

    return
