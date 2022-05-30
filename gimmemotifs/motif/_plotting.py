# Copyright (c) 2009-2021 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Scanning functions for Motif class"""

from warnings import warn


import logomaker as lm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_logo(
    self,
    fname=None,
    kind="information",
    title=True,
    ylabel=True,
    add_left=0,
    ax=None,
):
    """Plot motif logo

    Parameters
    ----------
    fname : str, optional
        If fname is set, the plot will be saved with fname as filename.
    kind : str, optional
        Type of logo to plot, can be 'information', 'frequency', 'energy' or
        'ensembl'.
    title : bool, optional
        Plot the motif id as the title.
    ylabel : bool, optional
        Plot the Y axis label.
    add_left : int, optional
        Add non-informative positions to the left (to align logo)
    """
    fig_height = 3
    fig_width = 0.45

    if add_left > 0:
        total = sum(self.pfm[0]) / 4
        pfm = np.vstack(([[total] * 4] * add_left, self.pfm))
    else:
        pfm = self.pfm

    matrix = pd.DataFrame(pfm, columns=["A", "C", "G", "T"])

    if kind == "ensembl":
        self.plot_ensembl_logo(fname=fname, title=title)
        return

    logo_params = {
        "information": {
            "df": lm.transform_matrix(
                matrix, from_type="probability", to_type="information"
            ),
            "figsize": (fig_width * matrix.shape[0], fig_height),
            "show_spines": False,
            "vpad": 0.02,
        },
        "frequency": {
            "df": lm.transform_matrix(
                matrix, from_type="probability", to_type="probability"
            ),
            "figsize": (fig_width * matrix.shape[0], fig_height),
            "show_spines": False,
            "vpad": 0.02,
            "font_name": "DejaVu Sans Mono",
        },
        "energy": {
            "df": lm.transform_matrix(
                lm.transform_matrix(matrix, from_type="probability", to_type="weight"),
                center_values=True,
            ),
            "figsize": (fig_width * matrix.shape[0], fig_height * 2),
            "fade_below": 0.7,
            "shade_below": 0.3,
            "flip_below": False,
            "show_spines": False,
        },
    }

    if ax is not None:
        logo_params[kind]["ax"] = ax
        del logo_params[kind]["figsize"]

    logo = lm.Logo(**logo_params[kind])

    if kind == "information":
        if ylabel:
            logo.ax.set_ylabel("Bits", fontsize=16)
        logo.ax.set_ylim(0, 2)
        logo.ax.set_yticks([0, 0.5, 1, 1.5, 2], minor=False)
    elif kind == "frequency":
        if ylabel:
            logo.ax.set_ylabel("Frequency", fontsize=16)
        logo.ax.set_ylim(0, 1)
    elif kind == "energy":
        if ylabel:
            logo.ax.set_ylabel(r"$\Delta \Delta G$/RT", labelpad=-1, fontsize=16)
    else:
        raise ValueError("Unknown motif visualization")

    logo.ax.set_xticks(range(matrix.shape[0]))
    logo.ax.set_xticklabels(range(1, matrix.shape[0] + 1))
    if title:
        logo.ax.set_title(self.id, fontsize=16)

    if fname:
        plt.savefig(fname, dpi=300)
        plt.close()
    else:
        return logo


def plot_ensembl_logo(self, fname=None, ic=True, title=True, letters=True, height=2):
    """Plot motif logo.

    This is an implementation of the logo presented here:
    http://www.ensembl.info/2018/10/15/new-ensembl-motif-features/

    Parameters
    ----------
    fname : str, optional
        If fname is set, the plot will be saved with fname as filename.
    ic : bool, optional
        Use the bit score. If this is set to False, the frequency
        will be used.
    title : bool, optional
        Plot the motif id as the title.
    letters : bool, optional
        Plot the nucleotides in the bars.
    height : float, optional
        Height of the plot.
    """
    width = 0.94

    ppm = self.ppm
    nucs = np.array(["A", "C", "G", "T"])

    neg_matrix = np.zeros((len(ppm), 4))

    pos_matrix = []
    nuc_ppm = []
    for row in ppm:
        if ic:
            ylabel = "bits"
            ic_row = []
            y_max = 2
            for p in row:
                if p < 0.25:
                    ic_row.append(0)
                else:
                    ic_row.append(p * np.log2((p) / 0.25))

        else:
            ic_row = row
            ylabel = "frequency"
            y_max = 1
        idx = np.argsort(ic_row)
        pos_matrix.append(np.array(ic_row)[idx])
        nuc_ppm.append(nucs[idx])

    colors = {
        "A": (0.308, 0.709, 0.280),
        "C": (0.145, 0.362, 0.6),
        "G": (0.969, 0.702, 0.172),
        "T": (0.841, 0.158, 0.224),
    }

    # x_max = np.max([np.sum(row) for row in pos_matrix])
    # x_min = -np.min([np.sum(row) for row in neg_matrix])
    # neg_matrix = neg_matrix / x_min * x_max

    plt.figure(figsize=(len(ppm) * 0.3, height))
    for (sign, matrix) in [(1, pos_matrix), (-1, neg_matrix)]:
        minbottom = np.zeros(len(matrix))
        alpha = 1
        if sign == -1:
            minbottom = np.array([np.sum(row) for row in matrix])
            alpha = 0.5

        # Print the bars
        for i in range(0, len(ppm[0])):

            pheight = [abs(r[i]) for r in matrix]
            bottom = minbottom + [sum(np.abs(r[:i])) for r in matrix]

            c = [colors[r[i]] for r in nuc_ppm]
            plt.bar(
                range(1, len(ppm) + 1),
                width=width,
                height=pheight,
                bottom=bottom,
                color=c,
                alpha=alpha,
            )

        if letters:
            # Print the letters
            for i in range(len(ppm)):
                for n in range(4):
                    x = i + 1
                    y = matrix[i][n] / 2 + sum(matrix[i][:n])
                    nuc = nuc_ppm[i][n]
                    c = "white"
                    if abs(matrix[i][n]) * height >= 0.5:
                        plt.text(
                            x,
                            y,
                            nuc,
                            horizontalalignment="center",
                            verticalalignment="center",
                            fontsize=8 + 10 * (matrix[i][n] / y_max),
                            color=c,
                        )

    # Remove axis lines
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_color("none")

    if title:
        plt.title(self.id)
    plt.xlim(0.47, len(self) + 0.5)
    plt.xticks(range(1, len(ppm) + 1))
    plt.ylabel(ylabel)

    if fname:
        plt.savefig(fname, dpi=300)
        plt.close()
    else:
        return ax


def to_img(self, fname, fmt="PNG", add_left=0, seqlogo=None, height=6):
    """Create a sequence logo using seqlogo.

    Create a sequence logo and save it to a file. Valid formats are: PNG,
    EPS, GIF and PDF.

    Parameters
    ----------
    fname : str
        Output filename.
    fmt : str , optional
        Output format (case-insensitive). Valid formats are PNG, EPS, GIF
        and PDF.
    add_left : int , optional
        Pad motif with empty positions on the left side.
    seqlogo : str
        Location of the seqlogo executable. By default the seqlogo version
        that is included with GimmeMotifs is used.
    height : float
        Height of the image
    """
    warn("Method to_img() is replaced by plot_logo()", DeprecationWarning)
    self.plot_logo(fname=fname, kind="information")
