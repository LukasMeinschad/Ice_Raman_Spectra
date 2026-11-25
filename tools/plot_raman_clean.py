#!/usr/bin/env python3
"""Clean presentation figure: weighted-average calculated Raman spectrum vs. experiment,
top unscaled, bottom scaled (stretch x s_OH, librations x s_lib, translations unscaled).

    python plot_raman_clean.py ../Ice_XIV/pbe0_d3_raman --weights 1:0.437,2:0.323,3:0.240 \
        --exp-stretch ../reference/salzmann_2006_iceXIV_fig1a_H2O.csv \
        --exp-lattice ../reference/salzmann_2006_iceXIV_fig1c_H2O.csv
"""
import argparse
import csv
import os

import numpy as np

from analyze_raman import exp_factor, parse_raman
from plot_raman_avg_vs_exp import GRID, norm, read_xy, spectrum


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("jobs")
    p.add_argument("--weights", required=True)
    p.add_argument("--exp-stretch", required=True)
    p.add_argument("--exp-lattice", required=True)
    p.add_argument("--exp-label", default="Experiment (Salzmann et al. 2006)")
    p.add_argument("--calc-label", default="PBE0-D3 (Boltzmann average)")
    p.add_argument("--phase", default="Ice XIV")
    p.add_argument("--s-oh", type=float, default=0.9615)
    p.add_argument("--s-lib", type=float, default=0.815)
    p.add_argument("--lib-min", type=float, default=350.0)
    p.add_argument("--laser", type=float, default=632.8)
    p.add_argument("--temp", type=float, default=80.0)
    p.add_argument("--fwhm", type=float, default=25.0)
    p.add_argument("-o", "--out")
    a = p.parse_args()

    w = {int(k): float(v) for k, v in (x.split(":") for x in a.weights.split(","))}
    models = {}
    for r in csv.DictReader(open(os.path.join(a.jobs, "jobs.csv"))):
        c = int(r["id"])
        if c in w:
            d = parse_raman(os.path.join(a.jobs, r["directory"], "freq", r["input"] + ".out"))
            ok = ((d["freq"] > 1) & (d["freq"] < 1300)) | (d["freq"] > 2500)
            f = d["freq"][ok]
            models[c] = (f, d["I_tot"][ok] * exp_factor(f, a.laser, a.temp) / d["n_H2O"])
    es, el = read_xy(a.exp_stretch), read_xy(a.exp_lattice)
    ys = [spectrum(models, w, a.fwhm), spectrum(models, w, a.fwhm, a.s_oh, a.s_lib, 1.0, a.lib_min)]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 12, "axes.linewidth": 1.1, "xtick.direction": "in",
                         "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), dpi=300, gridspec_kw={"width_ratios": [1, 1.4],
                                                                        "hspace": 0.12, "wspace": 0.06})
    for i, y in enumerate(ys):
        for j, ((lo, hi), exp) in enumerate((((3000, 3600), es), ((150, 1100), el))):
            x = ax[i, j]
            g = (GRID >= lo) & (GRID <= hi)
            se = (exp[:, 0] >= lo) & (exp[:, 0] <= hi)
            x.plot(exp[se, 0], norm(exp[se, 1]), color="#c0392b", lw=1.6, label=a.exp_label)
            x.plot(GRID[g], norm(y[g]), color="black", lw=1.6, label=a.calc_label)
            x.set_xlim(hi, lo)
            x.set_ylim(0, 1.12)
            x.set_yticks([])
            x.minorticks_on()
            if i == 0:
                x.tick_params(labelbottom=False)
            else:
                x.set_xlabel(r"Raman shift (cm$^{-1}$)")
    ax[0, 0].set_ylabel("harmonic", fontsize=13)
    ax[1, 0].set_ylabel("scaled", fontsize=13)
    ax[0, 0].set_title("O–H stretch", fontsize=13)
    ax[0, 1].set_title("lattice modes", fontsize=13)
    ax[1, 0].text(0.04, 0.93, f"×{a.s_oh:.3f}", transform=ax[1, 0].transAxes, fontsize=11, va="top")
    ax[1, 1].text(0.03, 0.93, f"librations ×{a.s_lib:.3f}", transform=ax[1, 1].transAxes, fontsize=11, va="top")
    ax[0, 1].legend(loc="upper left", frameon=False, fontsize=11)
    fig.suptitle(a.phase, fontsize=14, fontweight="bold", y=0.97)
    out = a.out or os.path.join(a.jobs, "analysis", "raman_clean_harmonic_vs_scaled.png")
    fig.savefig(out, bbox_inches="tight")
    print("written", out)


if __name__ == "__main__":
    main()
