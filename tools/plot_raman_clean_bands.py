#!/usr/bin/env python3
"""Clean presentation figure for one ordered structure: calculated Raman spectrum vs. experimental band
positions (markers), top harmonic, bottom scaled (stretch x s_OH, librations x s_lib).

    python plot_raman_clean_bands.py ../Ice_XI/pbe0_d3_raman --ref ../reference/shigenari_abe_2012_iceXI_raman.csv
"""
import argparse
import csv
import os

import numpy as np

from analyze_raman import exp_factor, lorentz, parse_raman

GRID = np.arange(100.0, 3800.0, 0.5)
HEIGHT = {"strongest": 1.0, "vs": 1.0, "s": 0.7, "m": 0.45, "mw": 0.3, "w": 0.25, "vw": 0.12}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("jobs")
    p.add_argument("--ref", required=True)
    p.add_argument("--exp-label", default="Experiment (Shigenari & Abe 2012)")
    p.add_argument("--calc-label", default="PBE0-D3")
    p.add_argument("--phase", default="Ice XI")
    p.add_argument("--s-oh", type=float, default=0.9615)
    p.add_argument("--s-lib", type=float, default=0.815)
    p.add_argument("--laser", type=float, default=514.5)
    p.add_argument("--temp", type=float, default=50.0)
    p.add_argument("--fwhm", type=float, default=15.0)
    p.add_argument("-o", "--out")
    a = p.parse_args()
    r = next(csv.DictReader(open(os.path.join(a.jobs, "jobs.csv"))))
    d = parse_raman(os.path.join(a.jobs, r["directory"], "freq", r["input"] + ".out"))
    ok = d["freq"] > 1
    f, I = d["freq"][ok], d["I_tot"][ok] * exp_factor(d["freq"][ok], a.laser, a.temp)
    ref = [x for x in csv.reader(l for l in open(a.ref) if not l.startswith("#"))][1:]
    ref = [(float(x[0]), x[2], x[3]) for x in ref]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 12, "axes.linewidth": 1.1, "xtick.direction": "in",
                         "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), dpi=300, gridspec_kw={"width_ratios": [1, 1.2],
                                                                        "hspace": 0.12, "wspace": 0.06})
    for i, sc in enumerate(((1.0, 1.0), (a.s_oh, a.s_lib))):
        fs = np.where(f > 2500, f * sc[0], np.where(f > 450, f * sc[1], f))
        y = lorentz(fs, I, GRID, a.fwhm)
        for j, ((lo, hi), reg) in enumerate((((2950, 3550), "stretch"), ((450, 1200), "libration"))):
            x = ax[i, j]
            g = (GRID >= lo) & (GRID <= hi)
            for k, (nu, rg, inten) in enumerate([q for q in ref if q[1] == reg]):
                h = HEIGHT.get(inten.split()[0], 0.3)
                x.plot([nu, nu], [0, h], color="#c0392b", lw=2.4, solid_capstyle="butt",
                       label=a.exp_label if k == 0 else None)
            x.plot(GRID[g], y[g] / y[g].max(), color="black", lw=1.6, label=a.calc_label)
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
    ax[0, 1].set_title("librations", fontsize=13)
    ax[1, 0].text(0.04, 0.93, f"×{a.s_oh:.3f}", transform=ax[1, 0].transAxes, fontsize=11, va="top")
    ax[1, 1].text(0.03, 0.93, f"×{a.s_lib:.3f}", transform=ax[1, 1].transAxes, fontsize=11, va="top")
    ax[0, 1].legend(loc="upper left", frameon=False, fontsize=11)
    fig.suptitle(a.phase, fontsize=14, fontweight="bold", y=0.97)
    out = a.out or os.path.join(a.jobs, "analysis", "raman_clean_harmonic_vs_scaled.png")
    fig.savefig(out, bbox_inches="tight")
    print("written", out)


if __name__ == "__main__":
    main()
