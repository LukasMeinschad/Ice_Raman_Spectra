#!/usr/bin/env python3
"""Calculated Raman spectrum of one ordered structure against experimental band positions (2 x 2 panels):
top unscaled, bottom with empirical scale factors (stretch x s_OH for modes > 2500 cm-1, librations
x s_lib for lib_min..1300 cm-1, translations unscaled).  Experimental bands from a reference csv
(H2O_cm-1,fwhm_cm-1,region,intensity,assignment,label); if a FWHM is given, the experimental band is
also drawn as a Lorentzian of that width (height from the intensity class), otherwise as a line.
Prints the calculated modes next to the nearest experimental band.

    python plot_raman_vs_bands.py ../Ice_XI/pbe0_d3_raman --ref ../reference/shigenari_abe_2012_iceXI_raman.csv \
        --s-oh 0.9615 --s-lib 0.815 --laser 532 --temp 50
"""
import argparse
import csv
import os

import numpy as np

from analyze_raman import exp_factor, lorentz, parse_raman

GRID = np.arange(100.0, 3800.0, 0.5)
HEIGHT = {"strongest": 1.0, "vs": 1.0, "s": 0.7, "m": 0.45, "mw": 0.3, "w": 0.25, "vw": 0.12}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("jobs")
    p.add_argument("--id", type=int, default=None, help="model id (default: the first finished one)")
    p.add_argument("--ref", required=True)
    p.add_argument("--ref-label", default="Shigenari & Abe 2012, H$_2$O ice XI, ~50 K")
    p.add_argument("--s-oh", type=float, default=0.9615)
    p.add_argument("--s-lib", type=float, default=0.815)
    p.add_argument("--lib-min", type=float, default=450.0)
    p.add_argument("--laser", type=float, default=514.5)
    p.add_argument("--temp", type=float, default=50.0)
    p.add_argument("--fwhm", type=float, default=15.0)
    p.add_argument("--title", default="Ice XI (Cmc2$_1$): PBE0-D3/POB-TZVP harmonic Raman vs. experiment")
    p.add_argument("-o", "--out")
    a = p.parse_args()

    rows = list(csv.DictReader(open(os.path.join(a.jobs, "jobs.csv"))))
    r = next(x for x in rows if a.id is None or int(x["id"]) == a.id)
    d = parse_raman(os.path.join(a.jobs, r["directory"], "freq", r["input"] + ".out"))
    ok = d["freq"] > 1
    f, I = d["freq"][ok], d["I_tot"][ok] * exp_factor(d["freq"][ok], a.laser, a.temp) / d["n_H2O"]
    irr = [d["irrep"][k] for k in np.where(ok)[0]]
    scale = lambda x: np.where(x > 2500, x * a.s_oh, np.where(x >= a.lib_min, x * a.s_lib, x))
    ref = [x for x in csv.reader(l for l in open(a.ref) if not l.startswith("#"))][1:]
    ref = [(float(x[0]), float(x[1]) if x[1] else None, x[2], x[3], x[4], x[5]) for x in ref]

    print(f"model {r['id']} ({r['space_group']}); scale factors: stretch {a.s_oh}, librations {a.s_lib} (> {a.lib_min:g})")
    print(f"{'calc':>8} {'irrep':5} {'rel. I':>7} {'scaled':>8} | nearest exp band")
    imax = {"st": I[f > 2500].max(), "lat": I[f < 2500].max()}
    for fi, ii, ir, fs in zip(f, I, irr, scale(f)):
        rel = ii / imax["st" if fi > 2500 else "lat"]
        if rel < 0.03 or 1300 < fi < 2500:
            continue
        cand = [x for x in ref if (x[2] == "stretch") == (fi > 2500)]
        nb = min(cand, key=lambda x: abs(x[0] - fs)) if cand else None
        txt = f"{nb[0]:.0f} ({nb[3]}, {nb[4]})  diff {fs - nb[0]:+.0f}" if nb and abs(nb[0] - fs) < 60 else "-"
        print(f"{fi:8.1f} {ir:5} {rel:7.2f} {fs:8.1f} | {txt}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.linewidth": 1.0, "xtick.direction": "in",
                         "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), dpi=200, gridspec_kw={"width_ratios": [1, 1.3]})
    regions = [((2950, 3650), "stretch"), ((450, 1250), "libration")]
    for i, (fs, lab) in enumerate(((f, "harmonic, unscaled"),
                                   (scale(f), f"empirically scaled (factors from ice XIV)\nstretch ×{a.s_oh:.3f}\n"
                                              f"librations ×{a.s_lib:.3f}\ntranslations unscaled"))):
        y = lorentz(fs, I, GRID, a.fwhm)
        for j, ((lo, hi), reg) in enumerate(regions):
            x = ax[i, j]
            g = (GRID >= lo) & (GRID <= hi)
            yexp = np.zeros_like(GRID)
            for nu, w, rg, inten, asg, lb in ref:
                if rg != reg:
                    continue
                h = HEIGHT.get(inten.split()[0], 0.3)
                x.plot([nu, nu], [0, h], color="#c0392b", lw=2.2, solid_capstyle="butt")
                if w:   # experimental FWHM as a horizontal bar at half height
                    x.plot([nu - w / 2, nu + w / 2], [h / 2, h / 2], color="#c0392b", lw=1.0)
                x.axvline(nu, color="#c0392b", ls=":", lw=0.7, zorder=0)
                if i == 0:
                    x.text(nu, 1.02, lb, transform=x.get_xaxis_transform(), rotation=90, fontsize=7,
                           ha="center", va="bottom", color="#7b241c")
            if yexp.any():
                x.plot(GRID[g], yexp[g], color="#c0392b", lw=1.5)
            x.plot(GRID[g], y[g] / y[g].max(), color="black", lw=1.5)
            x.set_xlim(hi, lo)
            x.set_ylim(0, 1.1)
            x.minorticks_on()
            if i == 1:
                x.set_xlabel(r"Raman shift (cm$^{-1}$)")
        ax[i, 0].set_ylabel("Raman intensity (normalised)")
        ax[i, 0].text(0.02, 0.95, f"({'ab'[i]}) {lab}", transform=ax[i, 0].transAxes, fontsize=9, va="top",
                      fontweight="bold")
    ax[0, 0].set_title("O-H stretch region", fontsize=10, pad=24)
    ax[0, 1].set_title("librations (translations below 450 cm$^{-1}$ not in this paper)", fontsize=10, pad=24)
    h = [plt.Line2D([], [], color="black", lw=1.5, label=f"calculated, ice XI, FWHM {a.fwhm:g} cm$^{{-1}}$"),
         plt.Line2D([], [], color="#c0392b", lw=2.2, label=a.ref_label + "\nband positions (stretch: Table X fit, bar = FWHM;\n"
                    "librations: Table IV); heights only show the\nintensity class given in the paper (qualitative)")]
    ax[0, 1].legend(handles=h, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0, frameon=True,
                    edgecolor="black", fancybox=False, fontsize=8.5)
    fig.suptitle(f"{a.title}\ncalculated intensities at {a.temp:g} K, {a.laser:g} nm (powder average)", fontsize=11)
    fig.tight_layout()
    out = a.out or os.path.join(a.jobs, "analysis", "raman_vs_shigenari_abe_2012.png")
    fig.savefig(out, bbox_inches="tight")
    print("written", out)


if __name__ == "__main__":
    main()
