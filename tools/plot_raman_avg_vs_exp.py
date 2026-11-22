#!/usr/bin/env python3
"""Weighted-average calculated Raman spectrum against digitised experimental spectra, 2 x 2 panels:

  top    : O-H stretch region | lattice region, harmonic frequencies as calculated
  bottom : the same with empirical frequency scaling per region
             stretch modes (> 2500 cm-1)       x s_OH
             librations   (lib_min ... 1300)   x s_lib
             translations (< lib_min)          unscaled (or x s_trans)
           s_OH and s_lib are fitted (maximum overlap = Pearson correlation of the broadened
           spectra with the experimental curve in each region) unless given explicitly.

Intensities are converted to the experimental conditions with the *unscaled* frequencies.

    python plot_raman_avg_vs_exp.py ../Ice_XIV/pbe0_d3_raman --weights 1:0.437,2:0.323,3:0.240 \
        --exp-stretch ../reference/salzmann_2006_iceXIV_fig1a_H2O.csv \
        --exp-lattice ../reference/salzmann_2006_iceXIV_fig1c_H2O.csv \
        --laser 632.8 --temp 80 --fwhm 25 --ref ../reference/salzmann_2006_iceXIV_raman.csv
"""
import argparse
import csv
import os

import numpy as np

from analyze_raman import exp_factor, lorentz, parse_raman

GRID = np.arange(100.0, 3900.0, 0.5)
STRETCH = (3000.0, 3650.0)
LATTICE = (150.0, 1100.0)


def read_xy(path):
    xy = np.array([[float(v) for v in l.split(",")[:2]] for l in open(path) if l[:1].isdigit()])
    return xy[np.argsort(xy[:, 0])]


def norm(y):
    y = y - y.min()
    return y / y.max()


def spectrum(models, w, fwhm, s_oh=1.0, s_lib=1.0, s_tr=1.0, lib_min=350.0):
    y = np.zeros_like(GRID)
    for c, (f, I) in models.items():
        fs = np.where(f > 2500, f * s_oh, np.where(f >= lib_min, f * s_lib, f * s_tr))
        y += w[c] * lorentz(fs, I, GRID, fwhm)
    return y / sum(w[c] for c in models)


def overlap(y, exp, lo, hi):
    s = (exp[:, 0] >= lo) & (exp[:, 0] <= hi)
    yc = np.interp(exp[s, 0], GRID, y)
    return np.corrcoef(yc, exp[s, 1])[0, 1]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("jobs")
    p.add_argument("--weights", required=True, help="id:weight,...")
    p.add_argument("--exp-stretch", required=True)
    p.add_argument("--exp-lattice", required=True)
    p.add_argument("--exp-label", default="Salzmann et al. 2006, H$_2$O ice XIV, 80 K")
    p.add_argument("--laser", type=float, default=632.8)
    p.add_argument("--temp", type=float, default=80.0)
    p.add_argument("--fwhm", type=float, default=25.0)
    p.add_argument("--s-oh", type=float, help="stretch scale factor (default: fitted)")
    p.add_argument("--s-lib", type=float, help="libration scale factor (default: fitted)")
    p.add_argument("--s-trans", type=float, default=1.0)
    p.add_argument("--lib-min", type=float, default=350.0, help="modes above this (and < 1300) are librations")
    p.add_argument("--ref", help="reference band list (H2O_cm-1,...,assignment,label)")
    p.add_argument("--skip-assign", default="b")
    p.add_argument("--title", default="Ice XIV: calculated (PBE0-D3, weighted average) vs. experiment")
    p.add_argument("--show-models", action="store_true", help="also draw the individual models (thin lines)")
    p.add_argument("-o", "--out", default=None)
    a = p.parse_args()

    w = {int(k): float(v) for k, v in (x.split(":") for x in a.weights.split(","))}
    models = {}
    for r in csv.DictReader(open(os.path.join(a.jobs, "jobs.csv"))):
        c = int(r["id"])
        if c not in w:
            continue
        d = parse_raman(os.path.join(a.jobs, r["directory"], "freq", r["input"] + ".out"))
        if d is None:
            raise SystemExit(f"model {c}: no finished Raman output")
        ok = (d["freq"] > 1) & (d["freq"] < 1300) | (d["freq"] > 2500)
        f = d["freq"][ok]
        models[c] = (f, d["I_tot"][ok] * exp_factor(f, a.laser, a.temp) / d["n_H2O"])
    es, el = read_xy(a.exp_stretch), read_xy(a.exp_lattice)

    s_oh, s_lib = a.s_oh, a.s_lib
    if s_oh is None:
        tr = np.arange(0.90, 1.0001, 0.0005)
        s_oh = tr[np.argmax([overlap(spectrum(models, w, a.fwhm, s_oh=s), es, *STRETCH) for s in tr])]
    if s_lib is None:
        tr = np.arange(0.70, 1.0001, 0.0025)
        lo = a.lib_min * 1.1
        s_lib = tr[np.argmax([overlap(spectrum(models, w, a.fwhm, s_lib=s, s_tr=a.s_trans, lib_min=a.lib_min),
                                      el, lo, LATTICE[1]) for s in tr])]
    y0 = spectrum(models, w, a.fwhm)
    y1 = spectrum(models, w, a.fwhm, s_oh, s_lib, a.s_trans, a.lib_min)
    one = lambda c, *sc: spectrum({c: models[c]}, {c: 1.0}, a.fwhm, *sc)
    ind = [{c: one(c) for c in models}, {c: one(c, s_oh, s_lib, a.s_trans, a.lib_min) for c in models}]
    mcol = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4"]
    print(f"scale factors: stretch {s_oh:.4f}, librations {s_lib:.4f} (> {a.lib_min:g} cm-1), "
          f"translations {a.s_trans:.4f}")
    for tag, y in (("unscaled", y0), ("scaled", y1)):
        print(f"  {tag:9s} overlap r: stretch {overlap(y, es, *STRETCH):.3f}, librations "
              f"{overlap(y, el, a.lib_min * 1.1, LATTICE[1]):.3f}, whole lattice region {overlap(y, el, *LATTICE):.3f}")
    for tag, y in (("unscaled", y0), ("scaled", y1)):
        g = (GRID >= STRETCH[0]) & (GRID <= STRETCH[1])
        print(f"  {tag:9s} stretch maximum {GRID[g][np.argmax(y[g])]:.0f} cm-1 (exp "
              f"{es[np.argmax(es[:, 1]), 0]:.0f})")

    ref = []
    if a.ref:
        for l in open(a.ref):
            x = l.rstrip("\n").split(",")
            if l.startswith("#") or l.startswith("H2O") or not x[0]:
                continue
            if x[2].strip() not in a.skip_assign.split(","):
                ref.append((float(x[0]), x[3]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.linewidth": 1.0, "xtick.direction": "in",
                         "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), dpi=200, gridspec_kw={"width_ratios": [1, 1.3]})
    rows = [(y0, "harmonic, unscaled"),
            (y1, f"empirically scaled\nstretch ×{s_oh:.3f}\nlibrations ×{s_lib:.3f}"
                 + (f"\ntranslations ×{a.s_trans:.3f}" if a.s_trans != 1 else "\ntranslations unscaled"))]
    for i, (y, lab) in enumerate(rows):
        for j, ((lo, hi), exp) in enumerate(((STRETCH, es), (LATTICE, el))):
            x = ax[i, j]
            g = (GRID >= lo) & (GRID <= hi)
            se = (exp[:, 0] >= lo) & (exp[:, 0] <= hi)
            x.plot(exp[se, 0], norm(exp[se, 1]), color="#c0392b", lw=1.5, label=a.exp_label)
            if a.show_models:
                for k, (c, yc) in enumerate(ind[i].items()):
                    x.plot(GRID[g], norm(yc[g]), color=mcol[k % 5], lw=0.9, alpha=0.8,
                           label=f"model {c} (weight {w[c]:.2f})")
            x.plot(GRID[g], norm(y[g]), color="black", lw=1.8, label="calculated, weighted average")
            for nu, lb in ref:
                if lo <= nu <= hi:
                    x.axvline(nu, color="#c0392b", ls=":", lw=0.8, zorder=0)
                    if i == 0:
                        x.text(nu, 1.02, lb.split(" ")[0], transform=x.get_xaxis_transform(), rotation=90,
                               fontsize=7, ha="center", va="bottom", color="#7b241c")
            x.set_xlim(hi, lo)
            x.set_ylim(0, 1.1)
            x.minorticks_on()
            if i == 1:
                x.set_xlabel(r"Raman shift (cm$^{-1}$)")
        ax[i, 0].set_ylabel("Raman intensity (normalised)")
        ax[i, 0].text(0.02, 0.95, f"({'ab'[i]}) {lab}", transform=ax[i, 0].transAxes, fontsize=9,
                      va="top", fontweight="bold")
    ax[0, 0].set_title("O-H stretch region", fontsize=10, pad=24)
    ax[0, 1].set_title("lattice region (librations, translations)", fontsize=10, pad=24)
    ax[0, 1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0, frameon=True,
                    edgecolor="black", fancybox=False, fontsize=8.5)
    fig.suptitle(f"{a.title}\nintensities at {a.temp:g} K, {a.laser:g} nm; Lorentzian FWHM {a.fwhm:g} cm$^{{-1}}$",
                 fontsize=11)
    fig.tight_layout()
    out = a.out or os.path.join(a.jobs, "analysis", "raman_avg_vs_exp_scaled.png")
    fig.savefig(out, bbox_inches="tight")
    print("written", out)


if __name__ == "__main__":
    main()
