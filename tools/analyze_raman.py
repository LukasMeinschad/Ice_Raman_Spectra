#!/usr/bin/env python3
"""Raman spectra from CRYSTAL FREQCALC runs with INTRAMAN (powder / polycrystalline intensities).

Reads every finished <jobs>/config_XXXXX/freq/<name>.out, takes the POLYCRYSTALLINE ISOTROPIC
INTENSITIES (I_tot, I_par, I_perp), optionally converts them to experimental conditions,

    I_exp ~ I * (nu_L - nu)^4 / nu * [1 + n(nu, T)],   n = 1/(exp(hc nu/kT) - 1)

(nu_L = laser wavenumber), broadens them and writes per-model and weighted-average spectra,
plots them in the O-H stretch and the lattice region and marks reference band positions.

    python analyze_raman.py ../Ice_XIV/pbe0_d3_raman --laser 632.8 --temp 80 \
        --ref ../reference/salzmann_2006_iceXIV_raman.csv --weights 1:0.36,2:0.48,3:0.16
"""
import argparse
import csv
import glob
import os
import re

import numpy as np

HC_K = 1.438777  # cm K


def parse_raman(path):
    txt = open(path, errors="replace").read()
    if "EEEEEEEEEE TERMINATION" not in txt or "POLYCRYSTALLINE ISOTROPIC INTENSITIES" not in txt:
        return None
    blk = txt[txt.index("POLYCRYSTALLINE ISOTROPIC INTENSITIES"):]
    blk = blk[:blk.index("SINGLE CRYSTAL")] if "SINGLE CRYSTAL" in blk else blk
    rows = re.findall(r"^\s*(\d+)-\s*(\d+)\s+([\d.]+)\s+\((\S+)\s*\)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", blk, re.M)
    fr = [float(x) for x in re.findall(r"^\s*\d+-\s*\d+\s+\S+\s+(-?[\d.]+)\s+-?[\d.]+\s+\(", txt, re.M)]
    nat = int(re.findall(r"N\. OF ATOMS PER CELL\s+(\d+)", txt)[-1])
    return {"freq": np.array([float(r[2]) for r in rows]), "irrep": [r[3] for r in rows],
            "I_tot": np.array([float(r[4]) for r in rows]), "I_par": np.array([float(r[5]) for r in rows]),
            "I_perp": np.array([float(r[6]) for r in rows]), "n_imag": sum(f < -1 for f in fr), "n_H2O": nat // 3}


def exp_factor(nu, laser_nm, T):
    if not laser_nm:
        return np.ones_like(nu)
    nuL = 1e7 / laser_nm
    bose = 1.0 / (np.exp(HC_K * nu / T) - 1.0) if T else 0.0
    return (nuL - nu) ** 4 / nu * (1.0 + bose)


def lorentz(f, I, grid, fwhm):
    g = fwhm / 2
    return (I[:, None] * (g / np.pi) / ((grid[None, :] - f[:, None]) ** 2 + g ** 2)).sum(0)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("jobs")
    p.add_argument("--laser", type=float, default=None, help="laser wavelength [nm] for the intensity factor")
    p.add_argument("--temp", type=float, default=None, help="temperature [K] for the Bose factor")
    p.add_argument("--fwhm", type=float, default=10.0, help="Lorentzian FWHM [cm-1]")
    p.add_argument("--scale-stretch", type=float, default=1.0, help="scale factor for modes > 2500 cm-1 (plot only)")
    p.add_argument("--weights", default=None, help="id:weight,... (default: multiplicity from jobs.csv)")
    p.add_argument("--ref", help="reference band list (csv with H2O_cm-1, assignment, label)")
    p.add_argument("--title", default="Raman spectrum")
    p.add_argument("--skip-assign", default="b", help="reference assignments not drawn (b = decoupled HOD stretches)")
    p.add_argument("--out-name", default="raman_spectra.png")
    a = p.parse_args()

    jobs = list(csv.DictReader(open(os.path.join(a.jobs, "jobs.csv"))))
    w = {int(r["id"]): float(r["multiplicity"]) for r in jobs}
    if a.weights:
        w = {int(k): float(v) for k, v in (x.split(":") for x in a.weights.split(","))}
    grid = np.arange(0.0, 4000.0, 0.5)
    res = {}
    for r in jobs:
        d = parse_raman(os.path.join(a.jobs, r["directory"], "freq", r["input"] + ".out"))
        if d is None:
            continue
        cid = int(r["id"])
        ok = d["freq"] > 1
        f = d["freq"][ok]
        I = d["I_tot"][ok] * exp_factor(f, a.laser, a.temp)
        I = I / d["n_H2O"]
        fs = np.where(f > 2500, f * a.scale_stretch, f)
        res[cid] = (r["space_group"], f, I, lorentz(fs, I, grid, a.fwhm), d)
        print(f"model {cid} ({r['space_group']}): {len(f)} Raman modes, {d['n_imag']} imaginary; strongest O-H "
              f"stretch {f[f > 2500][np.argmax(I[f > 2500])]:.1f} cm-1")
    if not res:
        raise SystemExit("no finished Raman calculations")
    out = os.path.join(a.jobs, "analysis")
    os.makedirs(out, exist_ok=True)
    ws = np.array([w.get(c, 0.0) for c in res])
    avg = sum(wc * res[c][3] for wc, c in zip(ws, res)) / ws.sum() if ws.sum() > 0 else None
    with open(os.path.join(out, "raman_modes.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["id", "freq_cm-1", "irrep", "I_tot_raw", "I_par_raw", "I_perp_raw", "I_exp_conditions_per_H2O"])
        for c, (sg, f, I, spec, d) in res.items():
            for k in range(len(d["freq"])):
                if d["freq"][k] > 1:
                    j = list(f).index(d["freq"][k])
                    wr.writerow([c, d["freq"][k], d["irrep"][k], d["I_tot"][k], d["I_par"][k], d["I_perp"][k],
                                 round(float(I[j]), 6)])

    ref = []
    if a.ref:
        for l in open(a.ref):
            if l.startswith("#") or l.startswith("H2O") or not l.strip():
                continue
            x = l.rstrip("\n").split(",")
            if x[0] and x[2].strip() not in a.skip_assign.split(","):
                ref.append((float(x[0]), x[2], x[3]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.linewidth": 1.0, "xtick.direction": "in",
                         "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), dpi=200, gridspec_kw={"width_ratios": [1, 1.25]})
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4"]
    regions = [(3800, 3000), (1100, 150)]
    for ax, (hi, lo) in zip(axes, regions):
        s = (grid >= lo) & (grid <= hi)
        curves = [(f"model {c} ({res[c][0]})", res[c][3], cols[k % 5], 1.2) for k, c in enumerate(res)]
        if avg is not None and len(res) > 1:
            curves.append(("weighted average", avg, "black", 1.8))
        top = max(cv[s].max() for _, cv, _, _ in curves)
        for lab, cv, col, lw in curves:
            ax.plot(grid[s], cv[s] / top, color=col, lw=lw, label=lab)
        for nu, asg, lab in ref:
            if lo <= nu <= hi:
                ax.axvline(nu, color="#c0392b", ls=":", lw=0.9, zorder=0)
                ax.text(nu, 1.02, lab.split(" ")[0], transform=ax.get_xaxis_transform(), rotation=90,
                        fontsize=7, ha="center", va="bottom", color="#7b241c")
        ax.set_xlim(hi, lo)
        ax.set_ylim(0, 1.08)
        ax.minorticks_on()
        ax.set_xlabel(r"Raman shift (cm$^{-1}$)")
    axes[0].set_ylabel("Raman intensity (normalised per region)")
    axes[0].set_title("O-H stretch region" + (f" (stretch ×{a.scale_stretch:g})" if a.scale_stretch != 1 else ""),
                      fontsize=10, pad=24)
    axes[1].set_title("lattice region (librations, translations)", fontsize=10, pad=24)
    axes[1].plot([], [], color="#c0392b", ls=":", label="experiment (reference)")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0, frameon=True, edgecolor="black",
                   fancybox=False, fontsize=8.5)
    cond = (f"intensities at {a.temp:g} K, {a.laser:g} nm" if a.laser else "raw CRYSTAL intensities") + \
        f"; Lorentzian FWHM {a.fwhm:g} cm$^{{-1}}$"
    fig.suptitle(f"{a.title} ({cond})", fontsize=11)
    fig.tight_layout()
    path = os.path.join(out, a.out_name)
    fig.savefig(path, bbox_inches="tight")
    print("written", path)


if __name__ == "__main__":
    main()
