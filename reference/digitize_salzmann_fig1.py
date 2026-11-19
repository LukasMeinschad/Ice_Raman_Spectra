#!/usr/bin/env python3
"""Digitise the H2O Raman spectra of ice XIV (Fig. 1a and 1c) of
C.G. Salzmann, A. Hallbrucker, J.L. Finney, E. Mayer, Chem. Phys. Lett. 429, 469 (2006):
recovered H2O ice XIV (9.0 mol% HOD), 80 K, ~10 mbar, 632.8 nm.

  a: coupled O-H stretching region, 3600 -> ~2950 cm-1
  c: low-frequency region (librations, translations), 1000 -> 200 cm-1
Intensity axis: arbitrary units (normalised to the maximum of each panel here).

Calibration from the major tick marks (pixel coordinates of the 1758 x 978 image):
  a: x = 348 -> 3600, 452 -> 3400, 555.5 -> 3200, 659.5 -> 3000 cm-1
  c: x = 1068 -> 1000, 1152 -> 800, 1237 -> 600, 1321 -> 400, 1405 -> 200 cm-1
  frame of the top row: y = 39 (top) ... 376 (bottom)

    python digitize_salzmann_fig1.py salzmann_2006_iceXIV_fig1.png [--check overlay.png]
"""
import argparse
import os

import numpy as np
from PIL import Image, ImageDraw

PANELS = {
    # name: (x_first, x_last, (x_px_ref, nu_ref), cm-1 per px (negative = decreasing to the right))
    "a_OH_stretch": (350, 680, (348.0, 3600.0), -200.0 / 103.83),
    "c_lattice": (1071, 1401, (1068.0, 1000.0), -200.0 / 84.25),
}
YLIM = (44, 360)
# labels, marker lines and the arrow in panel c that touch or approach the curve (x0, x1, y0, y1),
# blanked before tracing (read off the pixel map of the figure)
MASKS = [
    (1288, 1320, 175, 216),   # "488" / "462" labels
    (1303, 1332, 255, 275),   # "413" label
    (1313, 1317, 277, 298),   # 413 marker line
    (1342, 1373, 272, 290),   # "314" label
    (1354, 1359, 293, 313),   # 314 marker line
    (1354, 1384, 146, 163),   # "241" label
    (1387, 1390, 140, 166),   # 241 marker line
    (1372, 1394, 127, 142),   # "214" label
    (1395, 1400, 127, 234),   # 214 arrow
    (1401, 1402, 200, 212),   # 214 arrow head
]


def runs(col, thr):
    """(centre, top) of every dark pixel run in a column"""
    idx = np.where(col < thr)[0]
    if len(idx) == 0:
        return []
    out, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 1:
            cur.append(i)
        else:
            out.append((np.mean(cur), cur[0]))
            cur = [i]
    out.append((np.mean(cur), cur[0]))
    return out


def trace(img, x0, x1, thr=130, win=5):
    """follow the curve with the run centres; record the top of the run (= highest point of the
    curve in that column, so that narrow peaks drawn as vertical runs keep their height)"""
    c0 = [r for r in runs(img[:, x0], thr) if YLIM[0] <= r[0] <= YLIM[1]]
    yl = max(r[0] for r in c0)
    pts, xl, slope = {x0: max(c0)[1]}, x0, 0.0
    centres = {x0: yl}
    for x in range(x0 + 1, x1 + 1):
        ypred = yl + slope * (x - xl)
        w = min(win + 1.5 * (x - xl - 1), 14) + 2.5 * abs(slope)
        cand = [r for r in runs(img[:, x], thr) if YLIM[0] <= r[0] <= YLIM[1] and abs(r[0] - ypred) <= w]
        if not cand:
            # steep flanks of the narrow translational peaks: accept the only curve run in the column
            allr = [r for r in runs(img[:, x], thr) if YLIM[0] <= r[0] <= YLIM[1]]
            if len(allr) != 1:
                continue
            cand = allr
        c, top = min(cand, key=lambda r: abs(r[0] - ypred))
        pts[x] = top
        centres[x] = c
        xs = [k for k in sorted(centres) if k >= x - 6]
        if len(xs) >= 3:
            slope = float(np.polyfit(xs, [centres[k] for k in xs], 1)[0])
        xl, yl = x, c
    return pts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("image")
    p.add_argument("--check")
    a = p.parse_args()
    img = np.array(Image.open(a.image).convert("L")).astype(float)
    for x0, x1, y0, y1 in MASKS:
        img[y0:y1 + 1, x0:x1 + 1] = 255.0
    here = os.path.dirname(os.path.abspath(__file__))
    allpts = []
    for k, (name, (x0, x1, (xr, nr), slope)) in enumerate(PANELS.items()):
        pts = trace(img, x0, x1)
        ys = np.array([pts[x] for x in sorted(pts)])
        base, top = ys.max(), ys.min()           # y grows downwards
        out = os.path.join(here, f"salzmann_2006_iceXIV_fig1{name[0]}_H2O.csv")
        with open(out, "w") as fh:
            fh.write("# Digitised from Salzmann et al., Chem. Phys. Lett. 429, 469 (2006), Fig. 1"
                     f"{name[0]}: H2O ice XIV (9.0 mol% HOD), 80 K, 632.8 nm.\n"
                     "# intensity normalised to the panel maximum (0 = lowest point of the trace)\n")
            fh.write("raman_shift_cm-1,intensity\n")
            for x in sorted(pts):
                fh.write(f"{nr + slope * (x - xr):.1f},{(base - pts[x]) / (base - top):.4f}\n")
        print(f"{name}: {len(pts)} points, {nr + slope * (x0 - xr):.0f} -> {nr + slope * (x1 - xr):.0f} cm-1 -> {out}")
        allpts += [(x, y, k) for x, y in pts.items()]
    if a.check:
        im = Image.open(a.image).convert("RGB")
        W, H = im.size
        im = im.resize((2 * W, 2 * H), Image.NEAREST)
        d = ImageDraw.Draw(im)
        col = [(230, 0, 0), (0, 140, 255)]
        for x, y, k in allpts:
            d.ellipse([2 * x - 1.5, 2 * y - 1.5, 2 * x + 1.5, 2 * y + 1.5], fill=col[k])
        im.crop((2 * 330, 2 * 30, 2 * 1420, 2 * 390)).save(a.check)


if __name__ == "__main__":
    main()
