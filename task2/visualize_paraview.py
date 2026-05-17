#!/usr/bin/env python3
"""
Twisted Bifilar Solid — Visualisation Script
============================================
Produces two outputs:

1. ParaView programmatic state (requires pvpython):
     pvpython visualize_paraview.py [mesh.vtu]
   Saves: bifilar_perspective.png, bifilar_topdown.png

2. Matplotlib cross-section plots (always available):
     python3 visualize_paraview.py --matplotlib-only [surface.stl]
   Saves: bifilar_crosssections.png
"""

import sys
import os
import math
import argparse

# ────────────────────────────────────────────────────────────────────────────
# CLI args
# ────────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Bifilar solid visualisation")
parser.add_argument("mesh", nargs="?", default="output/task2/bifilar_tet_mesh.vtu",
                    help="Path to volume mesh VTU")
parser.add_argument("--surface", default="output/task2/bifilar_surface.stl",
                    help="Path to surface STL for cross-section plots")
parser.add_argument("--output-dir", default="output/task2",
                    help="Directory for output images")
parser.add_argument("--matplotlib-only", action="store_true",
                    help="Skip ParaView, only produce matplotlib plots")
args, _ = parser.parse_known_args()

os.makedirs(args.output_dir, exist_ok=True)

# ────────────────────────────────────────────────────────────────────────────
# 1. Matplotlib cross-section plot (no ParaView required)
# ────────────────────────────────────────────────────────────────────────────
def plot_crosssections():
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError:
        print("matplotlib/numpy not available, skipping cross-section plot.")
        return

    # Bifilar geometry parameters (must match task2_pipeline defaults)
    R         = 1.0
    r         = 0.44
    offset    = 0.48
    H         = 4.0
    revolutions = 2.0
    omega     = 2.0 * math.pi * revolutions / H

    n_z = 6  # number of cross-section slices to show
    z_vals = [i * H / (n_z - 1) for i in range(n_z)]

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.suptitle("Twisted Bifilar Solid — Cross-sections at different heights",
                 fontsize=14, fontweight="bold")
    fig.patch.set_facecolor("#1a1a2e")

    for ax, z in zip(axes.flat, z_vals):
        theta = omega * z

        ax.set_facecolor("#16213e")
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.set_aspect("equal")
        ax.set_title(f"z = {z:.2f}  (θ = {math.degrees(theta):.0f}°)",
                     color="white", fontsize=9)
        ax.tick_params(colors="gray")
        for spine in ax.spines.values():
            spine.set_edgecolor("#444")

        # Outer cylinder (blue outline)
        outer = plt.Circle((0, 0), R, fill=False, edgecolor="#4fc3f7",
                           linewidth=2, linestyle="--")
        ax.add_patch(outer)

        # Gap indicator line
        gap = R - offset - r
        ax.annotate("", xy=(R, 0), xytext=(offset + r, 0),
                    arrowprops=dict(arrowstyle="<->", color="#80cbc4", lw=1))
        ax.text((R + offset + r) / 2, 0.04, f"gap={gap:.2f}R",
                ha="center", color="#80cbc4", fontsize=7)

        # Inner cylinder 1 (green/teal)
        cx1 =  offset * math.cos(theta)
        cy1 =  offset * math.sin(theta)
        inner1 = plt.Circle((cx1, cy1), r, color="#26a69a", alpha=0.7)
        ax.add_patch(inner1)

        # Inner cylinder 2 (purple)
        cx2 = -offset * math.cos(theta)
        cy2 = -offset * math.sin(theta)
        inner2 = plt.Circle((cx2, cy2), r, color="#7e57c2", alpha=0.7)
        ax.add_patch(inner2)

        # Centre markers
        ax.plot([cx1, cx2], [cy1, cy2], "+", color="white", ms=5)

        ax.text(-1.2, 1.1, f"r = {r}R", color="#aaa", fontsize=7)

    fig.tight_layout()
    out = os.path.join(args.output_dir, "bifilar_crosssections.png")
    plt.savefig(out, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved: {out}")


# ────────────────────────────────────────────────────────────────────────────
# 2. ParaView programmatic visualisation
# ────────────────────────────────────────────────────────────────────────────
def paraview_visualize():
    try:
        from paraview.simple import (
            OpenDataFile, GetActiveViewOrCreate, Show, Hide,
            ColorBy, GetColorTransferFunction, GetOpacityTransferFunction,
            Threshold, Clip, SaveScreenshot, ResetCamera,
            GetActiveCamera, RenderAllViews, SetActiveSource,
            GetDisplayProperties, CreateRenderView
        )
        import paraview.simple as pv
    except ImportError:
        print("ParaView Python API not available. "
              "Run with: pvpython visualize_paraview.py")
        return

    mesh_file = args.mesh
    if not os.path.exists(mesh_file):
        print(f"Mesh file not found: {mesh_file}")
        return

    print(f"Loading: {mesh_file}")
    reader = OpenDataFile(mesh_file)

    # ── Render view ──────────────────────────────────────────────────────────
    view = GetActiveViewOrCreate("RenderView")
    view.Background = [0.15, 0.15, 0.2]
    view.ViewSize   = [1920, 1080]

    # ── Show full mesh coloured by cell_type ─────────────────────────────────
    display = Show(reader, view)
    ColorBy(display, ("CELLS", "cell_type"))
    ctf = GetColorTransferFunction("cell_type")
    ctf.RGBPoints = [0.0, 0.15, 0.65, 0.60,   # tets  → teal
                     1.0, 0.49, 0.34, 0.76]    # prisms → purple
    ctf.ColorSpace = "RGB"

    display.Opacity = 0.85
    display.PointSize = 1.0

    # ── Perspective view ─────────────────────────────────────────────────────
    ResetCamera()
    cam = GetActiveCamera()
    cam.SetPosition(4.0, -3.5, 3.0)
    cam.SetFocalPoint(0.0, 0.0, 2.0)
    cam.SetViewUp(0, 0, 1)
    RenderAllViews()

    persp_out = os.path.join(args.output_dir, "bifilar_perspective.png")
    SaveScreenshot(persp_out, view, ImageResolution=[1920, 1080])
    print(f"Saved: {persp_out}")

    # ── Top-down view ─────────────────────────────────────────────────────────
    cam.SetPosition(0.0, 0.0, 8.0)
    cam.SetFocalPoint(0.0, 0.0, 2.0)
    cam.SetViewUp(0, 1, 0)
    view.ViewSize = [1080, 1080]
    RenderAllViews()

    top_out = os.path.join(args.output_dir, "bifilar_topdown.png")
    SaveScreenshot(top_out, view, ImageResolution=[1080, 1080])
    print(f"Saved: {top_out}")

    # ── BL layer visualisation ────────────────────────────────────────────────
    prism_filter = Threshold(registrationName="Prisms", Input=reader)
    prism_filter.Scalars = ["CELLS", "cell_type"]
    prism_filter.LowerThreshold = 0.5
    prism_filter.UpperThreshold = 1.5

    SetActiveSource(prism_filter)
    disp2 = Show(prism_filter, view)
    Hide(reader, view)
    ColorBy(disp2, ("CELLS", "layer"))
    GetColorTransferFunction("layer").ApplyPreset("Cool to Warm", True)

    cam.SetPosition(3.0, -2.5, 2.5)
    cam.SetFocalPoint(0.0, 0.0, 2.0)
    view.ViewSize = [1920, 1080]
    RenderAllViews()

    bl_out = os.path.join(args.output_dir, "bifilar_bl_layers.png")
    SaveScreenshot(bl_out, view, ImageResolution=[1920, 1080])
    print(f"Saved: {bl_out}")


# ────────────────────────────────────────────────────────────────────────────
# Run
# ────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    plot_crosssections()
    if not args.matplotlib_only:
        paraview_visualize()
