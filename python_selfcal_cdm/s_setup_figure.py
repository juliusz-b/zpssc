"""s_setup_figure.py - publication-quality experimental-setup schematic of the
self-calibrating CDM FBG-array interrogator (project PN/01/0321/2022, WP4).

Draws, in one figure:
  * the electronic drive chain  : DC current source (PRO800 ITC) + coded RF from
    the AWG/STM32 -> bias-T -> directly-modulated HCG-VCSEL (BW10-1550-T-T7),
    plus the slow HCG tuning-voltage sweep input;
  * the optical path            : VCSEL -> (optional 90/10 tap) -> circulator
    (port 1->2) -> series chain of 3 identical FBGS DTG-A3A4 gratings (all
    ~1545 nm, separated only by code/delay) -> reflections back to the
    circulator (2->3) -> fiber-coupled InGaAs APD;
  * the Peltier/PID detail       : each grating clamped to an Al block on a TEC
    (TEC1-12706), PT100 + MAX31865 sensing, PID via a DRV8871 H-bridge (~10 pm/C);
    detailed once, referenced for the other two;
  * the processing block         : APD -> NUCLEO-H753ZI (internal ADC ~10 MSPS,
    equivalent-time sampling) -> despread + reference calibration + wavelength
    extraction -> LoRa telemetry;
  * the OPTIONAL MZI k-clock branch (dashed): 90/10 coupler tap -> unbalanced
    Mach-Zehnder (delta-L) -> second detector (PD_ref) -> ADC -> code-gated
    fringe ruler used to linearize the optical-frequency axis.

Optical fiber links are thick teal lines, electronic links thin dark arrows,
the optional branch is dashed, thermal control is dashed red. Output PNG:
figs/fig_setup_experiment.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------
# palette (matches the other schematic scripts in this folder)
# ----------------------------------------------------------------------------
INK      = "#1f4e79"     # box edge / labels
C_ELEC   = "#eaf1fb"     # electronic-block fill (light blue)
C_OPT    = "#fdf3e7"     # optical-block fill (warm)
C_PROC   = "#eafaf0"     # processing fill (light green)
C_PELT   = "#fbeaea"     # thermal/Peltier fill (light red)
C_OPTLN  = "#0b7a86"     # optical link colour (teal)
C_ELECLN = "#333333"     # electronic link colour
C_GREY   = "#8a8a8a"     # optional / auxiliary
C_THERM  = "#b03a3a"     # thermal control colour

OPT_LW   = 3.75           # optical line width
EL_LW    = 1.75           # electronic arrow width

fig, ax = plt.subplots(figsize=(13, 6.5))
ax.set_xlim(0, 25)
ax.set_ylim(0, 12.6)
ax.axis("off")


def box(x, y, w, h, text, fc=C_ELEC, ec=INK, lw=1.75, fs=11.2, ls="solid",
        weight="normal", tcol="black"):
    """Rounded box; returns a dict with centre and edge coordinates."""
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.03,rounding_size=0.14",
        fc=fc, ec=ec, lw=lw, linestyle=ls, zorder=3))
    if text:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fs, zorder=4, color=tcol, weight=weight)
    return dict(x=x, y=y, w=w, h=h, cx=x + w / 2, cy=y + h / 2,
                l=x, r=x + w, t=y + h, b=y)


def earrow(p1, p2, color=C_ELECLN, lw=EL_LW, ls="solid", rad=0.0,
           head=14, z=2):
    """Electronic / control signal arrow (thin, headed)."""
    ax.add_patch(FancyArrowPatch(
        p1, p2, arrowstyle="-|>", mutation_scale=head, lw=lw,
        color=color, linestyle=ls,
        connectionstyle=f"arc3,rad={rad}", zorder=z))


def oline(pts, color=C_OPTLN, lw=OPT_LW, ls="solid", z=1, arrow=False,
          head=18):
    """Optical fiber path as a thick poly-line (optional arrow head at end)."""
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.plot(xs, ys, color=color, lw=lw, ls=ls, solid_capstyle="round",
            solid_joinstyle="round", zorder=z)
    if arrow:
        ax.add_patch(FancyArrowPatch(
            pts[-2], pts[-1], arrowstyle="-|>", mutation_scale=head,
            lw=lw, color=color, linestyle=ls, zorder=z))


# title / caption (top, kept clear of all blocks)
ax.text(12.5, 12.45,
        "Experimental setup: self-calibrating CDM interrogation of a "
        "spectrally-overlapping FBG array",
        ha="center", va="top", fontsize=14.7, weight="bold")

# ============================================================================
# ELECTRONIC DRIVE CHAIN  (top-left)
# ============================================================================
b_idc = box(0.4, 9.7, 3.5, 1.45,
            "DC current source\n(PRO800 ITC)\n$I_\\mathrm{op}\\approx17$ mA",
            fc=C_ELEC, fs=10.9)
b_awg = box(0.4, 7.6, 3.5, 1.45,
            "Code generator\nAWG / STM32\nGold / Kasami $c_i(t)$",
            fc=C_ELEC, fs=10.9)
b_bt = box(5.2, 8.55, 2.1, 1.6, "bias-T\n(DC + RF)", fc=C_ELEC, fs=11.2)

# DC bias and RF code into the bias-T
earrow((b_idc["r"], b_idc["cy"]), (b_bt["l"], b_bt["cy"] + 0.45))
ax.text(4.55, b_bt["cy"] + 0.72, "DC", fontsize=9.8, color=C_ELECLN,
        ha="center", va="bottom")
earrow((b_awg["r"], b_awg["cy"]), (b_bt["l"], b_bt["cy"] - 0.45))
ax.text(4.55, b_bt["cy"] - 0.72, "RF", fontsize=9.8, color=C_ELECLN,
        ha="center", va="top")

# ============================================================================
# SOURCE: VCSEL
# ============================================================================
b_vcsel = box(8.6, 8.35, 3.6, 2.0,
              "Tunable HCG-VCSEL\nBW10-1550-T-T7\n1550 nm, integ. TEC\n"
              "(direct intensity mod.)",
              fc=C_OPT, weight="bold", fs=11.2)

# bias-T output -> VCSEL drive (current modulation carries the code)
earrow((b_bt["r"], b_bt["cy"]), (b_vcsel["l"], b_vcsel["cy"] + 0.35))
ax.text((b_bt["r"] + b_vcsel["l"]) / 2, b_vcsel["cy"] + 0.62,
        "$I_\\mathrm{LD}(t)$", fontsize=9.8, color=C_ELECLN,
        ha="center", va="bottom")

# slow HCG tuning sweep (separate input, from below)
b_tune = box(8.6, 6.15, 3.6, 1.4,
             "HCG tuning sweep\n$V_t$: 0..$-$12 V (slow)\n+ ESD clamp / RC",
             fc=C_ELEC, fs=10.6)
earrow((b_tune["cx"], b_tune["t"]), (b_vcsel["cx"], b_vcsel["b"]))
ax.text(b_vcsel["cx"] + 0.15, (b_tune["t"] + b_vcsel["b"]) / 2,
        "$\\lambda$ sweep", fontsize=9.5, color=C_ELECLN, ha="left",
        va="center")

# ============================================================================
# CIRCULATOR
# ============================================================================
b_circ = box(15.1, 8.35, 2.4, 2.0, "", fc=C_OPT)
ax.text(b_circ["cx"], b_circ["t"] - 0.4, "circulator", fontsize=11.2,
        ha="center", va="center", color="black")
ax.text(b_circ["l"] + 0.3, b_circ["cy"] + 0.35, "1", fontsize=10.5,
        color=INK, ha="center", va="center", weight="bold")
ax.text(b_circ["r"] - 0.3, b_circ["cy"] + 0.05, "2", fontsize=10.5,
        color=INK, ha="center", va="center", weight="bold")
ax.text(b_circ["l"] + 0.3, b_circ["b"] + 0.4, "3", fontsize=10.5,
        color=INK, ha="center", va="center", weight="bold")
# curved 1->2->3 hint inside
ax.add_patch(FancyArrowPatch((b_circ["cx"] - 0.55, b_circ["cy"] + 0.3),
                             (b_circ["cx"] + 0.5, b_circ["cy"] - 0.1),
                             arrowstyle="-|>", mutation_scale=10,
                             lw=1.25, color=C_OPTLN,
                             connectionstyle="arc3,rad=-0.5", zorder=4))

# ============================================================================
# OPTIONAL 90/10 TAP  (on VCSEL -> circulator link)
# ============================================================================
tap_x = 13.35
tap_y = b_vcsel["cy"]
# VCSEL -> tap -> circulator port 1  (main optical, solid teal)
oline([(b_vcsel["r"], b_vcsel["cy"]), (tap_x - 0.16, tap_y)])
oline([(tap_x + 0.16, tap_y), (b_circ["l"], b_circ["cy"] + 0.35)],
      arrow=True)
ax.add_patch(Circle((tap_x, tap_y), 0.16, fc="white", ec=C_OPTLN,
                    lw=2.50, zorder=5))
ax.text(tap_x, tap_y + 0.34, "90/10\ncoupler", fontsize=9.0, color=C_OPTLN,
        ha="center", va="bottom")

# ============================================================================
# FBG CHAIN (sensing fiber) + RETURN  (right, horizontal rails)
# ============================================================================
chain_y = 9.55                       # forward sensing rail height
ret_y   = 8.05                       # reflection-return rail height
fbg_x   = [19.6, 21.3, 23.0]
fbg_lab = ["FBG$_1$\nref", "FBG$_2$\nref", "FBG$_3$\nsensor"]
fbg_codes = ["$c_1$", "$c_2$", "$c_3$"]
rail_end = fbg_x[-1] + 0.95

# forward sensing-fiber rail: circ port 2 -> right along the chain
oline([(b_circ["r"], b_circ["cy"] + 0.35), (18.4, b_circ["cy"] + 0.35),
       (18.4, chain_y), (rail_end, chain_y)], lw=OPT_LW)
# return rail: from chain end back down to circulator port 2 (reflections)
oline([(rail_end, chain_y), (rail_end, ret_y), (18.4, ret_y),
       (18.4, b_circ["cy"] - 0.0), (b_circ["r"], b_circ["cy"] - 0.0)],
      lw=OPT_LW, arrow=True)
# reflection direction arrow on the return rail
ax.add_patch(FancyArrowPatch((rail_end - 0.6, ret_y), (19.4, ret_y),
                             arrowstyle="-|>", mutation_scale=17,
                             lw=OPT_LW, color=C_OPTLN, zorder=2))
ax.text((18.4 + rail_end) / 2, ret_y - 0.32, "reflections (sum of channels)",
        fontsize=9.5, color=C_OPTLN, ha="center", va="top")

# the three gratings as small fringed rectangles on the forward rail
for gx, lab, cc in zip(fbg_x, fbg_lab, fbg_codes):
    ax.add_patch(Rectangle((gx - 0.28, chain_y - 0.42), 0.56, 0.84,
                           fc="#f4c27a", ec=INK, lw=1.50, zorder=4))
    for fxoff in np.linspace(-0.18, 0.18, 5):
        ax.plot([gx + fxoff, gx + fxoff],
                [chain_y - 0.36, chain_y + 0.36],
                color=INK, lw=0.88, zorder=5)
    ax.text(gx, chain_y + 0.52, lab, fontsize=9.8, ha="center",
            va="bottom", color="black")
    ax.text(gx, chain_y - 0.52, cc, fontsize=10.1, ha="center", va="top",
            color=INK)

# spacing annotation (~4 m apart) above the chain
ax.annotate("", xy=(fbg_x[1], chain_y + 1.55), xytext=(fbg_x[0], chain_y + 1.55),
            arrowprops=dict(arrowstyle="<->", color="#666", lw=1.12))
ax.text((fbg_x[0] + fbg_x[1]) / 2, chain_y + 1.62, "$\\approx$4 m",
        fontsize=9.5, color="#666", ha="center", va="bottom")

# header note on the chain (kept below the title band)
ax.text((fbg_x[0] + fbg_x[-1]) / 2, chain_y + 2.05,
        "3$\\times$ FBGS DTG-A3A4, all $\\approx$1545 nm, FWHM $\\approx$250 pm, "
        "R = 10%\n(spectrally overlapping; separated only by code / delay)",
        fontsize=9.8, ha="center", va="bottom", color="black")

# ---- Peltier / PID detail (drawn fully under FBG_1, referenced for 2 & 3) ----
b_pel = box(18.55, 5.35, 3.6, 1.9, "", fc=C_PELT)
ax.text(b_pel["cx"], b_pel["t"] - 0.3, "Peltier stage (per grating)",
        fontsize=10.1, ha="center", va="center", weight="bold")
ax.text(b_pel["cx"], b_pel["cy"] - 0.22,
        "Al block + TEC1-12706\nPT100 + MAX31865\nPID, DRV8871 H-bridge\n"
        "$\\approx$10 pm/$^\\circ$C tuning",
        fontsize=9.7, ha="center", va="center")
# thermal coupling dashed up to FBG_1
earrow((b_pel["cx"] - 1.0, b_pel["t"]), (fbg_x[0], chain_y - 0.42),
       color=C_THERM, ls="dashed", lw=1.50, head=12)
ax.text(fbg_x[0] - 1.45, (b_pel["t"] + chain_y) / 2,
        "$T$ control", fontsize=9.2, color=C_THERM, ha="center", va="center")
# faint thermal stubs implying identical stages on FBG2, FBG3
for gx in fbg_x[1:]:
    ax.plot([gx, gx], [chain_y - 0.42, chain_y - 0.95], color=C_THERM,
            lw=1.25, ls="dashed", zorder=2)
ax.text((fbg_x[1] + fbg_x[2]) / 2, chain_y - 1.18,
        "$\\times$3 identical stages\n(one per grating)",
        fontsize=9.1, style="italic", color=C_THERM, ha="center", va="top")

# ============================================================================
# RETURN PATH -> APD -> PROCESSING -> LoRa  (bottom band)
# ============================================================================
b_apd = box(15.6, 0.6, 2.7, 1.7,
            "fiber-coupled\nInGaAs APD\nA-CUBE-I200-100-FC",
            fc=C_OPT, fs=10.4)
# circulator port 3 (bottom) down to APD
oline([(b_circ["cx"] - 0.45, b_circ["b"]), (b_circ["cx"] - 0.45, 1.45),
       (b_apd["r"], 1.45)], arrow=True)
ax.text(b_circ["cx"] - 0.6, (b_circ["b"] + 1.45) / 2 + 1.1,
        "2$\\rightarrow$3", fontsize=9.5, color=C_OPTLN, ha="right",
        va="center")

b_proc = box(8.4, 0.5, 6.1, 1.9,
             "NUCLEO-H753ZI\ninternal ADC $\\approx$10 MSPS (ETS)\n"
             "despread $+$ reference calibration\n"
             "$+$ wavelength extraction",
             fc=C_PROC, fs=10.9)
earrow((b_apd["l"], b_apd["cy"]), (b_proc["r"], b_proc["cy"]))
ax.text((b_apd["l"] + b_proc["r"]) / 2, b_proc["cy"] + 0.26,
        "analog", fontsize=9.2, color=C_ELECLN, ha="center", va="bottom")

# code reference feeding the despreader (left side, long dashed, clear route)
earrow((b_awg["l"] + 0.1, b_awg["b"]), (b_proc["l"] + 0.1, b_proc["t"]),
       color=C_ELECLN, ls=(0, (4, 3)), lw=1.25, rad=-0.30, head=12)
ax.text(4.55, 4.0, "$c_i(t)$ reference\n(despread / gate)", fontsize=9.2,
        color=C_ELECLN, ha="center", va="center", style="italic")

# LoRa telemetry out
b_lora = box(2.2, 0.6, 4.6, 1.0, "LoRa telemetry  ($\\lambda_i \\to$ host)",
             fc=C_ELEC, fs=10.9)
earrow((b_proc["l"], b_proc["cy"]), (b_lora["r"], b_lora["cy"]))

# ============================================================================
# OPTIONAL MZI k-CLOCK BRANCH  (lower-centre band, dashed grey)
# ============================================================================
b_mzi = box(11.0, 3.35, 3.0, 1.6,
            "unbalanced MZI\n($\\Delta L$ k-clock)\nFSR $\\approx$ 10 GHz",
            fc=C_OPT, ec=C_GREY, lw=1.62, ls="dashed", fs=10.4, tcol="#444")
b_pdref = box(14.7, 3.5, 1.9, 1.3, "PD$_\\mathrm{ref}$\n+ ADC",
              fc=C_PROC, ec=C_GREY, lw=1.62, ls="dashed", fs=10.4, tcol="#444")
b_ruler = box(17.2, 3.35, 5.0, 1.6,
              "code-gated fringe ruler\n$\\nu_k=\\nu_0+k\\cdot$FSR\n"
              "(linearizes optical-freq. axis)",
              fc=C_PROC, ec=C_GREY, lw=1.62, ls="dashed", fs=10.4, tcol="#444")

# optical tap (10%) down into MZI  (dashed teal-grey = optional optical)
oline([(tap_x, tap_y - 0.16), (tap_x, 6.7), (11.0 + 0.9, 6.7),
       (11.0 + 0.9, b_mzi["t"])], color=C_GREY, ls="dashed", lw=2.75,
      arrow=True, head=16)
ax.text(tap_x + 0.2, 7.1, "10%", fontsize=9.2, color=C_GREY, ha="left",
        va="center")
# MZI -> PD_ref (optical, dashed) -> ruler (electronic, dashed)
oline([(b_mzi["r"], b_mzi["cy"]), (b_pdref["l"], b_pdref["cy"])],
      color=C_GREY, ls="dashed", lw=2.75, arrow=True, head=16)
earrow((b_pdref["r"], b_pdref["cy"]), (b_ruler["l"], b_ruler["cy"]),
       color=C_GREY, ls="dashed")
# ruler axis info -> processing
earrow((b_ruler["l"], b_ruler["t"] - 0.15), (b_proc["r"] - 0.4, b_proc["b"]),
       color=C_GREY, ls="dashed", lw=1.25, rad=-0.2, head=12)
ax.text(11.0, 5.15, "optional axis-linearization branch", fontsize=9.8,
        style="italic", color=C_GREY, ha="left", va="bottom")

# ============================================================================
# LEGEND  (clear lower-left area)
# ============================================================================
legend_elems = [
    Line2D([0], [0], color=C_OPTLN, lw=OPT_LW, label="optical fiber path"),
    Line2D([0], [0], color=C_ELECLN, lw=EL_LW, marker=">", markersize=8,
           label="electronic signal"),
    Line2D([0], [0], color=C_GREY, lw=2.50, ls="dashed",
           label="optional MZI k-clock branch"),
    Line2D([0], [0], color=C_THERM, lw=1.50, ls="dashed",
           label="thermal (Peltier/PID) control"),
]
leg = ax.legend(handles=legend_elems, loc="lower left",
                bbox_to_anchor=(0.005, 0.18), fontsize=10.4,
                frameon=True, framealpha=0.95, borderpad=0.7,
                handlelength=2.4)
leg.get_frame().set_edgecolor(INK)

plt.tight_layout()
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "fig_setup_experiment.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print("saved", out_path)