# python_selfcal_cdm - simulation study of a coded FBG interrogator

Python scripts behind the paper *Limiting Factors for CDM Interrogation of Fiber Bragg Grating Arrays* (J. Bojarczuk, manuscript, 2026). Many fiber Bragg gratings on one fiber share a wavelength band and are told apart by the delay of their echo. A directly modulated, swept VCSEL launches a spreading code, the photodiode sees the sum of the delayed echoes, and correlation at every wavelength step recovers the spectrum of every grating. The scripts quantify what limits the number of gratings in a band and the accuracy of the recovered Bragg wavelengths: spectral shadowing, the critical reflectivity above which a co-tuned return flattens and splits, the order of gratings of unequal width, third-order ghosts, code leakage, source chirp, and beat noise between returns of one wavelength step.

Two models are used. The **power model** (`common.py`, scripts `s*.py`) adds the echoes in power with Gaussian lines and third-order ghosts. The **coherent model** (`s70_field.py`, folder `coherent/`) propagates the optical field through coupled-mode gratings with every order of multiple reflections and a rate-equation laser. Code and comments are in English.

## Requirements

Python 3.10 or newer with `numpy`, `scipy` and `matplotlib`:

```bash
pip install numpy scipy matplotlib
```

## Figures of the paper

| Fig. | Content | Script | Output |
|---|---|---|---|
| 1 | principle of the stepped-sweep CDM interrogation | `s_fig1_principle.py` | `figs/fig_s_fig1_principle` |
| 2 | shadowing shift and the two-grating rule | `s26_laws_concept.py` | `figs/fig_s26_laws_concept` |
| 3 | critical reflectivity, grating order, all 40320 orders (coherent model) | `coherent/fig_validation.py` | `figs/fig_validation_curvature` |
| 4 | third-order ghost paths | `s50_ghostpath.py` | `figs/fig_s50_ghostpath` |
| 5 | ghost collisions and the Golomb ruler | `s53_ruler.py` | `figs/fig_s53_ruler` |
| 6 | error mechanisms in the delay-wavelength map | `s16_principle.py` | `figs/fig_s16_mechanisms` |
| 7 | deshadowing | `s19_deshadow.py` | `figs/fig_s19_deshadow` |
| 8 | source chirp and reference correction | `s18_source.py` | `figs/fig_s18_source` |
| 9 | capacity of one band against reflectivity, ghost-limited curve | `s12_capacity.py` | `figs/fig_s12_capacity` |
| 10 | worst-case bound against both models | `s55_error_bound.py`, `coherent/fig_bound.py` | `figs/fig_bound_check` |
| 11 | shadowing shift and beat-limited SNR in the coherent model | `coherent/fig_validation.py` | `figs/fig_validation_mechanisms` |
| 12 | the 50-sensor array of Markowski et al. (JLT 2023) before and after the design rules | `coherent/fig_jlt.py` | `figs/fig_jlt_style` |

Each script writes PDF and PNG to `figs/`. `s12_capacity.py` recomputes about two hours of layouts, `S12_REPLOT=1 python s12_capacity.py` redraws from the cache `out/s12_results.npz`. The coherent-model figures are drawn from compact caches in `coherent/cache/`, see `coherent/README.md` for the record campaigns behind them.

## Other scripts

Supporting analyses that the paper cites in numbers or that fed the design rules. Each one is self-contained and documents its question in the docstring.

| Script | Topic |
|---|---|
| `s1_codes.py`, `s14_code_families.py`, `s30_beyond_mseq.py`, `s36_long_codes.py` | correlation side lobes, m-sequences against Gold, Kasami and Golay codes for a swept source that transmits one code at a time, leakage ceiling for long codes |
| `s13_chiprate.py`, `s27_ets.py`, `s61_crb_sampling.py` | delay resolution against chip rate, equivalent-time sampling, Cramer-Rao bound for the number of sweep steps |
| `s21_cdmwdm.py`, `s22_codelength.py`, `s23_theory.py` | CDM-WDM banding, code length as a design variable, closed-form rules and their numerical check |
| `s29_ghostpath.py`, `s31_golomb.py`, `s32_golomb_modular.py`, `s33_ncpc_robustness.py`, `s43_orders.py`, `s51_ghostcorr.py`, `s54_modular_ruler.py` | ghost paths, Golomb and modular rulers, reflection orders beyond the third |
| `s41_spectra.py`, `s42_taper.py`, `s44_additivity.py`, `s45_lineshape.py`, `s48_width_spread.py`, `s57_estimator_checks.py` | reconstructed lines, additivity of pairwise shifts, Gaussian line against coupled-mode spectra, unequal widths, the local Gaussian fit |
| `s20_inversion.py`, `s56_peel_depth.py`, `s58_decorrelator.py`, `s59_bidirectional.py`, `s60_two_end.py` | what deshadowing needs to know, its depth, decorrelating receivers, reading from both fiber ends |
| `s15_budget.py`, `s39_design_example.py`, `s40_temperature.py`, `s56_budget_e2e.py`, `s57_example_e2e.py` | error budgets and worked design examples, end-to-end runs |
| `s38_references.py`, `s54_refcal.py`, `s55_chirp_memory.py`, `s62_two_class.py`, `s70_diag.py` | reference gratings, tuning-axis calibration, pattern-dependent chirp, diagnostics of the rate-equation laser |
| `s0`-`s11`, `s17`, `s24`, `s28`, `s34`-`s37`, `s46`-`s47`, `s52` | earlier explanatory figures and the first round of the study (self-calibration by co-coded references, kept for the record) |
| `common.py`, `figstyle.py`, `s70_field.py` | shared physics of the power model, the figure style, the coherent model |
| `test_selfcal.py` | checks of the code bounds, estimators, ghost delay algebra and acquisition formulas (`python test_selfcal.py`) |

Cached results of the longer scripts are in `out/` (npz and txt).

## Hardware parameters

The default parameters follow the measurement setup of the project: FBGS DTG gratings of 250-pm width at 10 % reflectivity, a BW10 HCG-VCSEL at 1550 nm with its measured tuning curve, m-sequences of 127 to 511 chips at 25 to 100 Mchip/s. Units: optical frequency in GHz, 1 GHz is about 8 pm at 1550 nm.
