# coherent - the field-level model behind Figs. 3, 10, 11 and 12

Scripts that run the coherent model `../s70_field.py` (coupled-mode gratings, rate-equation VCSEL with chirp and a Lorentzian line, every order of multiple reflections, shot and thermal noise, Bessel receiver) and decode its records with the same chain as the power model of the paper. The raw records are large (about 400 MB for everything used in the paper) and are not kept in the repository. Instead `cache/` holds the decoded quantities that the figure scripts need, so every figure of this folder can be redrawn without recomputing anything. Running the campaign scripts fills `records/` and the figure scripts then use the records in place of the cache.

## Figures

| Figure | Script | Input |
|---|---|---|
| Fig. 3 (critical reflectivity, grating order, all 40320 orders) | `fig_validation.py` | `cache/order_perm.npz` from `permutations.py`, direct coupled-mode spectra computed on the fly |
| Fig. 10 (worst-case bound against both models) | `fig_bound.py` | `../out/s55_results.npz` (power model, `../s55_error_bound.py`) and `cache/bound_field.npz` decoded from the `capR`/`capRb` records of `campaign.py` |
| Fig. 11 (shadowing shift, beat-limited SNR) | `fig_validation.py` | `cache/ruleA_shift.npz` from `runs.py ruleA`, `cache/beat_floor_rows.npy` from `beat_floor.py` |
| Fig. 12 (50-sensor array of [JLT] and the designed array) | `fig_jlt.py` | `cache/jlt/jltc_*_dec.npz` from `jlt_chain.py` and `jlt_decode.py` |

Figures are written to `../figs/`.

## Files

| File | What it does |
|---|---|
| `decode.py` | readout of a record: 12-bit quantization, circular correlation with the bipolar replica, read at the nearest sample of the nominal delay, Gaussian fit, optional sequential deshadowing, reference correction with a drifted tuning table |
| `runs.py` | thin wrapper around `s70_field.run_array` with a record cache in `records/`, plus `ruleA` (two gratings, shift of the second against the detuning of the first) |
| `campaign.py` | records for the bound check: `capR` (300-MHz line), `capRb` (1-GHz line), `capRd` (direct paths only) |
| `permutations.py` | all orders of the eight gratings of Fig. 3(b): largest Omega_k, split returns, shift of the maximum |
| `jlt_chain.py` | the 50-sensor array of [JLT] and its variants (D0, J5, J5m, J5m50, D4, D4r, ...) as one record per source realization |
| `jlt_decode.py` | decoding of those records with the estimator of the paper (Gaussian fit within +-450 pm of the nominal wavelength, reference correction of offset and slope) |
| `beat_floor.py` | beat-noise floor of the correlogram for returns sharing one wavelength step, against the closed form (eq:beat) |
| `fig_validation.py`, `fig_bound.py`, `fig_jlt.py` | the figures listed above |

## Running

```bash
cd python_selfcal_cdm/coherent
python fig_validation.py      # Figs. 3 and 11 from the cache, a few seconds
python fig_bound.py           # Fig. 10 from the cache
python fig_jlt.py             # Fig. 12 from the cache

python runs.py ruleA          # recompute the shadowing shift (12 records of two gratings, minutes)
python permutations.py        # recompute the 40320 orders (about a minute)
python campaign.py capR       # 32 records, up to 32 gratings each, about an hour
python campaign.py capRb
python beat_floor.py          # 9 one-step records with 17 and 50 gratings
python jlt_chain.py J5 3      # one realization of the 50-sensor array, then
python jlt_decode.py          # decode every record in records/jlt/
```

Variants of `jlt_chain.py` used in the paper: `D0` (array as in [JLT]: co-tuned, uniform 2.5-m spacing, 31.25 Mchip/s, random code), `J5` (D0 with a +-5-cm spacing tolerance, seeds 1-8), `J5m` (m-sequence), `J5m50` (50 Mchip/s), `D4r` (irregular spacing, three bands, two references, seed 1), `D4rL3` (D4r with a 3-GHz line), `D4jr`, `D4yr`, `D4rR3`, `D4rR5` (position tolerance, detuned sensors, 3 % and 5 % reflectivity). Suffix `L<n>` sets the source line in GHz, `R<n>` the reflectivity in percent, `d` keeps the direct paths only.

[JLT]: K. Markowski et al., "Analysis of the performance of WDM-CDM Bragg grating interrogation system with high-contrast grating VCSEL," J. Lightwave Technol. 41(9), 2892-2903, 2023, doi:10.1109/JLT.2023.3237602.
