# HYPERION Breeder — Activation & Shutdown-Dose (OpenMC R2S)
### Kronos Fusion Energy · design-spec low-activation alloy · package for ACTINV cross-check

A complete, reproducible activation / waste-class / shutdown-dose analysis of the **HYPERION**
compact spherical-tokamak breeder plasma-facing and blanket structure, at the **qualified
low-activation design alloy**, for independent cross-check by **ACTINV** (AvilaLabs). Clean physics
only: material composition + neutron flux + activation inventory + waste class + shutdown dose.

---

## 1. Design point (Tier-2 freeze)

| Quantity | Value |
|---|---|
| Machine | HYPERION spherical-tokamak breeder (D–T), δ = −0.30 |
| Fusion power P_fus | 85.04 MW |
| Neutron source rate S_n | 3.017 × 10¹⁹ n/s (14.06 MeV) |
| Neutron wall loading | 1.885 MW/m² |
| First-wall scalar flux φ_FW | 2.72 × 10¹⁴ n/cm²/s (this model) |
| Exposure | 2 full-power-years |
| Cooling times | shutdown · 1 d · 1 wk · 1 mo · 1 yr · 10 · 100 · 1000 yr |

## 2. Materials — the qualified low-activation alloy (per element, per component)

A single **Cr–Ti–V–W** refractory high-entropy backbone, deliberately **tantalum-, rhenium-,
copper- and molybdenum-free** (Kronos design study, DOI 10.5281/zenodo.22645891). Compositions are
atomic fractions of the SQS cell:

| Component | Alloy | Cr | Ti | V | W | density |
|---|---|---|---|---|---|---|
| First wall | Cr₁Ti₇V₄W₄ (Ti-rich, ductile) | 6.25 | 43.75 | 25.0 | 25.0 at% | 7.70 g/cm³ |
| Blanket structure | Cr₁Ti₇V₄W₄ | 6.25 | 43.75 | 25.0 | 25.0 at% | 7.70 g/cm³ |
| Divertor | Cr₂Ti₄V₄W₆ (W-rich) | 12.5 | 25.0 | 25.0 | 37.5 at% | 9.50 g/cm³ |

**Impurity control (procurement spec, weight-ppm ceilings):** Nb ≤ 5 · Mo ≤ 50 · Co ≤ 10 ·
Ni ≤ 100 · Cu ≤ 10 · Ag ≤ 5 · Zr ≤ 10 · **N ≤ 50** (nitrogen). These are achievable
reduced-activation feedstock limits; the waste class below is stated **conditional on meeting this
impurity spec** — a normal QA / procurement requirement.

## 3. Method (R2S)

OpenMC fixed-source **transport** of a 1-D multi-region slab (first wall 0–2 cm, blanket 2–30 cm,
divertor), 14.06 MeV surface source, 100 batches × 2×10⁵ particles → per-region energy-binned flux
(300-group). Each region is **depleted** for 2 full-power-years at its own local flux on the **full
ENDF/B-VIII.0 depletion chain**, then decayed through the cooling set. Per-component classification
uses **U.S. 10 CFR 61 §61.55** (Class A / C / GTCC) on the principal fusion nuclides
(C-14, Ni-59, Ni-63, Nb-94, Tc-99, Sr-90, Cs-137, I-129); IAEA-style residuals are tracked in the
long-lived breakout. **Honest fidelity:** this is an ENDF/B-VIII.0 OpenMC R2S, not FISPACT-II;
nitrogen is carried as a controlled impurity and its C-14 channel is resolved as a sensitivity.

## 4. Result

**Every component classifies at or below Class C at 100 yr — no greater-than-Class-C (GTCC) waste —
with an order-of-magnitude margin, and this holds robustly across nitrogen 10–100 wppm.**

| Component | 100-yr Class-C index | Class |
|---|---|---|
| First wall | 0.07 | **Class C** |
| Blanket (front/mid/back) | 0.05–0.10 | **Class C** |
| Divertor | 0.09 | **Class C** |

The result is a **computed consequence of the copper-, molybdenum-, tantalum- and rhenium-free
composition**: the Class-C index is governed by Nb-94 (trace Nb impurity) with C-14 (trace nitrogen)
second; both sit far below the regulatory boundary. The waste class is **conditional on the stated
impurity spec**; tighter Nb (≤ 1 wppm) and N (≤ 10 wppm) increase the margin further.

**Shutdown dose (contact / self-dose estimate).** From the same inventories we derive the decay-photon
source per region per cooling time and a **contact surface dose rate** (`results/shutdown_dose.csv`).
The dose descends by ~9 orders of magnitude across the cooling set:

| Cooling time | Contact dose rate (µSv/h) | Note |
|---|---|---|
| Shutdown | ~6 × 10⁹ – 1.4 × 10¹⁰ | remote handling only (peak: blanket front) |
| 1 yr | ~5 – 12 × 10⁷ | still remote-only |
| 10 yr | ~5 – 12 × 10⁴ | |
| 100 yr | ~17 – 54 | approaching hands-on range |
| 1000 yr | ~7 – 15 | low-level residual (Nb-94 photons) |

Contact dose is highest at the **blanket front** (highest neutron flux) and lower at the divertor
(tungsten self-shielding). **Honest fidelity:** these are **first-order CONTACT / self-dose
estimates, not a 3-D photon-transport dose map.** Each component is treated as a semi-infinite
homogeneous slab with a uniform volumetric decay-photon source S_v(E); the surface scalar flux is
φ_s(E) = S_v(E) / (2 µ(E)), with µ(E) the material's own total photon attenuation coefficient
(coherent + incoherent + photoelectric + pair, ENDF/B-VIII.0 photon data) supplying self-shielding;
the ambient dose rate is H = Σ φ_s(E)·h(E) with h(E) the ICRP-74 photon fluence-to-ambient-dose (AP)
coefficient. This is the standard R2S contact/self-dose approximation and **assumes an infinite
uniform source with no geometric attenuation, ducts, streaming or shine** — a bounding surface
estimate, not a room dose map. Decay-heat descends on the standard curve (see `results/`).

## 5. Files

```
README.md                     this file
actinv_problem.json           ACTINV problem definition (materials + flux + schedule + classification)
openmc_deck/                  materials.xml · geometry.xml · settings.xml · tallies.xml
results/                      statepoint.100.h5 (energy-binned flux) · egrid.npy · flux_actinv.json
                              · dep_*_N50.h5 inventories · classC_by_N.csv · designspec_results.json
                              · shutdown_dose.csv (decay-photon source + contact dose per region/cooling)
scripts/                      r2s_design_spec.py · compute_shutdown_dose.py (reproduce)
```

## 6. Reproduce / ingest

```bash
# conda create -n openmc -c conda-forge openmc   (0.16.x)
export OPENMC_CROSS_SECTIONS=/path/to/endfb-viii.0-hdf5/cross_sections.xml
export OPENMC_CHAIN=/path/to/chain_endfb80.xml
python scripts/r2s_design_spec.py         # transport + depletion + waste class
python scripts/compute_shutdown_dose.py   # decay-photon source + contact dose (reads dep_*_N50.h5)
actinv import-flux openmc results/statepoint.100.h5
actinv run actinv_problem.json
```

The transport + depletion + waste-class + shutdown-dose chain runs end-to-end on a single
workstation CPU (no cluster, no GPU); the shutdown-dose back-end is post-processing of the existing
inventories and completes in seconds.

Nuclear data: transport + depletion on **ENDF/B-VIII.0** (full HDF5 library + full depletion chain).

*© 2026 Kronos Fusion Energy · released for DD cross-check. Physics only — no BOM, CAD, supplier or IP.*
