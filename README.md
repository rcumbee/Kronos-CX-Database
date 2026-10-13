# Kronos CX database v3.4 (work in progress)

Charge-exchange (CX) cross sections and radiative data for X-ray spectra from collisions of highly charged ions
with atoms and molecules (H, He, H2, H2O, CO, CO2, N2, O, OH, ...), and the Kronos code that turns them into spectra.

- `CXDatabase/Projectile_Ions/<El>/Charge/<q>/Targets/<Target>/`: cross sections (`*_sec_<method>.cs`), level list (`*.nist.in`)
  and radiative data (`radrathe.dat`) for each collision system. File header line 5 gives the reference for the data.
- `corrections.dat`: every change made to the database since the 2017 release (newest first).
- `changelog.txt`: versions of the database and code.
- `license.txt`: GNU General Public License v3.

## Kronos_GUI v1 (`kronos_gui.py`)

A single Python file: pick a collision system and collision energy, run the radiative cascade, get the spectrum (photon energy,
emission cross section) and the line list. Needs Python 3.8 or newer with `numpy`, `pandas` and `matplotlib`; the window also
needs `PySide6` (`pip install PySide6`).

The simplest way to run Kronos is by running the GUI (needs PySide6). In the terminal type:
```
python3 kronos_gui.py --gui
```
and choose your collision system.

Python / Jupyter (see `Kronos_example.ipynb`):

```python
from kronos_gui import Kronos, list_systems
list_systems()                                   # ions, charges, targets, methods
k = Kronos('O', 8, 'H', method='mclz_nres', energy=1000)   # energy in eV/u; nearest tabulated energy is used
k.spectrum                                       # Energy (eV), Intensity = emission cross section of each line (cm^2)
k.line_list(top=10)                              # strongest lines, upper and lower states
k.line_ratios()                                  # n -> 1 series relative to alpha
k.plot(fwhm=10)                                  # stick + Gaussian-convolved spectrum
k.save('results')                                # spectrum and line list as CSV
```

Terminal:

```
python3 kronos_gui.py O 8 H --energy 1000 --method mclz_nres --plot --fwhm 10
python3 kronos_gui.py --list Fe 25
python3 kronos_gui.py --gui
```

The emission cross section of a line u -> l is σ_u A_ul / Σ_k A_uk: the population of the upper level u (its capture cross
section plus all cascade feeding into it) times the line's branching ratio. Times the target density and the ion flux, it gives
the photon emission rate in that line.


### Default method

If no method is given (`-M` / `--method`, or `method=` in Python), Kronos uses the first one in this order whose tabulated
collision energies cover the requested energy (within 5%): RCMD → QMOCC → MOCC → AOCC → TC-BGM → CTMC → CDW → MCLZ
(MCLZ includes the n-resolved `mclz_nres` files). If none covers it, MCLZ is used. The nearest tabulated energy of the chosen
file is used, and the energy actually used is printed.

### l distributions for n-resolved cross sections

For bare projectiles (H-like products) the database gives cross sections resolved by n only (`*_mclz_nres.cs`). Kronos spreads
each n over l = 0 ... n-1 with a weight W(n, l), chosen with `--ldist` / `ldist=`:

| `ldist` | W(n, l) | Source | ACX counterpart |
|---|---|---|---|
| `low-energy` (default) | $(2l+1)\,\frac{[(n-1)!]^2}{(n+l)!\,(n-l-1)!}$ | Krasnopolsky et al. 2004 | none |
| `statistical` | $\frac{2l+1}{n^2}$ | Krasnopolsky et al. 2004; Janev & Winter 1985 | ACX 2, statistical (and 6) |
| `sl1` | the low-energy distribution shifted down by one l: W(n, l) = W_low-energy(n, l+1), with p also taking W_low-energy(n, 0) and W(n, n-1) = 0 | Mullen et al. 2016 | none |
| `even` | $\frac{1}{n}$ | Smith et al. 2014 | ACX 1, even (and 5) |
| `landau-zener` | $l(l+1)(2l+1)\,\frac{(n-1)!\,(n-2)!}{(n+l)!\,(n-l-1)!}$ (W = 0 for n = 1) | Smith et al. 2014; Janev & Winter 1985 (Eq. 3.50a) | ACX 3, Landau-Zener (and 7) |
| `separable` | $\frac{2l+1}{q}\,\exp\left[-\frac{l(l+1)}{q}\right]$ | Smith et al. 2014; Janev & Winter 1985 (Eq. 3.59) | ACX 4, separable (and 8) |

q is the projectile charge (for a bare projectile, q = Z). `low-energy`, `sl1` and `landau-zener` put most of the capture into low l
(for n = 5-13 they peak at l = 1-2, l = 1 and l = 2-4); `statistical` grows with l (largest at l = n-1) and is meant for high
collision velocities; `even` is flat. The `separable` weights are not normalized exactly (they sum to 0.99-1.01 for
n = 5-13) and are used as they are. ACX = the AtomDB charge-exchange model (Smith et al. 2014, ApJ 787, 77), whose
l-distribution options (`acxmodel`, ACX2 documentation) are 1 (even), 2 (statistical), 3 (Landau-Zener), 4 (separable);
ACX models 5-8 use the same four l distributions with the capture split over two n shells (the n distribution does not
apply here: the database gives n).

### Known limitations

- Fe8+ and Fe9+ have cross sections but no `radrathe.dat` yet: their cross sections can be plotted, but no cascade can be run.
- He-like products whose `radrathe.dat` is resolved by term only (LS terms, e.g. Fe25+): the 1s2p 3P term decays with a single
  set of rates, so the forbidden line z (1s2s 3S1 -> 1s2) may out too weak and x and y will appear as one line. This was tested
  for Fe25+; it is likely to affect higher-Z ions more than lower-Z ions.

## Earlier versions

The original Kronos code (v1 to v3.3: Python 2 and gfortran, with its PDF manual) is in this repository's history, one tag per
version (`v1`, `v2.5`, ..., `v3.3`).

## References

Please cite the papers given in the header of each file you use, and Mullen et al. (2016, ApJS 224, 31),
Cumbee et al. (2018, ApJ 852, 7), Lyons et al. (2017, ApJS 232, 27) and Mullen et al. (2017, ApJ 844, 7) for Kronos.
