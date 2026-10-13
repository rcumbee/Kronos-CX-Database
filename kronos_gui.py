#!/usr/bin/env python3
"""
Kronos_GUI v1 -- charge-exchange X-ray spectra from the Kronos CX database (v3.4).

One file: reads the database next to it (CXDatabase/Projectile_Ions/<El>/Charge/<q>/Targets/<Target>/), runs the radiative
cascade for one collision system at one collision energy, and returns the spectrum (photon energy, emission cross section) and
the line list. A small Qt window (PySide6) does the same with menus and a plot.

Use
  Python / Jupyter:
      from kronos_gui import Kronos, list_systems
      k = Kronos('O', 8, 'H', energy=1000)          # method: preferred one covering 1000 eV/u (see preferred_method)
      spec = k.spectrum                               # DataFrame: Energy (eV), Intensity (cm^2)
      lines = k.line_list()                           # strongest lines with upper/lower states
      k.plot(fwhm=10)                                 # stick spectrum + Gaussian-convolved spectrum
  Terminal:
      python kronos_gui.py O 8 H --energy 1000 [--method mclz_nres --ldist low-energy --fwhm 10 --plot --out results/]
      python kronos_gui.py --list O 8               # targets and methods available
      python kronos_gui.py --gui                    # the window (needs PySide6: pip install PySide6)
Needs numpy, pandas, matplotlib (and PySide6 for the window).

Physics, as in Kronos: the cross section of each capture state (nl term, or n spread over l for n-resolved files) is cascaded
through the radiative transitions in the folder's radrathe.dat (Einstein A coefficients); each state's population is shared
over its decays in proportion to A. Emission cross sections are in cm^2 (file units are taken from header line 7).
"""
import argparse
import heapq
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

__version__ = '1.0'

# ----------------------------------------------------------------------------------------------------------- constants

L_LETTERS = 'spdfghiklmnoqrtuvwxyz'                       # l = 0, 1, 2, ... (no j)
L_DICT = {c: i for i, c in enumerate(L_LETTERS)}
S_CODE = {'1': 0, '3': 1, '2': 2, '4': 4, '5': 5}          # 2S+1 in a label -> spin code used in radrathe.dat
CODE_TO_MULT = {v: k for k, v in S_CODE.items()}
ELEMENTS = {s: i + 1 for i, s in enumerate(
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn'.split())}
LDISTS = ['low-energy', 'statistical', 'sl1', 'even', 'landau-zener', 'separable']
# old Kronos names (the AtomDB ACX acxmodel numbers 1, 3, 4), still accepted
LDIST_ALIASES = {'ACX1': 'even', 'ACX3': 'landau-zener', 'ACX4': 'separable'}
# preference order when no method is given (the original Kronos v3 order + TC-BGM, CDW); 'mclz' covers mclz and mclz_nres;
# only methods whose energy range covers the requested energy are taken (preferred_method)
METHOD_ORDER = ['rcmd', 'qmocc', 'mocc', 'aocc', 'tcbgm', 'ctmc', 'cdw', 'mclz']
# H-like 2s -> 1s two-photon (2E1) rate: A = 8.2283 Z^6 s^-1 (non-relativistic; e.g. Parpia & Johnson 1982, Goldman & Drake 1981)
A_2E1_H = 8.2283
ENERGY_TOL = 0.05        # a method 'covers' a collision energy within 5% of its tabulated range (default method choice)
SERIES = {1: 'alpha', 2: 'beta', 3: 'gamma', 4: 'delta', 5: 'epsilon', 6: 'zeta', 7: 'eta'}
GREEK = {1: 'α', 2: 'β', 3: 'γ', 4: 'δ', 5: 'ε', 6: 'ζ', 7: 'η'}


def default_db():
    """The database next to this file (CXDatabase/), or $KRONOSDB."""
    here = Path(__file__).resolve().parent / 'CXDatabase'
    if here.is_dir():
        return here
    env = os.environ.get('KRONOSDB')
    if env:
        return Path(env)
    raise FileNotFoundError('No CXDatabase folder next to kronos_gui.py and KRONOSDB is not set; pass db=...')


# ----------------------------------------------------------------------------------------------------------- database

def _charge_dir(db, ion):
    for name in ('Charge', 'charge'):
        d = Path(db) / 'Projectile_Ions' / ion / name
        if d.is_dir():
            return d
    return Path(db) / 'Projectile_Ions' / ion / 'Charge'


def system_dir(db, ion, charge, target):
    return _charge_dir(db, ion) / str(int(charge)) / 'Targets' / target


def _method_of(fname, ion, charge):
    """'o8+h_sec_mclz_nres.cs' -> 'mclz_nres'; None for files Kronos does not use."""
    f = fname.lower()
    if not f.endswith('.cs') or not f.startswith(f'{ion.lower()}{int(charge)}+'):
        return None
    m = re.match(r'.*?_sec_(.+)\.cs$', f)
    return m.group(1) if m else None


def methods(db, ion, charge, target):
    """{method: file} for one system. n-resolved files are used for bare projectiles only (H-like products)."""
    d = system_dir(db, ion, charge, target)
    bare = int(charge) == ELEMENTS[ion]
    out = {}
    for f in sorted(os.listdir(d)) if d.is_dir() else []:
        m = _method_of(f, ion, charge)
        if m and (bare or 'nres' not in m) and not m.endswith('_j'):
            out[m] = d / f
    return out


def energy_range(path):
    """(min, max) tabulated collision energy (eV/u) of a .cs file."""
    e = read_cs(path).iloc[:, 0]
    return float(e.min()), float(e.max())


def preferred_method(db, ion, charge, target, energy, tol=ENERGY_TOL):
    """Default method for one system at one collision energy: the first method in METHOD_ORDER whose tabulated energy range
    covers the energy within a relative tolerance tol (e.g. 997 eV/u still uses a file that starts at 1000 eV/u; a file with one
    energy E0 serves E0 +- 5%); if none covers it, MCLZ; if there is no MCLZ, the first available in METHOD_ORDER."""
    avail = methods(db, ion, charge, target)
    if not avail:
        return None
    ranked = [m for want in METHOD_ORDER for m in sorted(avail) if m.split('_')[0] == want]
    ranked += [m for m in sorted(avail) if m not in ranked]
    e = float(energy)
    for m in ranked:
        lo, hi = energy_range(avail[m])
        if lo / (1 + tol) <= e <= hi * (1 + tol):
            return m
    return next((m for m in ranked if m.split('_')[0] == 'mclz'), ranked[0])


def list_systems(db=None):
    """DataFrame of every system in the database: ion, charge, target, methods, has radrathe.dat."""
    db = Path(db) if db else default_db()
    rows = []
    for ion in sorted((p.name for p in (db / 'Projectile_Ions').iterdir() if p.is_dir()), key=lambda s: ELEMENTS.get(s, 99)):
        cd = _charge_dir(db, ion)
        for q in sorted((p.name for p in cd.iterdir() if p.is_dir() and p.name.isdigit()), key=int):
            for t in sorted(p.name for p in (cd / q / 'Targets').iterdir() if p.is_dir()):
                ms = methods(db, ion, q, t)
                rows.append(dict(ion=ion, charge=int(q), target=t, methods=', '.join(ms),
                                 radrathe=(system_dir(db, ion, q, t) / 'radrathe.dat').exists()))
    return pd.DataFrame(rows)


def read_cs(path):
    """A Kronos .cs file -> DataFrame (first column '(eV/u)'), with attrs citation (header line 5) and units (line 7)."""
    lines = Path(path).read_text(encoding='latin1').splitlines()
    head = lines[:8]
    if len(head) < 8 or not all(l.startswith('#') for l in head):
        raise ValueError(f'{path}: expected an 8-line "#" header')
    cols = head[7].replace('#', ' ').split()
    rows = [[float(x) for x in l.split()] for l in lines[8:] if l.strip() and not l.startswith('#')]
    df = pd.DataFrame(rows, columns=cols[:len(rows[0])] if rows else cols)
    df.attrs['citation'] = head[4].lstrip('#').strip()
    m = re.search(r'\((.*?)\)', head[6])
    df.attrs['units'] = m.group(1) if m else 'cm^2'
    return df


def unit_factor(units):
    """Header units -> cm^2 (e.g. '10^-16 cm^2' -> 1e-16)."""
    m = re.search(r'10\^\s*(-?\d+)', units)
    return 10.0 ** int(m.group(1)) if m else 1.0


# ----------------------------------------------------------------------------------------------------------- l distributions

def l_weight(ldist, n, l, Z):
    """Fraction of an n-resolved cross section put into l (Krasnopolsky et al. 2004; Smith et al. 2014; Mullen et al. 2016).
    even, statistical, landau-zener, separable = AtomDB ACX l distributions 1, 2, 3, 4 (Smith et al. 2014)."""
    f = math.factorial
    ldist = LDIST_ALIASES.get(ldist, ldist)
    if ldist == 'statistical':
        return (2 * l + 1) / n ** 2
    if ldist == 'low-energy':
        return (2 * l + 1) * f(n - 1) ** 2 / (f(n + l) * f(n - 1 - l))
    if ldist == 'even':                                   # ACX 1
        return 1.0 / n
    if ldist == 'landau-zener':                            # ACX 3 (Janev & Winter 1985 Eq. 3.50a)
        return 0.0 if n == 1 else l * (l + 1) * (2 * l + 1) * f(n - 1) * f(n - 2) / (f(n + l) * f(n - l - 1))
    if ldist == 'separable':                               # ACX 4 (Janev & Winter 1985 Eq. 3.59)
        return (2 * l + 1) / Z * math.exp(-l * (l + 1) / Z)
    if ldist == 'sl1':
        if l == n - 1:
            return 0.0
        if l == 1:
            return l_weight('low-energy', n, 2, Z) + l_weight('low-energy', n, 0, Z)
        return l_weight('low-energy', n, l + 1, Z)
    raise ValueError(f'unknown l distribution {ldist!r}; choose from {LDISTS}')


# ----------------------------------------------------------------------------------------------------------- radiative data

def read_radrathe(path):
    """radrathe.dat -> DataFrame with integer state columns and E (wavelength, Angstrom), A (s^-1).

    12 columns: n l L S counter (upper) n l L S counter (lower) E A; 8 columns: n l S, n l S, E, A; 6 columns: n l, n l, E, A.
    Files that give 2S+1 instead of the spin code (a 3 anywhere) are converted."""
    df = pd.read_csv(path, sep=r'\s+', comment='#', header=None)
    names = {12: ['ni', 'li', 'Li', 'Si', 'ci', 'nf', 'lf', 'Lf', 'Sf', 'cf', 'E', 'A'],
             8: ['ni', 'li', 'Si', 'nf', 'lf', 'Sf', 'E', 'A'], 6: ['ni', 'li', 'nf', 'lf', 'E', 'A']}
    if len(df.columns) not in names:
        raise ValueError(f'{path}: cannot read a radrathe file with {len(df.columns)} columns')
    df.columns = names[len(df.columns)]
    if 'Si' in df and ((df.Si == 3).any() or (df.Sf == 3).any()):
        conv = {int(k): v for k, v in S_CODE.items()}
        conv[0] = 0
        for c in ('Si', 'Sf'):
            df[c] = df[c].map(lambda s: conv.get(int(s), int(s)))
    for c in df.columns[:-2]:
        df[c] = df[c].astype(int)
    return df


def _state_cols(df, side):
    order = ['n', 'l', 'L', 'S', 'c']
    return [f'{x}{side}' for x in order if f'{x}{side}' in df.columns]


def _cascade(rr, pop, extra=None):
    """Cascade the populations pop {state tuple: cross section} through radrathe rows rr.

    extra {state: A (s^-1)}: decays not in radrathe.dat that give no line (two-photon 2s -> 1s of H-like ions); they compete
    with the listed decays, and the flux they carry is returned (cross section per state).

    States are processed upper-first (topological order of the decay graph, ties by
    n descending; plain n order if the transitions form a loop), and each state's population is shared over its decays in
    proportion to A. Returns the transitions with an emitted cross section."""
    ci, cf = _state_cols(rr, 'i'), _state_cols(rr, 'f')
    up = [tuple(r) for r in rr[ci].to_numpy()]
    lo = [tuple(r) for r in rr[cf].to_numpy()]
    A = rr['A'].to_numpy(float); E = rr['E'].to_numpy(float)
    nup = rr['ni'].to_numpy()
    # topological order of the upper states
    n_of = {}
    for s, n in zip(up, nup):
        n_of.setdefault(s, n)
    for s, n in zip(lo, rr['nf'].to_numpy()):
        n_of.setdefault(s, n)
    feeds = {s: set() for s in n_of}
    indeg = dict.fromkeys(n_of, 0)
    for a, b in zip(up, lo):
        if a != b and b not in feeds[a]:
            feeds[a].add(b); indeg[b] += 1
    heap = [(-n_of[s], s) for s, d in indeg.items() if d == 0]
    heapq.heapify(heap)
    rank = {}
    while heap:
        _, s = heapq.heappop(heap)
        rank[s] = len(rank)
        for b in feeds[s]:
            indeg[b] -= 1
            if indeg[b] == 0:
                heapq.heappush(heap, (-n_of[b], b))
    if len(rank) < len(n_of):
        order = np.argsort(-nup, kind='stable')            # loop in the transitions: upper-n order
    else:
        order = sorted(range(len(up)), key=lambda i: (rank[up[i]], i))
    asum = {}
    for s, a in zip(up, A):
        asum[s] = asum.get(s, 0.0) + a
    extra = {s: a for s, a in (extra or {}).items() if s in asum}
    for s, a in extra.items():
        asum[s] += a
    cs = dict(pop)
    out = []
    done = set()
    for i in order:
        s, f = up[i], lo[i]
        if asum[s] <= 0:
            continue
        flux = cs.get(s, 0.0) * A[i] / asum[s]
        cs[f] = cs.get(f, 0.0) + flux
        done.add(s)
        if flux > 0 and E[i] > 0:
            out.append((s, f, 1.0e8 / E[i] / 8065.5, flux))   # photon energy in eV
    lineless = {s: cs.get(s, 0.0) * a / asum[s] for s, a in extra.items() if s in done}
    return out, ci, lineless


def _term_label(state, cols):
    """(n, l, [L,] [S, counter]) -> '2p(1P)' style label."""
    d = dict(zip([c[:-1] for c in cols], state))
    n, l = d['n'], d['l']
    lab = f'{n}{L_LETTERS[l] if l < len(L_LETTERS) else "?"}'
    if 'S' in d:
        L = d.get('L', l)
        Lc = L_LETTERS[L].upper() if L < len(L_LETTERS) else '?'
        lab += f'({CODE_TO_MULT.get(d["S"], "?")}{Lc})'
    return lab


# ----------------------------------------------------------------------------------------------------------- model

class Kronos:
    """One collision system at one collision energy.

    Kronos(ion, charge, target, method=None, energy=1000.0, ldist='low-energy', db=None, spin_dist=(0.25, 0.75))
      ion, charge, target: projectile element, projectile charge, target name as in the database (e.g. 'O', 8, 'H2O')
      method: file method name (see methods()); None = preferred_method(): the first in METHOD_ORDER whose energy range
              covers the energy (within 5%), else MCLZ
      energy: collision energy in eV/u; the nearest tabulated energy is used (energy_used)
      ldist: l distribution for n-resolved files (bare projectiles): low-energy (default), statistical, sl1, even, landau-zener,
             separable (the last three are ACX 1, 3, 4; the names ACX1, ACX3, ACX4 are also accepted)
      spin_dist: (singlet, triplet) split for spin-summed nl columns of He-like products
    Attributes after construction: cross_sections (the whole file), citation, units, energy_used, spectrum, transitions.
    """

    def __init__(self, ion, charge, target, method=None, energy=1000.0, ldist='low-energy', db=None, spin_dist=(0.25, 0.75),
                 verbose=True):
        self.db = Path(db) if db else default_db()
        self.ion, self.charge, self.target = ion, int(charge), target
        if ion not in ELEMENTS:
            raise ValueError(f'unknown element {ion!r}')
        self.Z = ELEMENTS[ion]
        self.folder = system_dir(self.db, ion, charge, target)
        if not self.folder.is_dir():
            raise FileNotFoundError(f'no folder for {ion}{charge}+ + {target}: {self.folder}')
        avail = methods(self.db, ion, charge, target)
        if not avail:
            raise FileNotFoundError(f'no cross-section files for {ion}{charge}+ + {target}')
        if method is None:
            method = preferred_method(self.db, ion, charge, target, energy)
        if method not in avail:
            raise ValueError(f'method {method!r} not available for {ion}{charge}+ + {target}; choose from {list(avail)}')
        self.method = method
        self.nres = 'nres' in method
        ldist = LDIST_ALIASES.get(ldist, ldist)
        if self.nres and ldist not in LDISTS:
            raise ValueError(f'unknown l distribution {ldist!r}; choose from {LDISTS}')
        self.ldist = ldist if self.nres else 'l-resolved'
        self.spin_dist = tuple(spin_dist)
        self.verbose = verbose
        self.cs_file = avail[method]
        self.cross_sections = read_cs(self.cs_file)
        self.citation = self.cross_sections.attrs['citation']
        self.units = self.cross_sections.attrs['units']
        Ecol = self.cross_sections.columns[0]
        idx = int((self.cross_sections[Ecol] - float(energy)).abs().idxmin())
        self.energy_requested = float(energy)
        self.energy_used = float(self.cross_sections.loc[idx, Ecol])
        self._row = self.cross_sections.loc[idx]
        self.run()

    # --- capture populations ---
    def _populations(self, rr):
        """{state tuple: cross section (cm^2)} in the state format of the radrathe file."""
        fac = unit_factor(self.units)
        cols = _state_cols(rr, 'i')
        have = set(tuple(r) for r in rr[cols].to_numpy()) | set(tuple(r) for r in rr[_state_cols(rr, 'f')].to_numpy())
        spins = {}
        for s in have:
            spins.setdefault((s[0], s[1]), set()).add(s[cols.index('Si')] if 'Si' in cols else None)
        he_like = self.charge == self.Z - 1
        pop, dropped = {}, []
        labels = [c for c in self.cross_sections.columns[1:] if c.lower() != 'total']

        def key(n, l, L, S, c=1):
            d = dict(n=n, l=l, L=L, S=S, c=c)
            return tuple(d[x[:-1]] for x in cols)

        for lab in labels:
            v = float(self._row[lab]) * fac
            if v <= 0:
                continue
            if self.nres:
                m = re.match(r'^n=?(\d+)$', lab)
                if not m:
                    dropped.append(lab); continue
                n = int(m.group(1))
                for l in range(n):
                    k = key(n, l, l, S_CODE['2'])
                    pop[k] = pop.get(k, 0.0) + v * l_weight(self.ldist, n, l, self.Z)
                continue
            m = re.match(r'^(\d+)([a-z])\((\d)([A-Za-z])\)(\d*)$', lab)
            if m:
                n, l, S, L = int(m.group(1)), L_DICT[m.group(2)], S_CODE[m.group(3)], L_DICT[m.group(4).lower()]
                k = key(n, l, L, S, int(m.group(5) or 1))
                if 'ci' in rr.columns and k not in have:
                    k = key(n, l, L, S, 1)                 # term file without counters
                pop[k] = pop.get(k, 0.0) + v
                continue
            m = re.match(r'^(\d+)([a-z])$', lab)          # spin-summed nl
            if m:
                n, l = int(m.group(1)), L_DICT[m.group(2)]
                sp = spins.get((n, l), set())
                if he_like and {0, 1} <= sp:
                    for S, f in ((0, self.spin_dist[0]), (1, self.spin_dist[1])):
                        k = key(n, l, l, S); pop[k] = pop.get(k, 0.0) + v * f
                elif len(sp) == 1:
                    k = key(n, l, l, next(iter(sp))); pop[k] = pop.get(k, 0.0) + v
                else:
                    dropped.append(lab)
                continue
            dropped.append(lab)
        missing = [k for k in pop if k not in have]
        if self.verbose and dropped:
            print(f'Note: {len(dropped)} cross-section columns not used (labels not understood): {dropped[:8]}')
        if self.verbose and missing:
            lost = sum(pop[k] for k in missing) / max(sum(pop.values()), 1e-300)
            print(f'Note: {len(missing)} capture states are not in radrathe.dat ({lost:.1%} of the capture, no emission).')
        return pop

    def run(self):
        """Run the cascade (done on construction); sets spectrum and transitions."""
        rr_path = self.folder / 'radrathe.dat'
        if not rr_path.exists():
            raise FileNotFoundError(f'{self.ion}{self.charge}+ + {self.target}: no radrathe.dat in {self.folder} '
                                    '(cross sections only; plot_cross_sections() still works)')
        rr = read_radrathe(rr_path)
        pop = self._populations(rr)
        self.capture_total = sum(pop.values())
        extra = {}
        if self.charge == self.Z:                          # H-like product: 2s -> 1s two-photon decay (2E1), not in radrathe.dat
            a2e1 = A_2E1_H * self.Z ** 6
            for s in set(tuple(r) for r in rr[_state_cols(rr, 'i')].to_numpy()):
                if s[0] == 2 and s[1] == 0:
                    extra[s] = a2e1
        out, cols, lineless = _cascade(rr, pop, extra)
        self.two_photon = sum(lineless.values())           # cross section (cm^2) into the 2s -> 1s two-photon continuum
        t = pd.DataFrame(out, columns=['upper_state', 'lower_state', 'Energy', 'Intensity'])
        t['upper'] = [_term_label(s, cols) for s in t.upper_state]
        t['lower'] = [_term_label(s, [c.replace('i', 'f', 1) if c[-1] == 'i' else c for c in cols]) for s in t.lower_state]
        t['n_upper'] = [s[0] for s in t.upper_state]
        t['n_lower'] = [s[0] for s in t.lower_state]
        prefix = {1: 'Ly', 2: 'K'}.get(self.Z - self.charge + 1, '')
        t['series'] = [f'{prefix}-{SERIES.get(nu - 1, f"n={nu}")}' if nl_ == 1 and nu > 1 and prefix else ''
                       for nu, nl_ in zip(t.n_upper, t.n_lower)]
        self.transitions = t.sort_values('Intensity', ascending=False, ignore_index=True)
        self.spectrum = (t.groupby('Energy', as_index=False).Intensity.sum().sort_values('Energy', ignore_index=True))
        self.spectrum.attrs.update(units='cm^2', citation=self.citation)
        e_tab = self.cross_sections.iloc[:, 0]
        if self.verbose and not (e_tab.min() <= self.energy_requested <= e_tab.max()):
            print(f'Note: {self.energy_requested:g} eV/u is outside the tabulated range of this file '
                  f'({e_tab.min():g}-{e_tab.max():g} eV/u); the nearest energy is used.')
        if self.verbose:
            print(f'{self.ion}{self.charge}+ + {self.target}: {self.method}'
                  + (f' ({self.ldist} l distribution)' if self.nres else '')
                  + f' at {self.energy_used:g} eV/u (requested {self.energy_requested:g}); '
                  f'{len(self.spectrum)} lines. Reference: {self.citation}')
        return self.spectrum

    # --- output ---
    def line_list(self, top=None, min_fraction=0.0):
        """Lines sorted by intensity: upper, lower, series, Energy (eV), Intensity (cm^2), fraction of the total emission."""
        t = self.transitions[['upper', 'lower', 'series', 'Energy', 'Intensity']].copy()
        t['fraction'] = t.Intensity / t.Intensity.sum()
        t = t[t.fraction >= min_fraction]
        return t.head(top).reset_index(drop=True) if top else t.reset_index(drop=True)

    def line_ratios(self, reference='alpha'):
        """Ground-state series lines (n -> 1) relative to the reference line (default the n = 2 -> 1 lines summed)."""
        s = self.transitions[self.transitions.series != '']
        if s.empty:
            raise ValueError('no n -> 1 series lines for this ion')
        g = s.groupby('series', sort=False).agg(n_upper=('n_upper', 'first'), Energy=('Energy', 'mean'), Intensity=('Intensity', 'sum'))
        g = g.sort_values('n_upper')
        ref = g[g.index.str.endswith(reference)]
        g['ratio'] = g.Intensity / (ref.Intensity.iloc[0] if len(ref) else g.Intensity.iloc[0])
        return g.reset_index()

    def convolve(self, fwhm=10.0, erange=None, npts=None):
        """Gaussian-convolved spectrum (Energy, Intensity per eV) with the given FWHM (eV) over erange (eV)."""
        e = self.spectrum.Energy.to_numpy(float); w = self.spectrum.Intensity.to_numpy(float)
        lo, hi = erange if erange else (e.min() - 3 * fwhm, e.max() + 3 * fwhm)
        grid = np.linspace(lo, hi, npts or max(int((hi - lo) / fwhm * 8), 200))
        sig = fwhm / (2 * math.sqrt(2 * math.log(2)))
        m = (e >= lo - 5 * fwhm) & (e <= hi + 5 * fwhm)
        y = (w[m][None, :] * np.exp(-0.5 * ((grid[:, None] - e[m][None, :]) / sig) ** 2)).sum(1) / (sig * math.sqrt(2 * math.pi))
        return pd.DataFrame({'Energy': grid, 'Intensity': y})

    def plot(self, fwhm=None, erange=None, label_lines=8, ax=None, show=True):
        """Stick spectrum (and the Gaussian-convolved spectrum when fwhm is given)."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(9, 4.5))
        s = self.spectrum
        if erange:
            s = s[(s.Energy >= erange[0]) & (s.Energy <= erange[1])]
        ax.vlines(s.Energy, 0, s.Intensity, color='tab:red', lw=0.8, label='lines')
        if fwhm:
            c = self.convolve(fwhm, erange)
            scale = s.Intensity.max() / c.Intensity.max() if c.Intensity.max() > 0 else 1
            ax.plot(c.Energy, c.Intensity * scale, 'k-', lw=1, label=f'Gaussian, FWHM {fwhm:g} eV (scaled)')
        for _, r in self.line_list(top=label_lines).iterrows() if label_lines else []:
            if not erange or erange[0] <= r.Energy <= erange[1]:
                ax.annotate(f'{r.upper}→{r.lower}', (r.Energy, r.Intensity), xytext=(0, 4), textcoords='offset points',
                            rotation=90, fontsize=7, ha='center', va='bottom')
        ax.set_xlabel('photon energy (eV)')
        ax.set_ylabel(r'emission cross section $\sigma_u A_{ul} / \Sigma_k A_{uk}$ (cm$^2$)')
        ax.set_title(f'{self.ion}$^{{{self.charge}+}}$ + {self.target}: {self.method}, {self.energy_used:g} eV/u', fontsize=10)
        ax.legend(fontsize=8)
        if show:
            plt.tight_layout(); plt.show()
        return ax

    def plot_cross_sections(self, top=10, n_resolved=False, ax=None, show=True):
        """Capture cross sections vs collision energy (top states; n_resolved sums nl columns over l)."""
        import matplotlib.pyplot as plt
        df = self.cross_sections.set_index(self.cross_sections.columns[0])
        df = df[[c for c in df.columns if c.lower() != 'total']] * unit_factor(self.units)
        if n_resolved and not self.nres:
            df = df.T.groupby(lambda c: int(re.match(r'(\d+)', c).group(1)) if re.match(r'\d+', c) else 0).sum().T
            df.columns = [f'n={c}' for c in df.columns]
        keep = df.max().nlargest(top).index
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 5))
        df[keep].replace(0, np.nan).plot(ax=ax, logx=True, logy=True)
        ax.set_xlabel('collision energy (eV/u)'); ax.set_ylabel('cross section (cm$^2$)')
        ax.set_title(f'{self.ion}$^{{{self.charge}+}}$ + {self.target}: {self.method} (top {top})', fontsize=10)
        ax.legend(fontsize=7, ncol=2)
        if show:
            plt.tight_layout(); plt.show()
        return ax

    def save(self, outdir='.'):
        """Write <system>_<method>_<E>_spectrum.csv and _lines.csv (with the reference in the first line)."""
        out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
        stem = f'{self.ion.lower()}{self.charge}+{self.target.lower()}_{self.method}' + (f'_{self.ldist}' if self.nres else '') + f'_{self.energy_used:g}eVu'
        paths = []
        for name, df in (('spectrum', self.spectrum), ('lines', self.line_list())):
            p = out / f'{stem}_{name}.csv'
            with open(p, 'w') as fh:
                fh.write(f'# {self.ion}{self.charge}+ + {self.target}, {self.method}, {self.energy_used:g} eV/u; '
                         f'energies in eV, intensities in cm^2; reference: {self.citation}\n')
                df.to_csv(fh, index=False)
            paths.append(p)
        return paths


# ----------------------------------------------------------------------------------------------------------- GUI

def gui(db=None, start=True):
    """Qt window (PySide6): choose a system, plot its cross sections, run the cascade, plot and save the spectrum.
    start=False returns the window without starting the event loop (for scripted tests)."""
    try:
        import matplotlib
        matplotlib.use('QtAgg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
        from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QFormLayout, QComboBox,
                                       QCheckBox, QLineEdit, QPushButton, QLabel, QSplitter, QMessageBox, QFileDialog,
                                       QScrollArea)
    except ImportError as ex:
        sys.exit(f'The window needs PySide6 (pip install PySide6): {ex}')
    db = Path(db) if db else default_db()
    systems = list_systems(db)

    class Plot(QWidget):
        def __init__(self):
            super().__init__()
            self.figure = plt.Figure(figsize=(7, 5))
            self.canvas = FigureCanvasQTAgg(self.figure)
            lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0)
            lay.addWidget(NavigationToolbar2QT(self.canvas, self)); lay.addWidget(self.canvas)

        def axes(self):
            self.figure.clear(); return self.figure.add_subplot(111)

        def draw(self):
            self.figure.tight_layout(); self.canvas.draw()

    class App(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle(f'Kronos_GUI v{__version__} -- {db}')
            self.resize(1100, 720)
            self.model = None; self.saved = None
            central = QWidget(); self.setCentralWidget(central); root = QHBoxLayout(central)
            controls = QWidget(); form = QFormLayout(controls)
            self.ion = QComboBox(); self.charge = QComboBox(); self.target = QComboBox(); self.method = QComboBox()
            self.ldist = QComboBox(); self.ldist.addItems(LDISTS)
            self.energy = QLineEdit('1000'); self.erange_info = QLabel('')
            self.gauss = QCheckBox('Gaussian convolution'); self.fwhm = QLineEdit('10')
            self.emin = QLineEdit(); self.emin.setPlaceholderText('auto'); self.emax = QLineEdit(); self.emax.setPlaceholderText('auto')
            self.labels = QCheckBox('Label lines'); self.labels.setChecked(True); self.nlabels = QLineEdit('8')
            self.nres_plot = QCheckBox('Cross sections summed to n')
            b_cs = QPushButton('Plot Cross Sections'); b_run = QPushButton('Run Cascade')
            b_ref = QPushButton('Show Reference'); b_save = QPushButton('Save Spectrum + Line List...')
            self.ratios = QLabel(''); self.ratios.setWordWrap(True)
            self.status = QLabel(''); self.status.setWordWrap(True)
            for lab, w in (('Element', self.ion), ('Charge', self.charge), ('Target', self.target), ('Method', self.method),
                           ('l distribution (n-resolved)', self.ldist), ('Energy (eV/u)', self.energy), ('', self.erange_info)):
                form.addRow(lab, w)
            form.addRow(self.gauss); form.addRow('FWHM (eV)', self.fwhm)
            form.addRow('Photon energy min (eV)', self.emin); form.addRow('Photon energy max (eV)', self.emax)
            form.addRow(self.labels); form.addRow('Max labels', self.nlabels); form.addRow(self.nres_plot)
            for b in (b_cs, b_run, b_ref, b_save):
                form.addRow(b)
            form.addRow(QLabel('<b>Line ratios (n → 1 series, relative to α)</b>')); form.addRow(self.ratios)
            form.addRow(self.status)
            scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(controls); scroll.setFixedWidth(360)
            self.plot = Plot()
            split = QSplitter(); split.addWidget(scroll); split.addWidget(self.plot); split.setStretchFactor(1, 1)
            root.addWidget(split)
            self.ion.currentTextChanged.connect(self.set_ion); self.charge.currentTextChanged.connect(self.set_charge)
            self.target.currentTextChanged.connect(self.set_target); self.method.currentTextChanged.connect(self.set_method)
            self.energy.editingFinished.connect(self.pick_method)
            b_cs.clicked.connect(self.do_cs); b_run.clicked.connect(self.do_run)
            b_ref.clicked.connect(self.do_ref); b_save.clicked.connect(self.do_save)
            self.fill(self.ion, list(dict.fromkeys(systems.ion)))

        @staticmethod
        def fill(combo, items):
            combo.blockSignals(True); combo.clear(); combo.addItems([str(i) for i in items]); combo.blockSignals(False)
            if items:
                combo.setCurrentIndex(0); combo.currentTextChanged.emit(combo.currentText())

        def sel(self):
            return self.ion.currentText(), int(self.charge.currentText() or 0), self.target.currentText(), self.method.currentText()

        def set_ion(self, ion):
            self.fill(self.charge, sorted(set(systems[systems.ion == ion].charge)))

        def set_charge(self, q):
            if q:
                s = systems[(systems.ion == self.ion.currentText()) & (systems.charge == int(q))]
                self.fill(self.target, list(s.target))

        def set_target(self, t):
            if t:
                ion, q, _, _ = self.sel()
                self.fill(self.method, list(methods(db, ion, q, t)))
                self.pick_method()

        def pick_method(self):
            ion, q, t, _ = self.sel()
            try:
                m = preferred_method(db, ion, q, t, float(self.energy.text() or 1000))
            except (ValueError, FileNotFoundError):
                return
            if m and self.method.findText(m) >= 0:
                self.method.setCurrentIndex(self.method.findText(m))

        def set_method(self, m):
            if not m:
                return
            ion, q, t, _ = self.sel()
            self.ldist.setEnabled('nres' in m)
            try:
                df = read_cs(methods(db, ion, q, t)[m])
                e = df.iloc[:, 0]
                self.erange_info.setText(f'tabulated: {e.min():g} - {e.max():g} eV/u ({len(e)} energies)')
                rr = (system_dir(db, ion, q, t) / 'radrathe.dat').exists()
                self.status.setText('' if rr else 'No radrathe.dat in this folder: cross sections only (no cascade).')
            except Exception as ex:
                self.erange_info.setText(f'could not read file: {ex}')

        def model_for(self, run=True):
            ion, q, t, m = self.sel()
            return Kronos(ion, q, t, m, float(self.energy.text() or 1000), self.ldist.currentText(), db=db, verbose=False) \
                if run else None

        def do_cs(self):
            ion, q, t, m = self.sel()
            try:
                df = read_cs(methods(db, ion, q, t)[m])
            except Exception as ex:
                QMessageBox.critical(self, 'Kronos', str(ex)); return
            dfi = df.set_index(df.columns[0]); dfi = dfi[[c for c in dfi.columns if c.lower() != 'total']] * unit_factor(df.attrs['units'])
            if self.nres_plot.isChecked() and 'nres' not in m:
                dfi = dfi.T.groupby(lambda c: int(re.match(r'(\d+)', c).group(1)) if re.match(r'\d+', c) else 0).sum().T
                dfi.columns = [f'n={c}' for c in dfi.columns]
            keep = dfi.max().nlargest(10).index
            ax = self.plot.axes()
            dfi[keep].replace(0, np.nan).plot(ax=ax, logx=True, logy=True)
            ax.set_xlabel('collision energy (eV/u)'); ax.set_ylabel('cross section (cm$^2$)')
            ax.set_title(f'{ion}$^{{{q}+}}$ + {t}: {m} (top 10 states)'); ax.legend(fontsize=7, ncol=2)
            self.plot.draw()
            self.saved = ('cs', dfi, df.attrs['citation'])
            self.status.setText('Cross sections plotted.')

        def do_run(self):
            try:
                k = self.model_for()
            except Exception as ex:
                QMessageBox.critical(self, 'Kronos', f'Cascade failed:\n{ex}'); return
            self.model = k
            def num(w, d):
                try:
                    return float(w.text()) if w.text().strip() else d
                except ValueError:
                    return d
            er = None
            if self.emin.text().strip() or self.emax.text().strip():
                er = (num(self.emin, k.spectrum.Energy.min() - 50), num(self.emax, k.spectrum.Energy.max() + 50))
            ax = self.plot.axes()
            k.plot(fwhm=num(self.fwhm, 10) if self.gauss.isChecked() else None, erange=er,
                   label_lines=int(num(self.nlabels, 8)) if self.labels.isChecked() else 0, ax=ax, show=False)
            self.plot.draw()
            try:
                r = k.line_ratios()
                self.ratios.setText('\n'.join(f'{row.series} ({row.Energy:.1f} eV): {row.ratio:.3g}' for row in r.itertuples()))
            except ValueError:
                self.ratios.setText('(no n → 1 series lines)')
            self.saved = ('cascade', k, k.citation)
            self.status.setText(f'Cascade at {k.energy_used:g} eV/u (requested {k.energy_requested:g}); {len(k.spectrum)} lines.')

        def do_ref(self):
            if not self.saved:
                QMessageBox.information(self, 'Kronos', 'Plot cross sections or run a cascade first.'); return
            QMessageBox.information(self, 'Kronos - Reference', self.saved[2] or 'No reference in the file header.')

        def do_save(self):
            if not self.saved:
                QMessageBox.information(self, 'Kronos', 'Plot cross sections or run a cascade first.'); return
            if self.saved[0] == 'cascade':
                d = QFileDialog.getExistingDirectory(self, 'Folder for the spectrum and line list')
                if d:
                    paths = self.saved[1].save(d)
                    self.status.setText('Saved: ' + ', '.join(p.name for p in paths))
            else:
                p, _ = QFileDialog.getSaveFileName(self, 'Save cross sections', 'cross_sections.csv', 'CSV (*.csv)')
                if p:
                    self.saved[1].to_csv(p); self.status.setText(f'Saved {p}')

    app = QApplication.instance() or QApplication(sys.argv)
    w = App(); w.show()
    if not start:
        return w
    app.exec()


# ----------------------------------------------------------------------------------------------------------- terminal

def main(argv=None):
    p = argparse.ArgumentParser(description='Kronos_GUI v1: charge-exchange X-ray spectra from the Kronos CX database.')
    p.add_argument('ion', nargs='?', help='projectile element, e.g. O')
    p.add_argument('charge', nargs='?', type=int, help='projectile charge, e.g. 8')
    p.add_argument('target', nargs='?', help='target as in the database, e.g. H, He, H2O')
    p.add_argument('-E', '--energy', type=float, default=1000.0, help='collision energy in eV/u (nearest tabulated is used)')
    p.add_argument('-M', '--method', help='method (file name after _sec_), e.g. mclz, mclz_nres, aocc; default: preferred')
    p.add_argument('-D', '--ldist', default='low-energy', choices=LDISTS + list(LDIST_ALIASES),
                   help='l distribution for n-resolved files (even, landau-zener, separable = ACX 1, 3, 4)')
    p.add_argument('--db', help='database folder (default: CXDatabase next to this file, or $KRONOSDB)')
    p.add_argument('--top', type=int, default=20, help='number of lines to print')
    p.add_argument('--fwhm', type=float, help='Gaussian FWHM (eV) for the plot')
    p.add_argument('--plot', action='store_true', help='plot the spectrum')
    p.add_argument('--out', help='folder for the spectrum and line-list CSV files')
    p.add_argument('--list', action='store_true', help='list systems (all, or for ION [CHARGE])')
    p.add_argument('--gui', action='store_true', help='open the window (needs PySide6)')
    a = p.parse_args(argv)
    if a.gui:
        return gui(a.db)
    if a.list or not a.target:
        s = list_systems(a.db)
        if a.ion:
            s = s[s.ion == a.ion]
        if a.charge:
            s = s[s.charge == a.charge]
        with pd.option_context('display.max_rows', None, 'display.width', 200, 'display.max_colwidth', 80):
            print(s.to_string(index=False))
        return
    try:
        k = Kronos(a.ion, a.charge, a.target, a.method, a.energy, a.ldist, db=a.db)
    except (FileNotFoundError, ValueError) as ex:
        sys.exit(f'Kronos: {ex}')
    with pd.option_context('display.width', 200):
        print(k.line_list(top=a.top).to_string(formatters={'Intensity': '{:.3e}'.format, 'fraction': '{:.4f}'.format,
                                                           'Energy': '{:.2f}'.format}))
        try:
            print('\nn -> 1 line ratios:\n' + k.line_ratios().to_string(index=False))
        except ValueError:
            pass
    if a.out:
        print('saved:', *k.save(a.out))
    if a.plot:
        k.plot(fwhm=a.fwhm)


if __name__ == '__main__':
    main()
