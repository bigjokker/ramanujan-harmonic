"""Recover each table entry's a from its CM point τ.

a is defined to be the Hauptmodul J_N(τ). This module recomputes J_N from
Eisenstein q-expansions and does not use the hypergeometric term in terms.py.

    level 1:  J_1 = j/1728 = E_4^3 / (E_4^3 - E_6^2)
    level 2, 3: J_N = F_2^4 / (F_2^4 - F_4^2)
    level 4:  J_4 = F_2^2 / (F_2^2 - G_2^2)

with F_2, F_4 and G_2 as in Cohen–Guillera, arXiv:2101.12592, Section 2.
"""

from __future__ import annotations

import mpmath as mp

from ramanujan_harmonic.table import Series

_SIGMA: dict[int, list[int]] = {}
_SIGMA_N = 0


def _sigma_table(power: int, nmax: int) -> list[int]:
    """sigma_power(n) for n = 0..nmax. Index 0 is unused."""
    global _SIGMA_N
    if nmax > _SIGMA_N:
        for key in (1, 3, 5):
            values = [0] * (nmax + 1)
            for d in range(1, nmax + 1):
                term = d**key
                for multiple in range(d, nmax + 1, d):
                    values[multiple] += term
            _SIGMA[key] = values
        _SIGMA_N = nmax
    return _SIGMA[power]


def _eisenstein(constant: int, power: int, q: mp.mpc, nmax: int) -> mp.mpc:
    sigma = _sigma_table(power, nmax)
    total = mp.mpc(0)
    qq = q
    for n in range(1, nmax + 1):
        total += sigma[n] * qq
        qq *= q
    return 1 + constant * total


def _tau(series: Series) -> mp.mpc:
    return mp.mpc(series.tau_p, mp.sqrt(series.tau_D)) / series.tau_q


def hauptmodul(series: Series, nmax: int = 80) -> mp.mpc:
    """J_N(τ) for this row. Caller sets the working precision."""
    return hauptmodul_at(series.level, _tau(series), nmax)


def hauptmodul_at(level: int, tau: mp.mpc, nmax: int = 80) -> mp.mpc:
    """J_N(τ) at any τ in the upper half-plane. Caller sets the working precision."""
    q = mp.exp(2 * mp.pi * mp.j * tau)
    if level == 1:
        e4 = _eisenstein(240, 3, q, nmax)
        e6 = _eisenstein(-504, 5, q, nmax)
        discriminant = e4**3 - e6**2
        return e4**3 / discriminant
    if level in (2, 3):
        e2 = _eisenstein(-24, 1, q, nmax)
        e2n = _eisenstein(-24, 1, q**level, nmax)
        e4 = _eisenstein(240, 3, q, nmax)
        e4n = _eisenstein(240, 3, q**level, nmax)
        f2 = (level * e2n - e2) / (level - 1)
        f4 = (level**2 * e4n - e4) / (level**2 - 1)
        return f2**4 / (f2**4 - f4**2)
    if level == 4:
        e2 = _eisenstein(-24, 1, q, nmax)
        e22 = _eisenstein(-24, 1, q**2, nmax)
        e24 = _eisenstein(-24, 1, q**4, nmax)
        f2 = (4 * e24 - e2) / 3
        g2 = 4 * e24 - 4 * e22 + e2
        return f2**2 / (f2**2 - g2**2)
    raise ValueError(f"unknown level {level}")


def hauptmodul_error(series: Series) -> mp.mpf:
    """Relative |J_N(τ) - a| / |a|, computed above the cancellation in j."""
    im = mp.sqrt(mp.mpf(series.tau_D)) / series.tau_q
    # j ~ 1/q, and q = exp(2πiτ), so large Im(τ) cancels leading 1's in E4 and E6.
    digits = int(mp.ceil(3 * im + 40))
    with mp.workdps(max(digits, 50)):
        value = hauptmodul(series)
        target = mp.mpf(series.a.numerator) / mp.mpf(series.a.denominator)
        return abs(value - target) / abs(target)
