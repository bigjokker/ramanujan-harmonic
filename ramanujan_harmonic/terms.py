"""Hypergeometric terms H_N(n) and the index derivative used in the post.

f(n) = H_N(n) / a^n is extended off the positive integers by the gamma
representation of H_N and the principal logarithm of a. For a < 0 and integer
n this still reproduces the ordinary real value a^{-n}. The principal
logarithm is log|a| + iπ. On every convergent row with a < 0 one has
Re(τ) = 1/2, so that branch is the shift τ -> τ - 1 in the identity
sum f'(n) = 2πi τ sum f(n).
"""

from __future__ import annotations

from fractions import Fraction

import mpmath as mp

from ramanujan_harmonic.table import Series


def _pack(parts: list[tuple[mp.mpf, mp.mpf]]) -> tuple[mp.mpf, mp.mpf]:
    log_u = mp.mpf(0)
    log_der = mp.mpf(0)
    for log_part, der_part in parts:
        log_u += log_part
        log_der += der_part
    return log_u, log_der


def _r2(n: int) -> tuple[mp.mpf, mp.mpf]:
    # 2^{-2n} (2n)! / (n!)^2
    log_r = -2 * n * mp.log(2) + mp.loggamma(2 * n + 1) - 2 * mp.loggamma(n + 1)
    der = -2 * mp.log(2) + 2 * mp.digamma(2 * n + 1) - 2 * mp.digamma(n + 1)
    return log_r, der


def _r3(n: int) -> tuple[mp.mpf, mp.mpf]:
    # 3^{-3n} (3n)! / (n!)^3
    log_r = -3 * n * mp.log(3) + mp.loggamma(3 * n + 1) - 3 * mp.loggamma(n + 1)
    der = -3 * mp.log(3) + 3 * mp.digamma(3 * n + 1) - 3 * mp.digamma(n + 1)
    return log_r, der


def _r4(n: int) -> tuple[mp.mpf, mp.mpf]:
    # 2^{-6n} (4n)! / ((2n)! (n!)^2)
    log_r = (
        -6 * n * mp.log(2)
        + mp.loggamma(4 * n + 1)
        - mp.loggamma(2 * n + 1)
        - 2 * mp.loggamma(n + 1)
    )
    der = (
        -6 * mp.log(2)
        + 4 * mp.digamma(4 * n + 1)
        - 2 * mp.digamma(2 * n + 1)
        - 2 * mp.digamma(n + 1)
    )
    return log_r, der


def _r6(n: int) -> tuple[mp.mpf, mp.mpf]:
    # 2^{-4n} 3^{-3n} (6n)! / ((3n)! (2n)! n!)
    log_r = (
        -4 * n * mp.log(2)
        - 3 * n * mp.log(3)
        + mp.loggamma(6 * n + 1)
        - mp.loggamma(3 * n + 1)
        - mp.loggamma(2 * n + 1)
        - mp.loggamma(n + 1)
    )
    der = (
        -4 * mp.log(2)
        - 3 * mp.log(3)
        + 6 * mp.digamma(6 * n + 1)
        - 3 * mp.digamma(3 * n + 1)
        - 2 * mp.digamma(2 * n + 1)
        - mp.digamma(n + 1)
    )
    return log_r, der


def log_u_and_derivative(level: int, n: int) -> tuple[mp.mpf, mp.mpf]:
    """Return log H_N(n) and d/dn log H_N(n). H_N(n) is positive."""
    if level == 1:
        return _pack([_r2(n), _r6(n)])
    if level == 2:
        return _pack([_r2(n), _r4(n)])
    if level == 3:
        return _pack([_r2(n), _r3(n)])
    if level == 4:
        return _pack([_r2(n), _r2(n), _r2(n)])
    raise ValueError(f"unknown level {level}")


def harmonic_log_derivative(level: int, n: int) -> mp.mpf:
    """The same derivative written with harmonic numbers, as in the post."""
    if n < 0:
        raise ValueError(n)
    if level == 4:
        # 3 * d/dn log R_2 = 6 (H_{2n} - H_n) - 6 log 2
        return 6 * (mp.harmonic(2 * n) - mp.harmonic(n)) - 6 * mp.log(2)
    if level == 2:
        # d/dn log(R_2 R_4) = 4 (H_{4n} - H_n) - 8 log 2
        return 4 * (mp.harmonic(4 * n) - mp.harmonic(n)) - 8 * mp.log(2)
    if level == 1:
        # d/dn log(R_2 R_6) = 3 (2 H_{6n} - H_{3n} - H_n) - log(2^6 * 3^3)
        return (
            3 * (2 * mp.harmonic(6 * n) - mp.harmonic(3 * n) - mp.harmonic(n))
            - mp.log(64 * 27)
        )
    if level == 3:
        # d/dn log(R_2 R_3) = 2(H_{2n}-H_n) + 3(H_{3n}-H_n) - 2 log 2 - 3 log 3
        return (
            2 * (mp.harmonic(2 * n) - mp.harmonic(n))
            + 3 * (mp.harmonic(3 * n) - mp.harmonic(n))
            - 2 * mp.log(2)
            - 3 * mp.log(3)
        )
    raise ValueError(f"unknown level {level}")


def _as_mp(value: Fraction) -> mp.mpf:
    return mp.mpf(value.numerator) / mp.mpf(value.denominator)


def term(series: Series, n: int) -> tuple[mp.mpc, mp.mpc]:
    """Return f(n) and f'(n), using the principal logarithm of a."""
    log_u, der_u = log_u_and_derivative(series.level, n)
    log_a = mp.log(_as_mp(series.a))
    f_n = mp.exp(log_u - n * log_a)
    f_prime = f_n * (der_u - log_a)
    return f_n, f_prime


def predicted_lambda(series: Series) -> mp.mpf:
    """Conjecture 1: λ = -2 π Im(τ), with Im(τ) = sqrt(tau_D) / tau_q."""
    im_sq = series.im_tau_squared()
    im_tau = mp.sqrt(_as_mp(im_sq))
    return -2 * mp.pi * im_tau


def direct_sums(series: Series, n_terms: int) -> tuple[mp.mpc, mp.mpc, mp.mpc, int]:
    """Partial sums of f, of f', and of (slope*n+intercept)*f.

    Stops early once terms are negligible compared with the working precision.
    Returns the sums and the last index included.
    """
    total_f = mp.mpc(0)
    total_fp = mp.mpc(0)
    total_pi = mp.mpc(0)
    cutoff = mp.mpf(10) ** (-(mp.mp.dps + 8))
    last = 0
    for n in range(n_terms):
        f_n, fp_n = term(series, n)
        weight = series.slope * n + series.intercept
        total_f += f_n
        total_fp += fp_n
        total_pi += weight * f_n
        last = n
        if n > 2 and abs(f_n) < cutoff and abs(fp_n) < cutoff:
            break
    return total_f, total_fp, total_pi, last


def nsum_sums(series: Series) -> tuple[mp.mpc, mp.mpc, mp.mpc]:
    """Sum a borderline series (|a| = 1) with mpmath's nsum."""

    def f_real(n):
        value, _ = term(series, int(n))
        return mp.re(value)

    def fp_real(n):
        _, deriv = term(series, int(n))
        return mp.re(deriv)

    def fp_imag(n):
        _, deriv = term(series, int(n))
        return mp.im(deriv)

    def pi_term(n):
        n_int = int(n)
        value, _ = term(series, n_int)
        return (series.slope * n_int + series.intercept) * mp.re(value)

    total_f = mp.nsum(f_real, [0, mp.inf])
    total_fp = mp.nsum(fp_real, [0, mp.inf]) + mp.j * mp.nsum(fp_imag, [0, mp.inf])
    total_pi = mp.nsum(pi_term, [0, mp.inf])
    return total_f, total_fp, total_pi


def terms_needed(series: Series) -> int:
    """Enough terms that a geometric tail is far below the working precision."""
    radius = abs(series.a)
    if radius <= 1:
        return 0
    # |term| decays like radius^{-n}. Ask for dps+12 guard digits.
    n_terms = int((mp.mp.dps + 12) * mp.log(10) / mp.log(radius)) + 5
    return max(n_terms, 8)


def sum_series(series: Series, n_terms: int | None = None):
    """Return (sum f, sum f', sum P f, terms_used).

    terms_used is None when the sum is performed by nsum.
    """
    if abs(series.a) <= 1:
        if not series.convergent:
            raise ValueError(f"series {series.id} diverges")
        total_f, total_fp, total_pi = nsum_sums(series)
        return total_f, total_fp, total_pi, None
    count = terms_needed(series) if n_terms is None else n_terms
    return direct_sums(series, count)
