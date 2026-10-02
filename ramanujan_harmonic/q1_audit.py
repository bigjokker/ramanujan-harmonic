"""Audit the MathOverflow answer to Q1 (series (2)).

The answer evaluates a Clausen/elliptic identity at Legendre's singular modulus

    k_3 = sqrt(2)/4 * (sqrt(3) - 1) = sin(π/12),

using K(k_3) = 3^{1/4} Γ(1/3)^3 / (2^{7/3} π) and K(k_3')/K(k_3) = sqrt(3).

This script checks each step that the answer uses, then checks that those steps
really produce series (2).

Run from the repository root:

    python -m ramanujan_harmonic.q1_audit
"""

from __future__ import annotations

from pathlib import Path

import mpmath as mp
import sympy as sp

from ramanujan_harmonic.console import configure_stdout

RESULTS = Path(__file__).resolve().parents[1] / "results"


def _record(rows: list[tuple[str, bool, str]], name: str, ok: bool, detail: str) -> None:
    rows.append((name, ok, detail))
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {name}")
    print(f"       {detail}")


def _symbolic(rows: list[tuple[str, bool, str]]) -> None:
    k = sp.sqrt(2) / 4 * (sp.sqrt(3) - 1)
    complement_square = sp.simplify(1 - k**2)
    product = sp.simplify(k**2 * (1 - k**2) - sp.Rational(1, 16))
    kk = sp.simplify(k * sp.sqrt(complement_square) - sp.Rational(1, 4))
    _record(
        rows,
        "singular modulus hits 256",
        product == 0 and kk == 0,
        "k_3^2 (1-k_3^2) = 1/16 and k_3 k_3' = 1/4, so the Nicholson sum is the 256-series",
    )

    # K^2 from the gamma evaluation quoted in the answer.
    gamma = sp.symbols("G", positive=True)
    pi = sp.pi
    K2 = sp.sqrt(3) * gamma**6 / (2 ** sp.Rational(14, 3) * pi**2)
    # At k_3 the identity says
    #   -(π/2) S = K^2 [ 1/sqrt(3) + (2/(3π)) log(k k'/4) ]
    # and log(k k'/4) = log(1/16) = -log(16).
    log_factor = -sp.log(16)
    S = sp.simplify(-(2 / pi) * K2 * (1 / sp.sqrt(3) + (2 / (3 * pi)) * log_factor))
    sigma = sp.sqrt(3) * gamma**6 / (4 * 2 ** sp.Rational(2, 3) * pi**4)
    claimed = (sp.log(256) - sp.sqrt(3) * pi) / 6 * sigma
    difference = sp.simplify(S - claimed)
    _record(
        rows,
        "identity at k_3 equals series (2)",
        difference == 0,
        "substituting K(k_3), K'/K = sqrt(3) and k k' = 1/4 recovers the closed form",
    )


def _numeric(rows: list[tuple[str, bool, str]]) -> None:
    mp.mp.dps = 40
    k = mp.sqrt(2) / 4 * (mp.sqrt(3) - 1)
    K = mp.ellipk(k**2)
    Kp = mp.ellipk(1 - k**2)
    K_formula = mp.power(3, mp.mpf("0.25")) * mp.gamma(mp.mpf(1) / 3) ** 3 / (
        mp.power(2, mp.mpf(7) / 3) * mp.pi
    )
    err_K = abs(K - K_formula)
    err_ratio = abs(Kp / K - mp.sqrt(3))
    _record(
        rows,
        "K(k_3) gamma value",
        err_K < mp.mpf("1e-30"),
        f"|ellipk(k_3^2) - formula| = {err_K}",
    )
    _record(
        rows,
        "K(k_3')/K(k_3) = sqrt(3)",
        err_ratio < mp.mpf("1e-30"),
        f"|K'/K - sqrt(3)| = {err_ratio}",
    )

    def nicholson_lhs(modulus: mp.mpf) -> mp.mpf:
        # -(π/2) ∑ binom(2n,n)^3 (H_{2n}-H_n) (k^2(1-k^2))^n / 16^n
        # binom(2n,n)^3 / 16^n grows like 4^n, so the ratio is 4 k^2 (1-k^2).
        total = mp.mpf(0)
        weight = modulus**2 * (1 - modulus**2)
        ratio = 4 * weight
        n_terms = int((mp.mp.dps + 8) * mp.log(10) / abs(mp.log(ratio))) + 10
        for n in range(1, n_terms + 1):
            term = (
                mp.binomial(2 * n, n) ** 3
                * (mp.harmonic(2 * n) - mp.harmonic(n))
                * mp.power(weight, n)
                / mp.power(16, n)
            )
            total += term
        return -(mp.pi / 2) * total

    def nicholson_rhs(modulus: mp.mpf) -> mp.mpf:
        KK = mp.ellipk(modulus**2)
        KKp = mp.ellipk(1 - modulus**2)
        return KK * KKp / 3 + (2 / (3 * mp.pi)) * mp.log(modulus * mp.sqrt(1 - modulus**2) / 4) * KK**2

    for sample, digits in ((mp.mpf("0.2"), 30), (mp.mpf("0.5"), 30), (k, 30)):
        gap = abs(nicholson_lhs(sample) - nicholson_rhs(sample))
        _record(
            rows,
            f"Nicholson identity at k={sample}",
            gap < mp.mpf(10) ** (-digits),
            f"|LHS - RHS| = {gap}",
        )

    def base_sum(kind: str, modulus: mp.mpf, n_terms: int = 80) -> mp.mpf:
        total = mp.mpf(0)
        for n in range(1, n_terms + 1):
            harmonic = mp.harmonic(2 * n) if kind == "H2n" else mp.harmonic(n)
            total += (
                mp.binomial(2 * n, n) ** 2
                * harmonic
                * mp.power(modulus**2, n)
                / mp.power(16, n)
            )
        return total

    sample = mp.mpf("0.4")
    Ks = mp.ellipk(sample**2)
    Ksp = mp.ellipk(1 - sample**2)
    rhs1 = Ksp + (1 / mp.pi) * Ks * mp.log(sample**2 / (16 * (1 - sample**2)))
    rhs2 = Ksp / 2 + (1 / mp.pi) * Ks * mp.log(sample / (4 * (1 - sample**2)))
    gap1 = abs(base_sum("Hn", sample) - rhs1)
    gap2 = abs(base_sum("H2n", sample) - rhs2)
    _record(
        rows,
        "generating function (1) at k=0.4",
        gap1 < mp.mpf("1e-30"),
        f"|series - elliptic form| = {gap1}",
    )
    _record(
        rows,
        "generating function (2) at k=0.4",
        gap2 < mp.mpf("1e-30"),
        f"|series - elliptic form| = {gap2}",
    )

    sigma = (
        mp.sqrt(3)
        * mp.gamma(mp.mpf(1) / 3) ** 6
        / (4 * mp.power(2, mp.mpf(2) / 3) * mp.pi**4)
    )
    claimed_S = (mp.log(256) - mp.sqrt(3) * mp.pi) / 6 * sigma
    partial_f = mp.mpf(0)
    partial_S = mp.mpf(0)
    for n in range(0, 60):
        f_n = mp.binomial(2 * n, n) ** 3 / mp.power(256, n)
        partial_f += f_n
        if n:
            partial_S += (mp.harmonic(2 * n) - mp.harmonic(n)) * f_n
    _record(
        rows,
        "sum f equals the gamma product",
        abs(partial_f - sigma) < mp.mpf("1e-30"),
        f"|partial sum - closed form| = {abs(partial_f - sigma)}",
    )
    _record(
        rows,
        "series (2) equals its closed form",
        abs(partial_S - claimed_S) < mp.mpf("1e-30"),
        f"|partial sum - closed form| = {abs(partial_S - claimed_S)}",
    )
    lam = (partial_S * 6 / partial_f) - mp.log(256)
    # From 6 S = (λ + log 256) sum f, with the post's sign convention λ = sum f' / sum f.
    gap_lambda = abs(lam - (-mp.sqrt(3) * mp.pi))
    _record(
        rows,
        "λ = -sqrt(3) π for this f",
        gap_lambda < mp.mpf("1e-30"),
        f"|λ + sqrt(3) π| = {gap_lambda}",
    )


def main() -> int:
    configure_stdout()
    rows: list[tuple[str, bool, str]] = []
    print("Q1 audit: User-Refolio's answer to series (2)")
    print(f"working precision for the numeric checks: 40 decimal places")
    _symbolic(rows)
    _numeric(rows)
    failed = [name for name, ok, _ in rows if not ok]
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "q1-audit.txt"
    lines = ["Q1 audit of https://mathoverflow.net/a/514778", ""]
    for name, ok, detail in rows:
        lines.append(f"{'PASS' if ok else 'FAIL'}  {name}")
        lines.append(f"      {detail}")
    lines.append("")
    if failed:
        lines.append("Result: the answer does not check out.")
        lines.append("Failed: " + ", ".join(failed))
    else:
        series_detail = next(
            detail for name, _, detail in rows if name.startswith("series (2)")
        )
        lines.append("Result: the answer checks out.")
        lines.append(
            "Legendre's modulus gives the 256^k series. K(k_3) and K'/K = sqrt(3) match."
        )
        lines.append(
            "Nicholson's identity is checked numerically at k=0.2, 0.5 and k_3. "
            "This script does not prove it."
        )
        lines.append(
            "Substituting that identity at k_3 produces the closed form of series (2). "
            + series_detail
        )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    start = next(i for i, line in enumerate(lines) if line.startswith("Result:"))
    for line in lines[start:]:
        print(line)
    print(f"wrote {out}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
