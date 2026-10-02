"""Test Q2 and Q3 on the Cohen–Guillera table.

The post answers Q2 by Conjecture 1. The identity checked here is stronger:

    sum f'(n) = 2πi τ sum f(n),

with τ replaced by τ - 1 when a < 0. That shift is the principal branch of
log a. Its real part is Conjecture 1, so λ = -2π Im(τ).

The same run checks J_N(τ) = a from q-expansions, and checks
sum (slope n + intercept) f(n) = sqrt(k)/π.

Divergent rows 37–44 are stored in the table but are not classical sums, so
they are skipped unless --include-divergent is passed.

Run from the repository root:

    python -m ramanujan_harmonic.conjecture
    python -m ramanujan_harmonic.conjecture --ids 35 23 11 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import mpmath as mp

from ramanujan_harmonic.console import configure_stdout
from ramanujan_harmonic.modular import hauptmodul_error
from ramanujan_harmonic.table import BY_ID, ENTRIES, Series
from ramanujan_harmonic.terms import harmonic_log_derivative, log_u_and_derivative, sum_series

RESULTS = Path(__file__).resolve().parents[1] / "results"


def check_harmonic_derivatives(dps: int = 30) -> None:
    """Lock the digamma derivative to the harmonic formulas in the post."""
    mp.mp.dps = dps
    for level in (1, 2, 3, 4):
        for n in (0, 1, 2, 5, 12):
            _, digamma_der = log_u_and_derivative(level, n)
            gap = abs(digamma_der - harmonic_log_derivative(level, n))
            if gap > mp.mpf(10) ** (-(dps - 6)):
                raise RuntimeError(f"level {level} n={n} derivative mismatch {gap}")


def _relative(got, expected) -> mp.mpf:
    scale = max(abs(expected), mp.mpf(1))
    return abs(got - expected) / scale


def branch_factor(series: Series) -> mp.mpc:
    """2πi τ, or 2πi(τ - 1) when a < 0.

    The shift is the principal logarithm: every convergent negative-a row has
    Re(τ) = 1/2, and log(-|a|) = log|a| + iπ = log|a| + 2πi(-1/2).
    """
    real = mp.mpf(series.tau_p) / series.tau_q
    imag = mp.sqrt(mp.mpf(series.tau_D)) / series.tau_q
    tau = mp.mpc(real, imag)
    if series.a < 0:
        tau -= 1
    return 2 * mp.pi * mp.j * tau


def check_one(series: Series, dps: int) -> dict:
    modular_err = hauptmodul_error(series)
    mp.mp.dps = dps
    total_f, total_fp, total_pi, used = sum_series(series)
    if series.k < 0:
        root = mp.sqrt(abs(mp.mpf(series.k.numerator) / mp.mpf(series.k.denominator)))
        target_pi = mp.sign(series.k) * root / mp.pi
    else:
        target_pi = mp.sqrt(mp.mpf(series.k.numerator) / mp.mpf(series.k.denominator)) / mp.pi
    factor = branch_factor(series)
    predicted = factor * total_f
    return {
        "id": series.id,
        "level": series.level,
        "label": series.label,
        "tau": series.tau_text(),
        "branch": "τ-1" if series.a < 0 else "τ",
        "a": series.a,
        "terms": used,
        "modular_err": modular_err,
        "pi_err": _relative(mp.re(total_pi), target_pi),
        "ratio": total_fp / total_f,
        "factor": factor,
        "complex_err": _relative(total_fp, predicted),
    }


def _format(row: dict, digits_needed: int) -> list[str]:
    ok = (
        row["modular_err"] < mp.mpf("1e-20")
        and row["pi_err"] < mp.mpf(10) ** (-digits_needed)
        and row["complex_err"] < mp.mpf(10) ** (-digits_needed)
    )
    label = f"  {row['label']}" if row["label"] else ""
    terms = "nsum" if row["terms"] is None else f"n<={row['terms']}"
    lines = [
        f"[{'PASS' if ok else 'FAIL'}] #{row['id']} level {row['level']}  τ={row['tau']}  a={row['a']}  ({terms}){label}",
        f"       J_{row['level']}(τ) relative error = {row['modular_err']}",
        f"       sqrt(k)/π relative error = {row['pi_err']}",
        f"       sum f' / sum f = {row['ratio']}",
        f"       2πi {row['branch']} = {row['factor']}",
        f"       complex relative error = {row['complex_err']}",
    ]
    return lines


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Test Conjecture 1 on the Cohen–Guillera table.")
    parser.add_argument("--ids", type=int, nargs="*", help="series numbers to test (default: all 36 convergent)")
    parser.add_argument("--dps", type=int, default=25, help="working decimal places (default 25)")
    parser.add_argument(
        "--digits",
        type=int,
        default=8,
        help="required correct digits in the relative error (default 8)",
    )
    parser.add_argument(
        "--include-divergent",
        action="store_true",
        help="also attempt rows 37–44; their ordinary sums diverge",
    )
    args = parser.parse_args(argv)
    configure_stdout()

    print("Checking harmonic-number derivatives against the digamma form...")
    check_harmonic_derivatives()
    print("Derivative convention matches the post.")
    print()

    if args.ids:
        chosen = [BY_ID[i] for i in args.ids]
    else:
        chosen = [entry for entry in ENTRIES if entry.convergent or args.include_divergent]

    rows = []
    lines = [
        "sum f'(n) = 2πi τ sum f(n), with τ replaced by τ-1 when a<0",
        "Real part: Re(sum f') = -2π Im(τ) sum f   (Conjecture 1)",
        f"dps={args.dps}; a row fails if the complex or sqrt(k)/π error loses {args.digits} digits,",
        "or if J_N(τ) misses a by more than 1e-20.",
        "",
    ]
    print(lines[0])
    print(lines[1])
    failed = []
    for series in chosen:
        if not series.convergent and not args.include_divergent:
            continue
        print(f"summing #{series.id} ...", flush=True)
        try:
            row = check_one(series, args.dps)
        except Exception as exc:  # a divergent sum should be reported, not hidden
            failed.append(series.id)
            message = f"[FAIL] #{series.id} raised {type(exc).__name__}: {exc}"
            print(message)
            lines.append(message)
            continue
        block = _format(row, args.digits)
        print("\n".join(block))
        lines.extend(block)
        ok = (
            row["modular_err"] < mp.mpf("1e-20")
            and row["pi_err"] < mp.mpf(10) ** (-args.digits)
            and row["complex_err"] < mp.mpf(10) ** (-args.digits)
        )
        if not ok:
            failed.append(series.id)
        rows.append(row)

    lines.append("")
    if not rows and failed:
        summary = "No series was summed."
    elif failed:
        summary = "Failed for series " + ", ".join(str(i) for i in failed) + "."
    else:
        worst_complex = max(row["complex_err"] for row in rows)
        worst_modular = max(row["modular_err"] for row in rows)
        worst_pi = max(row["pi_err"] for row in rows)
        summary = (
            f"All {len(rows)} tested series satisfy sum f' = 2πi τ_branch sum f, "
            f"J_N(τ) = a, and sum P(n) f(n) = sqrt(k)/π. "
            f"Worst relative errors at dps={args.dps}: "
            f"complex {worst_complex}, Hauptmodul {worst_modular}, sqrt(k)/π {worst_pi}."
        )
    lines.append(summary)
    print()
    print(summary)

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "conjecture.txt"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
