"""Exact proofs of the modular identities in step (2d) of docs/proof.md.

Each identity is an equality of two holomorphic modular forms of the same
weight k and character on Γ0(N). By Sturm's theorem (Ono, The Web of
Modularity, Thm 2.58, p. 40) two such forms are equal
once their q-expansions agree through q^B, B = floor(k m / 12) with
m = [SL2(Z) : Γ0(N)]. This script computes every q-expansion with exact integer
or rational arithmetic and compares the coefficients through q^B. It also
compares them further, up to q^(PREC-1), but that is only a sanity check.

Every form records why it lies in its space:

- Eisenstein: E4(τ) and E6(τ) are in M4 and M6 of SL2(Z) (Diamond–Shurman §1.1,
  pp. 5–6), and
  E2(τ) - N E2(Nτ) is in M2(Γ0(N)) (Diamond–Shurman §1.2, p. 18).
- Eta quotient: the Gordon–Hughes–Newman conditions (Ono Thm 1.64) and
  Ligozat's cusp orders (Ono Thm 1.65) are checked here, not assumed.
- Theta series: a(q) = sum q^(m^2+mn+n^2) is in M1(Γ0(3), (-3/.))
  (Diamond–Shurman Thm 4.11.3, p. 159, with N = 1 and χ trivial). It is cited,
  not checked.
- f(Mτ) for f in M_k(Γ0(N), χ) is in M_k(Γ0(NM), χ), proved in docs/proof.md,
  [C6]. Sums and products follow the usual rules.

The F2, F4 and G2 below are the Eisenstein combinations in modular.py.

Run from the repository root:

    python -m ramanujan_harmonic.sturm

Writes results/sturm.txt. Exits 1 if any check fails, including the negative
controls, which confirm that a false identity fails and that bad eta quotients
are rejected.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path

from ramanujan_harmonic.console import configure_stdout

PREC = 60  # coefficients q^0 .. q^(PREC-1); every Sturm bound below is far smaller
RESULTS = Path(__file__).resolve().parents[1] / "results"


# --- exact power series -----------------------------------------------------

def _mul(a: list, b: list) -> list:
    c = [0] * PREC
    for i, x in enumerate(a):
        if x:
            for j in range(PREC - i):
                c[i + j] += x * b[j]
    return c


def _inv(a: list) -> list:
    if a[0] != 1:
        raise ValueError("series must start with 1")
    r = [Fraction(0)] * PREC
    r[0] = Fraction(1)
    for n in range(1, PREC):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1))
    return r


def _pow(a: list, e: int) -> list:
    r = [1] + [0] * (PREC - 1)
    for _ in range(e):
        r = _mul(r, a)
    return r


def _sigma(k: int, n: int) -> int:
    return sum(d**k for d in range(1, n + 1) if n % d == 0)


def _index(N: int) -> int:
    """[SL2(Z) : Γ0(N)] = N prod_{p | N} (1 + 1/p)."""
    m, n, p = N, N, 2
    while p * p <= n:
        if n % p == 0:
            m = m // p * (p + 1)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        m = m // n * (n + 1)
    return m


def _squarefree(x: Fraction) -> int:
    """Squarefree integer D with x = D * (rational square)."""
    num, den = x.numerator * x.denominator, 1  # x and num/den^2 differ by a square
    sign = -1 if num < 0 else 1
    num = abs(num)
    d, p = 1, 2
    while p * p <= num:
        while num % (p * p) == 0:
            num //= p * p
        if num % p == 0:
            d *= p
            num //= p
        p += 1
    return sign * d * num


# --- modular forms with provenance ---------------------------------------------

@dataclass
class Form:
    coeffs: list
    weight: int
    level: int
    char: int  # squarefree D: the character is the Kronecker symbol (D/.)
    why: str

    def _same_space(self, other: "Form") -> int:
        if self.weight != other.weight or self.char != other.char:
            raise ValueError(f"cannot add {self.why} and {other.why}: different spaces")
        return self.level * other.level // gcd(self.level, other.level)

    def __add__(self, other: "Form") -> "Form":
        level = self._same_space(other)
        return Form([x + y for x, y in zip(self.coeffs, other.coeffs)], self.weight, level, self.char, "sum")

    def __sub__(self, other: "Form") -> "Form":
        return self + other.scale(-1)

    def __mul__(self, other: "Form") -> "Form":
        level = self.level * other.level // gcd(self.level, other.level)
        char = _squarefree(Fraction(self.char * other.char))
        return Form(_mul(self.coeffs, other.coeffs), self.weight + other.weight, level, char, "product")

    def __pow__(self, e: int) -> "Form":
        result = self
        for _ in range(e - 1):
            result = result * self
        return result

    def scale(self, c) -> "Form":
        return Form([c * x for x in self.coeffs], self.weight, self.level, self.char, self.why)

    def at_multiple(self, M: int) -> "Form":
        """f(Mτ), which lies in M_k(Γ0(N M), χ)."""
        coeffs = [0] * PREC
        for n in range(0, PREC, M):
            coeffs[n] = self.coeffs[n // M]
        return Form(coeffs, self.weight, self.level * M, self.char, f"({self.why})(M={M})")


def e4(M: int = 1) -> Form:
    """E4(Mτ) in M4(Γ0(M))."""
    c = [0] * PREC
    c[0] = 1
    for n in range(1, PREC):
        if M * n < PREC:
            c[M * n] += 240 * _sigma(3, n)
    return Form(c, 4, M, 1, f"E4({M}τ)")


def e6() -> Form:
    """E6 in M6(SL2(Z))."""
    c = [0] * PREC
    c[0] = 1
    for n in range(1, PREC):
        c[n] = -504 * _sigma(5, n)
    return Form(c, 6, 1, 1, "E6(τ)")


def e2_minus(N: int) -> Form:
    """E2(τ) - N E2(Nτ) in M2(Γ0(N))."""
    c = [0] * PREC
    c[0] = 1 - N
    for n in range(1, PREC):
        c[n] += -24 * _sigma(1, n)
        if N * n < PREC:
            c[N * n] += 24 * N * _sigma(1, n)
    return Form(c, 2, N, 1, f"E2(τ) - {N} E2({N}τ)")


def _euler(d: int) -> list:
    """prod_{n>=1} (1 - q^(d n))."""
    r = [1] + [0] * (PREC - 1)
    for n in range(1, PREC):
        if d * n >= PREC:
            break
        factor = [0] * PREC
        factor[0], factor[d * n] = 1, -1
        r = _mul(r, factor)
    return r


def eta_quotient(N: int, exps: dict[int, int], const: int = 1, allow_poles: bool = False) -> Form:
    """const * prod eta(δτ)^r_δ, after checking it lies in M_k(Γ0(N), χ).

    Gordon–Hughes–Newman: if k = sum r_δ / 2 is an integer, sum δ r_δ ≡ 0 and
    sum (N/δ) r_δ ≡ 0 (mod 24), then f transforms under Γ0(N) with weight k
    and character ((-1)^k s / .), s = prod δ^r_δ. Ligozat: the order at the cusp
    c/d (d | N) is (N/24) sum gcd(d,δ)^2 r_δ / (gcd(d, N/d) d δ), and f is
    holomorphic there iff this is >= 0.
    """
    name = f"{const}·" * (const != 1) + "".join(f"η({d}τ)^{r}" for d, r in exps.items())
    if any(N % d for d in exps):
        raise ValueError(f"{name}: every δ must divide N={N}")
    total = sum(exps.values())
    if total % 2:
        raise ValueError(f"{name}: odd weight sum")
    k = total // 2
    if sum(d * r for d, r in exps.items()) % 24 or sum((N // d) * r for d, r in exps.items()) % 24:
        raise ValueError(f"{name}: fails the Gordon–Hughes–Newman congruences at level {N}")
    orders = {}
    for d in (d for d in range(1, N + 1) if N % d == 0):
        order = Fraction(N, 24) * sum(
            Fraction(gcd(d, delta) ** 2 * r, gcd(d, N // d) * d * delta) for delta, r in exps.items()
        )
        if order < 0 and not allow_poles:
            raise ValueError(f"{name}: pole at cusps with denominator {d}")
        orders[d] = order
    s = Fraction(1)
    for d, r in exps.items():
        s *= Fraction(d) ** r
    char = _squarefree((-1) ** k * s)
    offset24 = sum(d * r for d, r in exps.items())
    series = [1] + [0] * (PREC - 1)
    for d, r in exps.items():
        e = _euler(d)
        series = _mul(series, _pow(e, r) if r >= 0 else _pow(_inv(e), -r))
    shift = offset24 // 24
    coeffs = [0] * shift + [const * x for x in series[: PREC - shift]]
    cusps = ", ".join(
        f"{'∞' if d == N else '0' if d == 1 else f'1/{d}'}: {o}" for d, o in sorted(orders.items(), key=lambda x: -x[0])
    )
    space = f"M{k}(Γ0({N})" + (f", ({char}/.))" if char != 1 else ")")
    if allow_poles:
        space = f"A{k}(Γ0({N}))"  # weight-k modular function, poles allowed at cusps
    record(f"eta quotient {name} in {space}",
           True,
           f"Σr_δ = {total}, Σδr_δ = {offset24} ≡ 0, Σ(N/δ)r_δ = {sum((N // d) * r for d, r in exps.items())} ≡ 0 (mod 24); "
           f"Ligozat orders at the cusps (c/d grouped by d) {cusps}")
    form = Form(coeffs, k, N, char, name)
    form.orders = orders
    return form


def theta_a() -> Form:
    """a(q) = sum q^(m^2+mn+n^2), Berndt–Bhargava–Garvan (2.2). Cited: M1(Γ0(3), (-3/.))."""
    c = [0] * PREC
    r = isqrt(4 * PREC) + 2
    for m in range(-r, r + 1):
        for n in range(-r, r + 1):
            e = m * m + m * n + n * n
            if e < PREC:
                c[e] += 1
    return Form(c, 1, 3, -3, "a(q) [theta series; Diamond–Shurman Thm 4.11.3]")


def series_b() -> list:
    """b(q) = sum ω^(m-n) q^(m^2+mn+n^2), Berndt–Bhargava–Garvan (2.3), as an exact series."""
    c = [Fraction(0)] * PREC
    r = isqrt(4 * PREC) + 2
    for m in range(-r, r + 1):
        for n in range(-r, r + 1):
            e = m * m + m * n + n * n
            if e < PREC:
                c[e] += 1 if (m - n) % 3 == 0 else Fraction(-1, 2)  # Re ω^j; imaginary parts cancel
    return c


# --- the identities -----------------------------------------------------------------

rows: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str) -> None:
    rows.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    print(f"       {detail}")


def compare(lhs: Form, rhs: Form) -> tuple[bool, str]:
    """Sturm comparison of two forms. Returns (holds, detail)."""
    if (lhs.weight, lhs.char) != (rhs.weight, rhs.char):
        return False, f"different spaces: weight {lhs.weight}/{rhs.weight}, character {lhs.char}/{rhs.char}"
    N = lhs.level * rhs.level // gcd(lhs.level, rhs.level)
    bound = lhs.weight * _index(N) // 12
    upto = next((n for n in range(PREC) if lhs.coeffs[n] != rhs.coeffs[n]), None)
    ok = upto is None or upto > bound
    space = f"M{lhs.weight}(Γ0({N})" + (f", ({lhs.char}/.))" if lhs.char != 1 else ")")
    agree = PREC - 1 if upto is None else upto - 1
    detail = (f"{space}; Sturm bound: agreement needed through q^{bound}; "
              f"compared q^0..q^{PREC - 1} ({PREC} coefficients), agree through q^{agree}")
    return ok, detail


def identity(name: str, lhs: Form, rhs: Form) -> None:
    record(name, *compare(lhs, rhs))


def _cusp_count(N: int) -> int:
    """Number of cusps of Γ0(N): sum over d | N of φ(gcd(d, N/d)) (Diamond–Shurman §3.8)."""
    def phi(n: int) -> int:
        return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)
    return sum(phi(gcd(d, N // d)) for d in range(1, N + 1) if N % d == 0)


def step3_hauptmoduln() -> None:
    """docs/proof.md (3a) and (3b): the fundamental-domain inequality and the degree of t and λ."""
    # (3a): the integer bounds, proved in the text for all c, d; checked here on a box as a control.
    worst = []
    for N in (1, 2, 3, 4):
        for c in range(-200, 201):
            for d in range(-200, 201):
                if c and c % N == 0 and gcd(c, d) == 1 and Fraction(c * c, N) - abs(c * d) + d * d < 1:
                    worst.append((N, "Γ0", c, d))
                if N > 1 and c and gcd(c, N * d) == 1 and c * c - N * abs(c * d) + N * d * d < 1:
                    worst.append((N, "W", c, d))
    record("(3a) fundamental-domain bounds on |c|, |d| <= 200 (proved in the text for all c, d)",
           not worst, "c²/N - |cd| + d² >= 1 on Γ0(N); c² - N|cd| + N d² >= 1 on the Fricke coset"
           + ("" if not worst else f"; violations {worst[:5]}"))

    # (3b): each function has weight 0, trivial character, one cusp per divisor d of N,
    # and exactly one pole, which is simple.
    for N, exps, const, label in (
        (2, {2: 24, 1: -24}, 1, "t = (η(2τ)/η(τ))^24"),
        (3, {3: 12, 1: -12}, 1, "t = (η(3τ)/η(τ))^12"),
        (4, {1: 8, 4: 16, 2: -24}, 16, "λ = 16 η(τ)^8 η(4τ)^16 / η(2τ)^24"),
    ):
        f = eta_quotient(N, exps, const, allow_poles=True)
        divisors = [d for d in range(1, N + 1) if N % d == 0]
        one_cusp_per_d = _cusp_count(N) == len(divisors)
        poles = {d: o for d, o in f.orders.items() if o < 0}
        ok = f.weight == 0 and f.char == 1 and one_cusp_per_d and list(poles.values()) == [-1] \
            and sum(f.orders.values()) == 0
        record(f"(3b) {label}: one simple pole on X0({N})", ok,
               f"weight {f.weight}, character {f.char}, {_cusp_count(N)} cusps, "
               "orders by cusp denominator d: " + ", ".join(f"d={d}: {o}" for d, o in sorted(f.orders.items(), key=lambda x: -x[0])) + " (d = N is ∞, d = 1 is 0)")


def negative_controls(A: Form, B: Form, F2: Form, weight2: Form) -> None:
    """The checks above must be able to fail. Each control passes only if it is caught."""
    holds, detail = compare(F2 * F2, A + B.scale(63))
    record("control: a false identity fails", not holds, f"F2^2 = A + 63 B: {detail}")
    for N, exps, why in (
        (3, {1: 3, 3: -1}, "η(τ)^3/η(3τ) at level 3 fails Σ(N/δ)r_δ ≡ 0 (mod 24)"),
        (4, {1: -8, 2: 4}, "η(2τ)^4/η(τ)^8 has a pole at the cusp 0"),
        (2, {1: 24, 2: -24}, "(η(τ)/η(2τ))^24 has a pole at the cusp ∞"),
    ):
        before = len(rows)
        try:
            eta_quotient(N, exps)
            rejected, reason = False, "accepted"
        except ValueError as err:
            rejected, reason = True, str(err)
        del rows[before:]  # drop the row an accepted quotient would have added
        record(f"control: rejected, {why}", rejected, f"level {N}: {reason}")
    try:
        _ = A + weight2
        mixed = False
    except ValueError:
        mixed = True
    record("control: adding forms of different weight is refused", mixed, "M4(Γ0(2)) + M2(Γ0(4))")


def main() -> int:
    configure_stdout()

    # Theta functions at nome q = e^{2πiτ}, as eta quotients on Γ0(4).
    th3 = eta_quotient(4, {2: 20, 1: -8, 4: -8})        # θ3^4 = η(2τ)^20 / (η(τ)^8 η(4τ)^8)
    th4 = eta_quotient(4, {1: 8, 2: -4})                 # θ4^4 = η(τ)^8 / η(2τ)^4
    th2 = eta_quotient(4, {4: 8, 2: -4}, const=16)       # θ2^4 = 16 η(4τ)^8 / η(2τ)^4

    # The theta sums themselves. Ono Thm 1.60 proves they equal these eta quotients;
    # the comparison here only checks the transcription.
    direct = {3: [0] * PREC, 4: [0] * PREC}
    for n in range(-isqrt(PREC) - 1, isqrt(PREC) + 2):
        if n * n < PREC:
            direct[3][n * n] += 1
            direct[4][n * n] += (-1) ** abs(n)
    half = [0] * PREC  # sum q^(n^2+n) over all n, so θ2 = q^(1/4) * half and θ2^4 = q * half^4
    for n in range(-isqrt(PREC) - 1, isqrt(PREC) + 2):
        if n * n + n < PREC:
            half[n * n + n] += 1
    th2_direct = [0] + _pow(half, 4)[: PREC - 1]
    record("θ3^4, θ4^4, θ2^4 eta quotients match the theta sums (Ono Thm 1.60; transcription check)",
           th3.coeffs == _pow(direct[3], 4) and th4.coeffs == _pow(direct[4], 4) and th2.coeffs == th2_direct,
           f"exact through q^{PREC - 1}")

    # Level 4: F2 = (4 E2(4τ) - E2)/3, G2 = 4 E2(4τ) - 4 E2(2τ) + E2 (modular.py).
    F2_4 = e2_minus(4).scale(Fraction(-1, 3))
    # G2 = 2 (E2 - 2 E2(2τ)) - (E2 - 4 E2(4τ)): two E2(τ) - N E2(Nτ) forms, no f(Mτ) needed.
    G2_4 = e2_minus(2).scale(2) - e2_minus(4)
    identity("level 4: F2 = θ3^4", F2_4, th3)
    identity("level 4: G2 = θ4^4 - θ2^4", G2_4, th4 - th2)
    identity("Jacobi: θ3^4 = θ2^4 + θ4^4", th3, th2 + th4)

    # Level 2: F2 = 2 E2(2τ) - E2, F4 = (4 E4(2τ) - E4)/3; A = η^16/η(2τ)^8, B = η(2τ)^16/η^8.
    F2_2 = e2_minus(2).scale(-1)
    F4_2 = (e4(2).scale(4) - e4()).scale(Fraction(1, 3))
    A = eta_quotient(2, {1: 16, 2: -8})
    B = eta_quotient(2, {2: 16, 1: -8})
    identity("level 2: F2^2 = A + 64 B", F2_2 * F2_2, A + B.scale(64))
    identity("level 2: F4 = A - 64 B", F4_2, A - B.scale(64))

    # Level 2, the λ form of Theorem 9.3, with τ -> 2τ so λ = θ2^4/θ3^4 at nome q:
    # 16 λ^2 (1-λ) / (2-λ)^4 = 256 t / (1+64t)^2, t = (η(4τ)/η(2τ))^24, cleared of denominators.
    D2 = eta_quotient(4, {2: 24})
    D4 = eta_quotient(4, {4: 24})
    identity("level 2: Theorem 9.3 λ form = eta form (τ -> 2τ)",
             (th2 * th2 * th4 * th3 * (D2 + D4.scale(64)) ** 2).scale(16),
             (D4 * D2 * (th3 + th4) ** 4).scale(256))

    # Level 1: modular.py has J1 = E4^3/(E4^3 - E6^2). The discriminant identity
    # E4^3 - E6^2 = 1728 η^24 turns that into J1 = j/1728 with j = E4^3/η^24.
    identity("level 1: E4^3 - E6^2 = 1728 η^24", e4() ** 3 - e6() ** 2, eta_quotient(1, {1: 24}).scale(1728))

    # Level 1: j = E4^3/η^24 starts q^-1 + 744 + ..., the expansion of Apostol's 1728 J
    # (Thm 1.20, p. 21). With Apostol Thm 2.6 this gives j = 1728 J (docs/proof.md, (3b)).
    prod24 = _pow(_euler(1), 24)  # η^24 / q = prod (1 - q^n)^24
    q_times_j = _mul(_mul(e4().coeffs, e4().coeffs), _mul(e4().coeffs, _inv(prod24)))
    lead = [int(c) for c in q_times_j[:4]]
    record("level 1: j = E4^3/η^24 = q^-1 + 744 + ... (Apostol Thm 1.20)",
           lead[:2] == [1, 744], "first coefficients of q·j: " + ", ".join(map(str, lead)))

    # Level 1: j = 256 (1-λ+λ^2)^3 / (λ^2 (1-λ)^2), with τ -> 2τ: j(2τ) = E4(2τ)^3 / η(2τ)^24.
    P = th3 * th3 - th2 * th3 + th2 * th2
    identity("level 1: j = 256(1-λ+λ²)³/(λ²(1-λ)²) (τ -> 2τ)",
             (P ** 3 * D2).scale(256),
             e4(2) ** 3 * th2 * th2 * th3 * th3 * th4 * th4)

    # Level 3: F2 = (3 E2(3τ) - E2)/2, F4 = (9 E4(3τ) - E4)/8; B3 = η^9/η(3τ)^3, C3 = 27 η(3τ)^9/η^3.
    F2_3 = e2_minus(3).scale(Fraction(-1, 2))
    F4_3 = (e4(3).scale(9) - e4()).scale(Fraction(1, 8))
    B3 = eta_quotient(3, {1: 9, 3: -3})
    C3 = eta_quotient(3, {3: 9, 1: -3}, const=27)
    a = theta_a()
    b_eta = eta_quotient(9, {1: 3, 3: -1})
    record("b = (3 a(q^3) - a(q))/2 as q-series (proved by the lattice argument in docs/proof.md)",
           series_b() == [Fraction(3 * x - y, 2) for x, y in zip(a.at_multiple(3).coeffs, a.coeffs)],
           f"exact through q^{PREC - 1}")
    identity("level 3: a^2 = F2", a * a, F2_3)
    identity("level 3: b = η(τ)^3/η(3τ)", (a.at_multiple(3).scale(3) - a).scale(Fraction(1, 2)), b_eta)
    identity("level 3: F2^3 = (B3 + C3)^2", F2_3 ** 3, (B3 + C3) ** 2)
    identity("level 3: F4^2 = F2 (B3 - C3)^2", F4_3 * F4_3, F2_3 * (B3 - C3) ** 2)

    step3_hauptmoduln()

    negative_controls(A, B, F2_2, th3)

    failed = [name for name, ok, _ in rows if not ok]
    RESULTS.mkdir(exist_ok=True)
    lines = ["Sturm-bound proofs of the modular identities in docs/proof.md, step (2d)", ""]
    for name, ok, detail in rows:
        lines.append(f"{'PASS' if ok else 'FAIL'}  {name}")
        lines.append(f"      {detail}")
    lines.append("")
    lines.append("All identities hold." if not failed else "Failed: " + ", ".join(failed))
    out = RESULTS / "sturm.txt"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print(lines[-1])
    print(f"wrote {out}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
