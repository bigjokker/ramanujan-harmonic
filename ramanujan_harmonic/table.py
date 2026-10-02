"""The 36 convergent and 8 divergent rational hypergeometric 1/π series.

Source: Henri Cohen and Jesús Guillera, "Rational Hypergeometric Ramanujan
Identities for 1/π^c: Survey and Generalizations", arXiv:2101.12592v1,
Section 3. Rows 1–36 are the table on p. 7 (the table the MathOverflow
post calls "page 7"). Rows 37–44 are the divergent table on p. 8. The
sign column on p. 8 is not stored: it is not the sign of k (row 40 has k < 0 and sign +).
For row 37 the paper prints k = 3; the analytic continuation of its series
gives sqrt(3)/π, which confirms the printed value.

Each series is

    sum_{n>=0} (slope * n + intercept) * H_N(n) / a^n = sqrt(k) / π,

with H_N the hypergeometric product for level N:

    N = 1:  R_n(2) R_n(6)
    N = 2:  R_n(2) R_n(4)
    N = 3:  R_n(2) R_n(3)
    N = 4:  R_n(2)^3

and τ = (tau_p + sqrt(-tau_D)) / tau_q a CM point with J_N(τ) = a.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


def rat(spec: str) -> Fraction:
    """Parse a signed prime factorization such as '-2^9*3^-2*5^3'."""
    text = spec.replace(" ", "")
    sign = 1
    if text.startswith("-"):
        sign = -1
        text = text[1:]
    elif text.startswith("+"):
        text = text[1:]
    value = Fraction(1)
    for part in text.split("*"):
        if part == "":
            raise ValueError(spec)
        if "^" in part:
            base, exp = part.split("^")
            value *= Fraction(int(base)) ** int(exp)
        else:
            value *= int(part)
    return sign * value


@dataclass(frozen=True)
class Series:
    id: int
    level: int
    tau_p: int
    tau_q: int
    tau_D: int
    a: Fraction
    slope: int
    intercept: int
    k: Fraction
    convergent: bool
    label: str = ""

    def tau_text(self) -> str:
        disc = f"sqrt(-{self.tau_D})" if self.tau_D != 1 else "sqrt(-1)"
        if self.tau_p == 0 and self.tau_q == 1:
            return disc
        if self.tau_p == 0:
            return f"{disc}/{self.tau_q}"
        if self.tau_q == 1:
            return f"({self.tau_p}+{disc})"
        return f"({self.tau_p}+{disc})/{self.tau_q}"

    def im_tau_squared(self) -> Fraction:
        """Im(τ)^2 = tau_D / tau_q^2. Im(τ) itself is the positive square root."""
        return Fraction(self.tau_D, self.tau_q * self.tau_q)


def _s(
    id_: int,
    level: int,
    tau_p: int,
    tau_q: int,
    tau_D: int,
    a: str,
    slope: int,
    intercept: int,
    k: str,
    convergent: bool = True,
    label: str = "",
) -> Series:
    return Series(
        id_,
        level,
        tau_p,
        tau_q,
        tau_D,
        rat(a),
        slope,
        intercept,
        rat(k),
        convergent,
        label,
    )


# Numbers 1–36 are the convergent table. Numbers 37–44 are the divergent table
# (|a| < 1); those sums are not classical and are excluded from Conjecture 1.
ENTRIES: tuple[Series, ...] = (
    _s(1, 1, 1, 2, 7, "-2^-6*5^3", 63, 8, "3*5^3", label="complex example in the note"),
    _s(2, 1, 1, 2, 11, "-2^9*3^-3", 154, 15, "2^11"),
    _s(3, 1, 1, 2, 19, "-2^9", 342, 25, "2^11*3"),
    _s(4, 1, 1, 2, 27, "-2^9*3^-2*5^3", 506, 31, "2^11*3^-3*5^3"),
    _s(5, 1, 1, 2, 43, "-2^12*5^3", 5418, 263, "2^14*3^-1*5^3"),
    _s(6, 1, 1, 2, 67, "-2^9*5^3*11^3", 261702, 10177, "2^11*3*5^3*11^3"),
    _s(7, 1, 1, 2, 163, "-2^12*5^3*23^3*29^3", 545140134, 13591409, "2^14*3*5^3*23^3*29^3"),
    _s(8, 1, 0, 1, 2, "3^-3*5^3", 28, 3, "5^3"),
    _s(9, 1, 0, 1, 3, "2^-2*5^3", 11, 1, "2^-2*3^-1*5^3"),
    _s(10, 1, 0, 1, 4, "2^-3*11^3", 63, 5, "2^-4*3*11^3"),
    _s(11, 1, 0, 1, 7, "2^-6*5^3*17^3", 133, 8, "2^-2*3^-5*5^3*17^3", label="MO series 3"),
    _s(12, 2, 1, 2, 5, "-2^2", 20, 3, "2^6"),
    _s(13, 2, 1, 2, 7, "-2^-8*3^4*7^2", 65, 8, "3^4*7"),
    _s(14, 2, 1, 2, 9, "-2^4*3", 28, 3, "2^8*3^-1"),
    _s(15, 2, 1, 2, 13, "-2^2*3^4", 260, 23, "2^6*3^4"),
    _s(16, 2, 1, 2, 25, "-2^6*3^4*5", 644, 41, "2^10*3^4*5^-1"),
    _s(17, 2, 1, 2, 37, "-2^2*3^4*7^4", 21460, 1123, "2^6*3^4*7^4"),
    _s(18, 2, 0, 1, 1, "2^-5*3^4", 7, 1, "2^-2*3^4"),
    _s(19, 2, 0, 2, 6, "3^2", 8, 1, "2^2*3"),
    _s(20, 2, 0, 2, 10, "3^4", 10, 1, "2^-3*3^4"),
    _s(21, 2, 0, 2, 18, "7^4", 40, 3, "3^-3*7^4"),
    _s(22, 2, 0, 2, 22, "3^4*11^2", 280, 19, "2^2*3^4*11"),
    _s(23, 2, 0, 2, 58, "3^8*11^4", 26390, 1103, "2^-3*3^8*11^4", label="MO series 2"),
    _s(24, 3, 3, 6, 27, "-2^4*3^-2", 5, 1, "2^4*3^-1"),
    _s(25, 3, 3, 6, 51, "-2^4", 51, 7, "2^4*3^3"),
    _s(26, 3, 3, 6, 75, "-2^4*5", 9, 1, "2^4*3*5^-1"),
    _s(27, 3, 3, 6, 123, "-2^10", 615, 53, "2^10*3^3"),
    _s(28, 3, 3, 6, 147, "-2^4*3^3*7", 165, 13, "2^4*3^6*7^-1"),
    _s(29, 3, 3, 6, 267, "-2^4*5^6", 14151, 827, "2^4*3^3*5^6"),
    # a = J_3(sqrt(-6)/3) = 2, not 4. The ar5iv table prints this cell as "2".
    _s(30, 3, 0, 3, 6, "2", 6, 1, "3^3"),
    _s(31, 3, 0, 3, 12, "2^-1*3^3", 15, 2, "2^-4*3^6"),
    _s(32, 3, 0, 3, 15, "2^-2*5^3", 33, 4, "2^-2*3^3*5^2"),
    _s(33, 4, 1, 2, 2, "-1", 4, 1, "2^2"),
    _s(34, 4, 1, 2, 4, "-2^3", 6, 1, "2^3"),
    _s(35, 4, 0, 2, 3, "2^2", 6, 1, "2^4", label="MO series 1, Q1"),
    _s(36, 4, 0, 2, 7, "2^6", 42, 5, "2^8"),
    _s(37, 2, -1, 2, 3, "-2^-4*3^2", 5, 1, "3", convergent=False),  # p. 8 prints k = 3
    _s(38, 2, 1, 4, 7, "2^-8*3^4", 35, 8, "-2^2*3^4", convergent=False),
    _s(39, 3, 2, 6, 2, "2*3^-3", 10, 3, "-2^2*5^2", convergent=False),
    _s(40, 3, 1, 6, 11, "2^4*3^-3", 11, 3, "-2^4*3^2", convergent=False),
    _s(41, 3, 3, 6, 15, "-2^-2", 15, 4, "3^3", convergent=False),
    _s(42, 4, 1, 2, 1, "-2^-3", 3, 1, "1", convergent=False),
    _s(43, 4, 3, 8, 7, "2^-6", 21, 8, "-2^4", convergent=False),
    _s(44, 4, 1, 4, 3, "2^-2", 3, 1, "-2^2", convergent=False),
)

BY_ID = {entry.id: entry for entry in ENTRIES}


def _check_table() -> None:
    if len(ENTRIES) != 44:
        raise RuntimeError(f"expected 44 series, found {len(ENTRIES)}")
    if sum(entry.convergent for entry in ENTRIES) != 36:
        raise RuntimeError("expected 36 convergent series")
    ids = [entry.id for entry in ENTRIES]
    if ids != list(range(1, 45)):
        raise RuntimeError("series ids are not 1..44 in order")
    # The three worked examples in the MathOverflow post, plus the complex one.
    s35, s23, s11, s1 = BY_ID[35], BY_ID[23], BY_ID[11], BY_ID[1]
    if s35.a != 4 or s35.im_tau_squared() != Fraction(3, 4):
        raise RuntimeError("series 35 does not match MO series 1")
    if s23.a != Fraction(3**8 * 11**4) or s23.im_tau_squared() != Fraction(58, 4):
        raise RuntimeError("series 23 does not match MO series 2")
    if s11.a != Fraction(5**3 * 17**3, 2**6) or s11.im_tau_squared() != 7:
        raise RuntimeError("series 11 does not match MO series 3")
    if s1.a != Fraction(-125, 64) or s1.im_tau_squared() != Fraction(7, 4):
        raise RuntimeError("series 1 does not match the complex example")
    if BY_ID[30].a != 2:
        raise RuntimeError("series 30 should have a = 2")


_check_table()
