"""Numerical checks of the intermediate claims in docs/proof.md.

- Step 2a: the connection formula [C1].
- Step 2b: Clausen [C2], and Y1/y0 = Phi/F.
- Step 2c: the nome [C3]. q_r(x) from Berndt–Bhargava–Garvan (1.6), read as
  q_r = e^{2πiτ}, gives 4x(1-x) = 1/J_N(τ) with J_N from modular.py.
- Step 2c sources: the inversion statements cited for [C3], as transcribed from
  Berndt–Bhargava–Garvan, and the algebra that turns each into 1/J_N:
  Lemma 2.9 / Theorem 2.10 (signature 3), Theorem 9.3 (signature 4) and
  Theorem 11.3 (signature 6), checked on the proof's ray v > 1/sqrt(N).
- Step 3: X_N = 1/J_N is 1 at i/sqrt(N) and infinite at the corner. It is real
  and decreasing on the imaginary axis inside (0,1), and real, negative and
  increasing on the vertical side.

Run from the repository root:

    python -m ramanujan_harmonic.proof_checks

Exits 1 if any check fails.
"""

from __future__ import annotations

import mpmath as mp

import sympy as sp

from ramanujan_harmonic.console import configure_stdout
from ramanujan_harmonic.modular import hauptmodul_at

I = mp.mpc(0, 1)
DENOM = {1: 6, 2: 4, 3: 3, 4: 2}  # s = 1/DENOM[N]


def levels():
    """(N, s) pairs, with s computed at the current working precision."""
    return [(N, mp.mpf(1) / d) for N, d in DENOM.items()]


failures: list[str] = []


def check(name: str, ok: bool, detail: str) -> None:
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


def A2(s, n):  # (s)_n (1-s)_n / n!^2, Gamma extension
    return mp.rf(s, n) * mp.rf(mp.mpf(1) - s, n) / mp.gamma(n + 1) ** 2


def A3(s, n):  # (s)_n (1/2)_n (1-s)_n / n!^3
    return mp.rf(s, n) * mp.rf(mp.mpf(1) / 2, n) * mp.rf(mp.mpf(1) - s, n) / mp.gamma(n + 1) ** 3


def frob(coef, z, M=400):  # (y0, Y1) with Y1 = d/dx sum coef(n+x) z^(n+x) at x=0
    y0 = mp.fsum(coef(n) * z**n for n in range(M))
    Y1 = mp.fsum((mp.diff(coef, n) + coef(n) * mp.log(z)) * z**n for n in range(M))
    return y0, Y1


def _pulled_back_L_on_products(coeff) -> "sp.Expr":
    """L, pulled back along X = 4x(1-x), applied to y*z, where y and z solve
    x(1-x)F'' + (1-2x)F' - coeff*F = 0. Returns the result after eliminating y'', z''."""
    x, s = sp.symbols("x s")
    y, z = sp.Function("y"), sp.Function("z")
    X = 4 * x * (1 - x)
    theta = lambda f: X / sp.diff(X, x) * sp.diff(f, x)  # X d/dX in the variable x
    u = y(x) * z(x)
    a = theta(u) + s * u
    b = theta(a) + sp.Rational(1, 2) * a
    c = theta(b) + (1 - s) * b
    expr = theta(theta(theta(u))) - X * c
    p = (1 - 2 * x) / (x * (1 - x))
    q = -coeff(s) / (x * (1 - x))
    for f in (y, z):
        d1, d2, d3 = (sp.Derivative(f(x), (x, k)) for k in (1, 2, 3))
        second = -p * d1 - q * f(x)
        third = sp.diff(second, x).subs(d2, second)
        expr = expr.subs(d3, third).subs(d2, second)
    return sp.simplify(sp.together(sp.expand(expr)))


def step2b_symmetric_square() -> None:
    """docs/proof.md (2b): products of solutions of E_s solve L pulled back by X = 4x(1-x)."""
    res = _pulled_back_L_on_products(lambda s: s * (1 - s))
    check("2b products of E_s solutions solve L(4x(1-x))", res == 0,
          f"exact, s symbolic: L[y z] reduces to {res}")
    control = _pulled_back_L_on_products(lambda s: s**2)
    check("2b control: a wrong E_s fails", control != 0, "with s^2 in place of s(1-s) the result is nonzero")


def step2ab() -> None:
    mp.mp.dps = 30
    tol = mp.mpf(10) ** -27
    for N, s in levels():
        x = mp.mpf("0.2")
        F, Phi = frob(lambda n: A2(s, n), x)
        conn = abs(mp.hyp2f1(s, 1 - s, 1, 1 - x) + mp.sin(mp.pi * s) / mp.pi * Phi)
        y0, Y1 = frob(lambda n: A3(s, n), 4 * x * (1 - x))
        clausen = abs(F**2 - y0)
        ratio = abs(Y1 / y0 - Phi / F)
        check(f"2a connection N={N}", conn < tol, f"|2F1(s,1-s;1;1-x) + sin(πs)/π Φ| = {mp.nstr(conn, 3)}")
        check(f"2b Clausen N={N}", clausen < tol, f"|F^2 - y0(4x(1-x))| = {mp.nstr(clausen, 3)}")
        check(f"2b Y1/y0 = Φ/F N={N}", ratio < tol, f"|Y1/y0 - Φ/F| = {mp.nstr(ratio, 3)}")


def step2c() -> None:
    mp.mp.dps = 40
    for N, s in levels():
        worst = mp.mpf(0)
        for x in ("0.5", "0.35", "0.2", "0.1", "0.02"):
            x = mp.mpf(x)
            # BBG (1.6): q_r = exp(-π csc(π/r) F(1-x)/F(x)), and q_r = e^{2πiτ} with τ = iv.
            v = mp.csc(mp.pi * s) / 2 * mp.hyp2f1(s, 1 - s, 1, 1 - x) / mp.hyp2f1(s, 1 - s, 1, x)
            with mp.workdps(70):  # E4^3 - E6^2 cancels about log10(1/|q|) digits
                inv_j = 1 / hauptmodul_at(N, I * v)
            worst = max(worst, abs(inv_j - 4 * x * (1 - x)) / (4 * x * (1 - x)))
        check(f"2c nome N={N}", worst < mp.mpf(10) ** -35,
              f"worst relative |1/J_N(τ) - 4x(1-x)| over 5 values of x = {mp.nstr(worst, 3)}")


def _nome(r: int, x):
    """Berndt–Bhargava–Garvan (1.6)."""
    s = mp.mpf(1) / r
    F = lambda z: mp.hyp2f1(s, 1 - s, 1, z)
    return mp.exp(-mp.pi * mp.csc(mp.pi * s) * F(1 - x) / F(x))


def _eta(t):
    return mp.exp(2 * mp.pi * I * t / 24) * mp.qp(mp.exp(2 * mp.pi * I * t))


def step2c_sources() -> None:
    mp.mp.dps = 40
    tol = mp.mpf(10) ** -35

    # The algebra, exactly.
    lam, al = sp.symbols("lambda alpha", positive=True)
    x4 = (lam / (2 - lam)) ** 2  # inverse of lambda = 2 sqrt(x)/(1 + sqrt(x))
    sig4 = sp.simplify(4 * x4 * (1 - x4) - 16 * lam**2 * (1 - lam) / (2 - lam) ** 4) == 0
    j = 256 * (1 - al + al**2) ** 3 / (al**2 * (1 - al) ** 2)
    sig6 = sp.simplify(sp.Rational(27, 4) * al**2 * (1 - al) ** 2 / (1 - al + al**2) ** 3 - 1728 / j) == 0
    check("2c algebra, signature 4", sig4, "x = (λ/(2-λ))^2 gives 4x(1-x) = 16λ^2(1-λ)/(2-λ)^4")
    check("2c algebra, signature 6", sig6, "27/4 α^2(1-α)^2/(1-α+α^2)^3 = 1728/j(α)")

    # Signature 4, Theorem 9.3: q4(x) = q(λ)^2 with λ = 2 sqrt(x)/(1+sqrt(x)), q the classical nome.
    thm, link = mp.mpf(0), mp.mpf(0)
    for x in ("0.02", "0.04", "0.1", "0.2", "0.36", "0.45"):
        x = mp.mpf(x)
        lam_x = 2 * mp.sqrt(x) / (1 + mp.sqrt(x))
        q_classical = mp.exp(-mp.pi * mp.ellipk(1 - lam_x) / mp.ellipk(lam_x))
        thm = max(thm, abs(_nome(4, x) / q_classical**2 - 1))
        tau = mp.log(_nome(4, x)) / (2 * mp.pi * I)
        with mp.workdps(70):
            link = max(link, abs(16 * lam_x**2 * (1 - lam_x) / (2 - lam_x) ** 4 * hauptmodul_at(2, tau) - 1))
    check("2c Theorem 9.3", thm < tol, f"worst |q4(x)/q(λ)^2 - 1| = {mp.nstr(thm, 3)}")
    check("2c level 2 link", link < tol, f"worst |16λ^2(1-λ)/(2-λ)^4 · J_2(τ) - 1| = {mp.nstr(link, 3)}")

    # Signature 3: cubic theta series (2.2)-(2.4), eta products, Lemma 2.9 q3(c^3/a^3) = q.
    worst = {"eta": mp.mpf(0), "thm22": mp.mpf(0), "lem29": mp.mpf(0), "alg": mp.mpf(0), "link": mp.mpf(0)}
    for v in ("0.6", "0.8", "1.0", "1.5", "2.5"):  # all above 1/sqrt(3)
        tau = I * mp.mpf(v)
        q = mp.exp(2 * mp.pi * I * tau)
        third = mp.mpf(1) / 3
        rng = range(-40, 40)
        a = mp.fsum(q ** (m * m + m * n + n * n) for m in rng for n in rng)
        c = mp.fsum(q ** ((m + third) ** 2 + (m + third) * (n + third) + (n + third) ** 2) for m in rng for n in rng)
        b_eta, c_eta = _eta(tau) ** 3 / _eta(3 * tau), 3 * _eta(3 * tau) ** 3 / _eta(tau)
        worst["eta"] = max(worst["eta"], abs(c / c_eta - 1))
        worst["thm22"] = max(worst["thm22"], abs((b_eta**3 + c**3) / a**3 - 1))
        x = mp.re(c**3 / a**3)
        worst["lem29"] = max(worst["lem29"], abs(_nome(3, x) / q - 1))
        t = (_eta(3 * tau) / _eta(tau)) ** 12
        worst["alg"] = max(worst["alg"], abs(4 * x * (1 - x) / (108 * t / (1 + 27 * t) ** 2) - 1))
        with mp.workdps(70):
            worst["link"] = max(worst["link"], abs(108 * t / (1 + 27 * t) ** 2 * hauptmodul_at(3, tau) - 1))
    check("2c c = 3η(3τ)^3/η(τ)", worst["eta"] < tol, f"worst relative error {mp.nstr(worst['eta'], 3)}")
    check("2c Theorem 2.2 with b = η(τ)^3/η(3τ)", worst["thm22"] < tol, f"worst |(b^3+c^3)/a^3 - 1| = {mp.nstr(worst['thm22'], 3)}")
    check("2c Lemma 2.9", worst["lem29"] < mp.mpf(10) ** -30, f"worst |q3(c^3/a^3)/q - 1| = {mp.nstr(worst['lem29'], 3)}")
    check("2c level 3 algebra", worst["alg"] < tol, f"worst |4x(1-x) / (108t/(1+27t)^2) - 1| = {mp.nstr(worst['alg'], 3)}")
    check("2c level 3 link", worst["link"] < tol, f"worst |108t/(1+27t)^2 · J_3(τ) - 1| = {mp.nstr(worst['link'], 3)}")

    # Signature 6, Theorem 11.3: q6(β) = q(α)^2, α = λ(τ), on the ray v > 1 only.
    # Below v = 1 the root β <= 1/2 belongs to -1/τ instead.
    thm, link = mp.mpf(0), mp.mpf(0)
    for v in ("1.05", "1.5", "2.5", "4.0"):
        # Near v = 1, Y -> 1 and sqrt(1 - Y) cancels digits, so work at 70 places.
        with mp.workdps(70):
            tau = I * mp.mpf(v)
            qc = mp.exp(I * mp.pi * tau)
            alpha = mp.re((mp.jtheta(2, 0, qc) / mp.jtheta(3, 0, qc)) ** 4)
            Y = mp.mpf(27) / 4 * alpha**2 * (1 - alpha) ** 2 / (1 - alpha + alpha**2) ** 3
            beta = (1 - mp.sqrt(1 - Y)) / 2
            thm = max(thm, abs(_nome(6, beta) / qc**2 - 1))
            link = max(link, abs(Y * hauptmodul_at(1, tau) - 1))
    check("2c Theorem 11.3", thm < tol, f"worst |q6(β)/q(α)^2 - 1| = {mp.nstr(thm, 3)}")
    check("2c level 1 link", link < tol, f"worst |4β(1-β) · J_1(τ) - 1| = {mp.nstr(link, 3)}")


def step3() -> None:
    mp.mp.dps = 30

    def X(N, t):
        return 1 / hauptmodul_at(N, t)

    def real(u):
        return abs(mp.im(u)) < mp.mpf(10) ** -20 * max(1, abs(u))

    for N in (1, 2, 3, 4):
        fix = X(N, I / mp.sqrt(N))
        vc = mp.sqrt(mp.mpf(1) / N - mp.mpf(1) / 4)
        corner = X(N, mp.mpf(1) / 2 + I * max(vc, mp.mpf("0.02")))  # N=4: corner is the cusp 1/2
        with mp.workdps(60):
            if N == 4:  # the cusp 1/2 needs many terms near |q| = 1
                corner = 1 / hauptmodul_at(N, mp.mpf(1) / 2 + I * mp.mpf("0.02"), nmax=4000)
        vs = [mp.mpf(1) / mp.sqrt(N) * (1 + mp.mpf(k) / 20) for k in range(1, 60)]
        axis = [X(N, I * v) for v in vs]
        side = [X(N, mp.mpf(1) / 2 + I * (vc + mp.mpf(k) / 20)) for k in range(1, 60)]
        mono_axis = all(real(u) and 0 < mp.re(u) < 1 for u in axis) and all(
            mp.re(axis[i + 1]) < mp.re(axis[i]) for i in range(len(axis) - 1)
        )
        mono_side = all(real(u) and mp.re(u) < 0 for u in side) and all(
            mp.re(side[i + 1]) > mp.re(side[i]) for i in range(len(side) - 1)
        )
        check(f"3 X(i/√N) = 1, N={N}", abs(fix - 1) < mp.mpf(10) ** -25, f"X = {mp.nstr(mp.re(fix), 15)}")
        check(f"3 corner is a pole, N={N}", abs(corner) > mp.mpf(10) ** 20, f"|X(corner)| = {mp.nstr(abs(corner), 5)}")
        check(f"3 axis in (0,1), decreasing, N={N}", mono_axis, "59 sample points")
        check(f"3 side negative, increasing, N={N}", mono_side, "59 sample points")


def step3b() -> None:
    """docs/proof.md (3b): factorizations, Fricke relations, and values at the corner."""
    x, y = sp.symbols("x y")
    facts = (
        ("X2", 256 * x / (1 + 64 * x) ** 2 - 256 * y / (1 + 64 * y) ** 2, (x - y) * (1 - 4096 * x * y)),
        ("X3", 108 * x / (1 + 27 * x) ** 2 - 108 * y / (1 + 27 * y) ** 2, (x - y) * (1 - 729 * x * y)),
        ("X4", 4 * x * (1 - x) - 4 * y * (1 - y), 4 * (x - y) * (1 - x - y)),
    )
    for name, diff, factor in facts:
        num = sp.factor(sp.numer(sp.together(diff)))
        quotient = sp.simplify(num / factor)  # must be a nonzero constant
        ok = quotient.is_number and quotient != 0
        check(f"3b factorization {name}", bool(ok), f"numerator of the difference = {num}")

    mp.mp.dps = 30
    eta = lambda tau: mp.exp(2 * mp.pi * I * tau / 24) * mp.qp(mp.exp(2 * mp.pi * I * tau))
    tau0 = mp.mpc("0.137", "0.811")
    for N, r in ((2, 24), (3, 12)):
        tfun = lambda tau: (eta(N * tau) / eta(tau)) ** r
        W = lambda tau: -1 / (N * tau)
        c = mp.mpf(N) ** (-mp.mpf(r) / 2)
        err = abs(tfun(W(tau0)) * tfun(tau0) / c - 1)
        corner = mp.mpf(1) / 2 + I * mp.sqrt(mp.mpf(1) / N - mp.mpf(1) / 4)
        tc = tfun(corner)
        target = -mp.mpf(N) ** (-mp.mpf(6) / (N - 1))
        check(f"3b Fricke t(W τ) = c/t(τ), N={N}", err < mp.mpf(10) ** -25, f"relative error {mp.nstr(err, 3)}")
        check(f"3b corner value N={N}", abs(tc - target) < mp.mpf(10) ** -25, f"t(corner) = {mp.nstr(mp.re(tc), 15)}, expected {mp.nstr(target, 15)}")
    lam = lambda tau: 16 * eta(tau) ** 8 * eta(4 * tau) ** 16 / eta(2 * tau) ** 24
    err = abs(lam(-1 / (4 * tau0)) - (1 - lam(tau0)))
    check("3b Fricke λ(W4 τ) = 1 - λ(τ)", err < mp.mpf(10) ** -25, f"error {mp.nstr(err, 3)}")


def main() -> int:
    configure_stdout()
    step2b_symmetric_square()
    step2ab()
    step2c()
    step2c_sources()
    step3()
    step3b()
    print()
    print("All checks pass." if not failures else "Failed: " + ", ".join(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
