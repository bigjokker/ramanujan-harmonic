# Proof of Conjecture 1 (MathOverflow 507693, Q3)

## Statement

Fix a level $N\in\{1,2,3,4\}$ and put

| $N$ | 1 | 2 | 3 | 4 |
| --- | --- | --- | --- | --- |
| $s=s_N$ | $1/6$ | $1/4$ | $1/3$ | $1/2$ |
| $C_N$ | $1728$ | $256$ | $108$ | $64$ |

Let

$$A(x)=\frac{\Gamma(s+x)\,\Gamma(\tfrac12+x)\,\Gamma(1-s+x)}{\Gamma(s)\,\Gamma(\tfrac12)\,\Gamma(1-s)\,\Gamma(1+x)^3}.$$

By Gauss's multiplication formula this is the same analytic function as $H_N(x)$ in `ramanujan_harmonic/terms.py`, for example $A(n)=\binom{2n}{n}^3/64^n$ when $N=4$. For a real $a$ with $|a|\ge1$ put $f(x)=A(x)\,e^{-x\operatorname{Log}a}$, with the principal logarithm. So $f(n)$ is the term of the post and $f'(n)$ is its index derivative. Write $J_N$ for the Hauptmodul of `ramanujan_harmonic/modular.py` and $X_N=1/J_N$, and put $v_N=\sqrt{1/N-1/4}$.

**Theorem.** Let $N\in\{1,2,3,4\}$, and let $\tau$ satisfy one of the following:

- (a) $\tau=iv$ with $v>1/\sqrt N$;
- (b) $\tau=\tfrac12+iv$ with $v>v_N$ and $-1\le X_N(\tau)<0$.

Put $a=J_N(\tau)$. Then $a$ is real, with $a>1$ in case (a) and $a\le-1$ in case (b), and

$$\sum_{n\ge0}f'(n)=2\pi i\,\tilde\tau\sum_{n\ge0}f(n),\qquad
\tilde\tau=\begin{cases}\tau,&a>0,\\ \tau-1,&a<0.\end{cases}$$

Complex multiplication is not used: $\tau$ need not be a CM point.

**Corollary (Conjecture 1).** For every convergent row of the Cohen–Guillera table (rows 1–36), the printed $\tau$ and $a$ satisfy the hypotheses of the Theorem with $a=J_N(\tau)$. So the identity holds for each row, and taking real parts gives $\operatorname{Re}\sum f'(n)=-2\pi\operatorname{Im}(\tau)\sum f(n)$. Step 5 checks the hypotheses. The only appeal to the table is the pairing $J_N(\tau)=a$, which is [C5].

The imaginary part depends on how $a^{-x}$ is extended, but the real part does not. For row 1, writing the term with $(-15)^{-3x}$ and $\operatorname{Log}(-15)$ instead of $(-3375)^{-x}$ and $\operatorname{Log}(-3375)$ changes $f'(n)$ by $-2\pi i\,f(n)$. The Theorem fixes $e^{-x\operatorname{Log}a}$ with $a=J_N(\tau)$, which for the table rows is the printed $a$.

The proof uses the classical inputs marked **[C1]** and **[C3]–[C6]**, listed at the end. The label [C2] (Clausen's identity and the symmetric square) and the fundamental domain in [C4] are proved in the text. Everything else is proved here.

## Notation

Write $z=1/a$, $\theta=z\,d/dz$, and

$$L=\theta^3-z(\theta+s)(\theta+\tfrac12)(\theta+1-s).$$

Let $D'=\{|z|<1\}\setminus(-1,0]$ be the slit disk. On $D'$ define

$$y_0(z)=\sum_{n\ge0}A(n)z^n,\qquad h(z)=\sum_{n\ge0}A'(n)z^n,\qquad Y_1(z)=y_0(z)\operatorname{Log}z+h(z).$$

Since $A(n)\asymp n^{-3/2}$ and $A'(n)/A(n)=\psi(s+n)+\psi(\tfrac12+n)+\psi(1-s+n)-3\psi(1+n)=O(1/n)$, both $y_0$ and $h$ converge absolutely and uniformly on the closed disk $|z|\le1$.

At $z=1/a$:

- If $a>1$, then $\operatorname{Log}(1/a)=-\log a$, so $y_0(1/a)=\sum f(n)$ and $Y_1(1/a)=\sum f'(n)$.
- If $a\le-1$, then $-\operatorname{Log}a=\log|1/a|-i\pi$. That is the limit of $\operatorname{Log}z$ as $z\to1/a$ from $\operatorname{Im}z<0$, so
$$\sum f'(n)=\lim_{z\to1/a,\ \operatorname{Im}z<0}Y_1(z).\tag{0}$$

## Step 1. $Y_1$ solves $L$ (Frobenius)

From $\Gamma(w+1)=w\Gamma(w)$,

$$(n+1+x)^3A(n+1+x)=(n+x+s)(n+x+\tfrac12)(n+x+1-s)\,A(n+x).$$

Put $Y(x,z)=\sum_{n\ge0}A(n+x)z^{n+x}$. Since $\theta z^{w}=wz^{w}$, the recurrence makes $LY$ telescope down to its $n=0$ term:

$$L\,Y(x,z)=x^3A(x)\,z^x .$$

The series converges locally uniformly for $x$ near $0$ and $z$ in compact subsets of $D'$, so we may differentiate in $x$ at $x=0$. The right side has a triple zero at $x=0$, so $L\,\partial_xY|_{x=0}=0$. Since $\partial_xY|_{x=0}=Y_1$, both $y_0$ and $Y_1$ solve $L$ on $D'$.

## Step 2. The mirror map, on the imaginary axis

Let $F(x)={}_2F_1(s,1-s;1;x)$ and $A_2(x)=(s)_x(1-s)_x/\Gamma(1+x)^2$, using the Gamma extension. Let

$$\Phi(x)=\sum_{n\ge0}\bigl(A_2'(n)+A_2(n)\log x\bigr)x^n ,$$

which is the Frobenius log-solution of the second-order equation $E_s$ satisfied by $F$. Take $0<x<\tfrac12$ and $X=4x(1-x)\in(0,1)$.

**(2a)** $F(1-x)=-\dfrac{\sin\pi s}{\pi}\,\Phi(x)$.

This is **[C1]**, Abramowitz–Stegun 15.3.10, p. 559, taking $a=s$, $b=1-s$ (so $c=a+b=1$) and $z=1-x$. The formula holds for $|\arg(1-z)|<\pi$ and $|1-z|<1$, which covers $0<x<\tfrac12$. The bracket $2\psi(n+1)-\psi(s+n)-\psi(1-s+n)-\log x$ there equals $-(A_2'(n)/A_2(n)+\log x)$.

**(2b)** $y_0(X)=F(x)^2$ and $Y_1(X)=F(x)\,\Phi(x)$.

**[C2], proved here.** Let $L^*$ be $L$ pulled back along $X=4x(1-x)$. On $0<x<\tfrac12$ we have $X'(x)=4(1-2x)\ne0$, so $Y(X)$ solves $L$ exactly when $Y(4x(1-x))$ solves $L^*$. Here $E_s$ is $x(1-x)F''+(1-2x)F'-s(1-s)F=0$.

- *Products of solutions of $E_s$ solve $L^*$.* Take any solutions $y,z$ of $E_s$. Substitute $u=yz$ into $L^*$ and eliminate $y'',z''$ (and their derivatives) using $E_s$. The result is identically $0$, as a rational function of $x$ and $s$. `proof_checks.py` does this computation exactly in sympy, with $s$ symbolic. As a control, replacing $s(1-s)$ by $s^2$ gives a nonzero result.
- *They span the solutions.* If $F,\Phi$ are independent solutions of $E_s$, then $F^2,F\Phi,\Phi^2$ are linearly independent. A relation $\alpha F^2+\beta F\Phi+\gamma\Phi^2=0$ would make $\Phi/F$ a root of a fixed quadratic, hence constant. $L^*$ has order $3$, so these three span its solutions.
- *Clausen.* $y_0(X)$ and $F(x)^2$ both solve $L^*$, are analytic at $x=0$, and equal $1$ there. All three indicial roots of $L$ at $X=0$ are $0$, so its analytic solutions are the multiples of $y_0$. Hence $y_0(X)=F(x)^2$, which is Clausen's identity in this case.

So $Y_1(X)=\alpha F^2+\beta F\Phi+\gamma\Phi^2$ for some constants. Now compare behaviour as $x\to0^+$:

- $Y_1$ contains $\log x$ only to the first power, so $\gamma=0$.
- $Y_1(X)=\log X+A'(0)+o(1)=\log x+2\log2+A'(0)+o(1)$.
- $F\Phi=\log x+A_2'(0)+o(1)$.

Hence $\beta=1$ and $\alpha=2\log2+A'(0)-A_2'(0)=2\log2+\psi(\tfrac12)-\psi(1)=0$.

**(2c)** This is **[C3]**, Ramanujan's inversion theorem in signature $r=1/s$. Berndt–Bhargava–Garvan [BBG], *Trans. AMS* 347 (1995), pp. 4165–4166, define the nome $q_r(x)$ in (1.6) by the formula below with $s=1/r$; (1.7)–(1.9) are the cases $r=3,4,6$. For each level, a theorem expresses $x$ through theta functions of the nome, and $4x(1-x)$ then follows by algebra:

- $N=4$, $r=2$: Jacobi. $x=k^2=\lambda=\theta_2^4/\theta_3^4$ at the classical nome $q=e^{-\pi K'/K}$. Abramowitz–Stegun define the nome in §16.27, p. 576. Their 16.38.5 and 16.38.7, p. 579, give $(2K/\pi)^{1/2}=\theta_3(0,q)$ and $(2m^{1/2}K/\pi)^{1/2}=\theta_2(0,q)$ with $m=k^2$, so $\theta_2^4/\theta_3^4=m$.
- $N=3$, $r=3$: [BBG] Lemma 2.6, (2.26), ${}_2F_1(\tfrac13,\tfrac23;1;c^3/a^3)=a(q)$, with the cubic theta functions (2.2)–(2.4) and $a^3=b^3+c^3$ (Theorem 2.2). Lemma 2.9 gives $q_3(c^3/a^3)=q$, and Theorem 2.10, (2.32), packages both. So $x=c^3/a^3$ and $4x(1-x)=4b^3c^3/a^6$.
- $N=2$, $r=4$: [BBG] Theorem 9.3, (9.4), $q_4(x)=q(\lambda)^2$ with $\lambda=2\sqrt x/(1+\sqrt x)$ and $q$ the classical nome. So $x=(\lambda/(2-\lambda))^2$ and $4x(1-x)=16\lambda^2(1-\lambda)/(2-\lambda)^4$. Equation (9.4) is stated for $0<x<1$. On $x\in(0,\tfrac12]$, $\lambda$ runs through $(0,2(\sqrt2-1)]$ and $v$ from $\infty$ down to $1/\sqrt2$, which is the ray.
- $N=1$, $r=6$: [BBG] Theorem 11.3, (11.11), $q_6(\beta)=q(\alpha)^2$ with $\alpha,\beta$ as in (11.1). This gives $4\beta(1-\beta)=\tfrac{27}{4}\alpha^2(1-\alpha)^2/(1-\alpha+\alpha^2)^3$. The hypothesis is $0<p<1$, and the proof of Theorem 11.1 shows $\alpha$ and $\beta$ increasing from $0$ to $1$. At $\alpha=\tfrac12$ the right side is $1$, so $\beta=\tfrac12$ there as well. The ray $v\ge1$ is therefore the part $\alpha\le\tfrac12$, $\beta\le\tfrac12$. On it, $\beta$ is the smaller root of the quadratic, which is the root `proof_checks.py` uses. By the $j$–$\lambda$ formula proved in (2d), this is $1728/j=1/J_1$.

**(2d)** The expressions in (2c) equal $1/J_N$, where $J_N$ is defined by the Eisenstein formulas of `modular.py`. Those formulas give $1/J_1=1728/j$, $1/J_N=1-F_4^2/F_2^4$ for $N=2,3$, and $1/J_4=1-G_2^2/F_2^2$. Each step below is an identity between holomorphic modular forms of one weight and character on $\Gamma_0(N)$. Each is proved by **[C6]**: Sturm's theorem plus a finite exact computation in `ramanujan_harmonic/sturm.py`, which compares coefficients through the Sturm bound $\lfloor k\,[\mathrm{SL}_2(\mathbb Z):\Gamma_0(N)]/12\rfloor$. Every eta quotient used is placed in its space by **[C6]**: the Gordon–Hughes–Newman congruences give weight and character on $\Gamma_0(N)$, and Ligozat's orders at the cusps are nonnegative. `results/sturm.txt` lists, for each quotient, the level, weight, character, the two congruence sums and the order at every cusp. At nome $q=e^{2\pi i\tau}$, the theta sums are eta quotients by [Ono] Thm 1.60, p. 17: $\theta_3=\sum q^{n^2}=\eta(2\tau)^5/(\eta(\tau)^2\eta(4\tau)^2)$ and $\theta_4=\sum(-1)^nq^{n^2}=\eta(\tau)^2/\eta(2\tau)$. For $\theta_2=\sum q^{(n+1/2)^2}$, Ono gives $\eta(16z)^2/\eta(8z)=\sum_{n\ge0}q^{(2n+1)^2}$, and $z\mapsto z/4$ gives $\theta_2=2\eta(4\tau)^2/\eta(2\tau)$. Hence $\theta_3^4=\eta(2\tau)^{20}/(\eta(\tau)^8\eta(4\tau)^8)$, $\theta_4^4=\eta(\tau)^8/\eta(2\tau)^4$ and $\theta_2^4=16\eta(4\tau)^8/\eta(2\tau)^4$, all in $M_2(\Gamma_0(4))$ with trivial character. The script also compares these quotients with the theta sums, as a check of the transcription.

- $N=4$. $F_2=\theta_3^4$, $G_2=\theta_4^4-\theta_2^4$ and $\theta_3^4=\theta_2^4+\theta_4^4$, all in $M_2(\Gamma_0(4))$. Hence $1/J_4=\bigl(\theta_3^8-(\theta_4^4-\theta_2^4)^2\bigr)/\theta_3^8=4\theta_2^4\theta_4^4/\theta_3^8=4\lambda(1-\lambda)$.
- $N=2$. Put $A=\eta(\tau)^{16}/\eta(2\tau)^8$ and $B=\eta(2\tau)^{16}/\eta(\tau)^8$. Then $F_2^2=A+64B$ and $F_4=A-64B$ in $M_4(\Gamma_0(2))$. So $1/J_2=256AB/(A+64B)^2=256t/(1+64t)^2$ with $t=(\eta(2\tau)/\eta(\tau))^{24}$. To meet Theorem 9.3, where $\lambda$ has the classical nome $e^{i\pi\tau}$, replace $\tau$ by $2\tau$. After that substitution, $\lambda=\theta_2^4/\theta_3^4$ and $t=(\eta(4\tau)/\eta(2\tau))^{24}$, with all theta functions at nome $q$. Using $1-\lambda=\theta_4^4/\theta_3^4$ and $2-\lambda=(\theta_3^4+\theta_4^4)/\theta_3^4$, the claim $16\lambda^2(1-\lambda)/(2-\lambda)^4=256t/(1+64t)^2$ is equivalent to

$$16\,\theta_2^8\theta_4^4\theta_3^4\,\bigl(\eta(2\tau)^{24}+64\,\eta(4\tau)^{24}\bigr)^2=256\,\eta(2\tau)^{24}\eta(4\tau)^{24}\,\bigl(\theta_3^4+\theta_4^4\bigr)^4 .$$

$t$ itself has a pole at a cusp, which is why the denominators are cleared first. Each factor is a holomorphic form on $\Gamma_0(4)$ with trivial character: $\theta_i^4\in M_2$, and $\eta(2\tau)^{24},\eta(4\tau)^{24}\in M_{12}$, with Ligozat orders $2,2,2$ and $4,1,1$ at $\infty,\tfrac12,0$. So both sides lie in $M_{32}(\Gamma_0(4))$, and the Sturm bound is $32\cdot6/12=16$. An identity that holds at $2\tau$ for every $\tau$ holds at $\tau$.
- $N=1$. `modular.py` defines $J_1=E_4^3/(E_4^3-E_6^2)$. The discriminant identity $E_4^3-E_6^2=1728\,\eta(\tau)^{24}$ holds in $M_{12}(\mathrm{SL}_2(\mathbb Z))$: $E_6\in M_6(\mathrm{SL}_2(\mathbb Z))$ by [DS] §1.1, and $\eta^{24}$ passes the eta-quotient checks at level 1. The Sturm bound is $q^1$. So $1/J_1=1728/j$ with $j=E_4^3/\eta^{24}$. Then, after $\tau\mapsto2\tau$, $j=256(1-\lambda+\lambda^2)^3/(\lambda^2(1-\lambda)^2)$ clears to $256\,(\theta_3^8-\theta_2^4\theta_3^4+\theta_2^8)^3\,\eta(2\tau)^{24}=E_4(2\tau)^3\,\theta_2^8\theta_3^8\theta_4^8$ in $M_{24}(\Gamma_0(4))$.
- $N=3$. Put $B_3=\eta(\tau)^9/\eta(3\tau)^3$ and $C_3=27\eta(3\tau)^9/\eta(\tau)^3$. The following hold:
  - $a^2=F_2$ in $M_2(\Gamma_0(3))$. Here $a\in M_1(\Gamma_0(3),(\tfrac{-3}{\cdot}))$ by [C6], so $a^2$ has trivial character.
  - $b=(3a(q^3)-a(q))/2$. This is elementary. The terms with $m\equiv n\pmod 3$ sum to $a(q^3)$: put $m=n+3k$, so $m^2+mn+n^2=3(n^2+3kn+3k^2)$, and $n=n'-k$ turns $n^2+3kn+3k^2$ into $n'^2+n'k+k^2$. The classes $m-n\equiv1,2$ are swapped by $(m,n)\mapsto(n,m)$, and $\omega+\omega^2=-1$.
  - $b=\eta(\tau)^3/\eta(3\tau)$ in $M_1(\Gamma_0(9),(\tfrac{-3}{\cdot}))$, so $b^3=B_3$.
  - $F_2^3=(B_3+C_3)^2$ in $M_6(\Gamma_0(3))$, so $a^6=(B_3+C_3)^2$. Since $a^3$ and $B_3+C_3$ are holomorphic with constant term $1$, $a^3=B_3+C_3$.

  Theorem 2.2 then gives $c^3=a^3-b^3=C_3$, so no product formula for $c$ is needed. So $4x(1-x)=4B_3C_3/(B_3+C_3)^2$. Finally, $F_4^2=F_2(B_3-C_3)^2$ in $M_8(\Gamma_0(3))$, so $1/J_3=1-(B_3-C_3)^2/(B_3+C_3)^2=4B_3C_3/(B_3+C_3)^2$.

Put $q=\exp\bigl(-\pi\csc(\pi s)\,F(1-x)/F(x)\bigr)$ and $q=e^{2\pi i\tau}$, so that $\tau=iv$. Then

$$X_N(\tau)=4x(1-x),$$

where $X_N=1/J_N$ is the Hauptmodul computed in `modular.py`, with $X_N=C_Nq+O(q^2)$. As $x$ runs over $(0,\tfrac12]$, $v$ runs over $[1/\sqrt N,\infty)$. At $x=\tfrac12$ we get $q=e^{-\pi\csc\pi s}=e^{-2\pi/\sqrt N}$.

**Consequence (M).** Combining 2a–2d: for every $v>1/\sqrt N$, $X_N(iv)\in(0,1)$ and

$$\frac{Y_1(X_N(iv))}{y_0(X_N(iv))}=\frac{\Phi(x)}{F(x)}=-\pi\csc(\pi s)\frac{F(1-x)}{F(x)}=\log q=2\pi i\,(iv).$$

As a consistency check of the constants, $C_N=e^{-A'(0)}$. For example, for $N=1$, $A'(0)=-\log1728$.

## Step 3. The geometry of $X_N$

Let $D_N=\{\tau:\ |\operatorname{Re}\tau|<\tfrac12,\ |\tau|>1/\sqrt N\}$. Let $D_N^\pm$ be its parts with $\operatorname{Re}\tau\gtrless0$, and put $v_N=\sqrt{1/N-1/4}$, so $v_4=0$. The corner of $D_N$ is $\tfrac12+iv_N$: it is $\rho+1=e^{\pi i/3}$ for $N=1$ (with $\rho=e^{2\pi i/3}$), a point of $\mathbb H$ for $N=2,3$, and the cusp $\tfrac12$ for $N=4$.

The group is $\Gamma_0(N)^+=\Gamma_0(N)\cup W_N\Gamma_0(N)$ with $W_N=\tfrac1{\sqrt N}\left[\begin{smallmatrix}0&-1\\N&0\end{smallmatrix}\right]$, and $\Gamma_0(1)^+=\mathrm{SL}_2(\mathbb Z)$. Every element has determinant $1$, so $\operatorname{Im}g\tau=\operatorname{Im}\tau/|\gamma\tau+\delta|^2$ for its bottom row $(\gamma,\delta)$. For $g\in\Gamma_0(N)$ the bottom row is $(c,d)$ with $N\mid c$ and $\gcd(c,d)=1$. The coset $W_N\Gamma_0(N)$ consists of the matrices $\tfrac1{\sqrt N}\left[\begin{smallmatrix}Nx&y\\Nc&Nd\end{smallmatrix}\right]$ with $Nxd-yc=1$, so its bottom rows are $\sqrt N(c,d)$ with $\gcd(c,Nd)=1$. For $N>1$ this forces $c\neq0$. For $N=1$ the coset is $\mathrm{SL}_2(\mathbb Z)$ itself, so it needs no separate treatment below.

**(3a) The domain.** This replaces a citation; for $N=1$ it is Serre, *A Course in Arithmetic*, Ch. VII, Thm 1, pp. 77–78.

*(i) Every $\tau\in\mathbb H$ is equivalent to a point of $\overline{D_N}$.* The bottom rows form a discrete subset of $\mathbb R^2$ and $\tau\notin\mathbb R$, so $|\gamma\tau+\delta|$ has a positive minimum over nonzero rows and $\operatorname{Im}g\tau$ attains a maximum. Take $g$ with $\operatorname{Im}g\tau$ maximal and translate by a power of $T$ so that $|\operatorname{Re}g\tau|\le\tfrac12$. If $|g\tau|<1/\sqrt N$, then $\operatorname{Im}W_Ng\tau=\operatorname{Im}g\tau/(N|g\tau|^2)>\operatorname{Im}g\tau$, a contradiction. So $g\tau\in\overline{D_N}$.

*(ii) If $\tau\in D_N$ and $g\in\Gamma_0(N)^+$ is not $\pm T^k$, then $\operatorname{Im}g\tau<\operatorname{Im}\tau$.* For $\tau\in D_N$ and integers $c\neq0$, $d$,
$$|c\tau+d|^2=c^2|\tau|^2+2cd\operatorname{Re}\tau+d^2>c^2/N-|cd|+d^2,$$
because $|\tau|^2>1/N$ and $|\operatorname{Re}\tau|<\tfrac12$. It remains to show the right side is at least $1$, or at least $1/N$ for the coset, where $|\gamma\tau+\delta|^2=N|c\tau+d|^2$.

- On $\Gamma_0(N)$, write $c=Nc'$ with $c'\ne0$. The bound $c^2/N-|cd|+d^2$ is $c^2-|cd|+d^2$ for $N=1$, $c'^2+(|c'|-|d|)^2$ for $N=2$, $(|d|-\tfrac32|c'|)^2+\tfrac34c'^2$ for $N=3$, and $(2|c'|-|d|)^2$ with $d$ odd for $N=4$. Each is a positive integer.
- On the coset, $N(c^2/N-|cd|+d^2)=c^2-N|cd|+Nd^2$. This is $(|c|-|d|)^2+d^2$ for $N=2$, $(|c|-\tfrac32|d|)^2+\tfrac34d^2$ (equal to $c^2$ if $d=0$) for $N=3$, and $(|c|-2|d|)^2$ with $c$ odd for $N=4$. Again each is a positive integer.

So $|\gamma\tau+\delta|^2>1$. Two consequences follow:

- No non-translation maps a point of $D_N$ into $D_N$; a translation cannot either, because $D_N$ has width $1$. So distinct points of $D_N$ are inequivalent, and no point of $D_N$ is an elliptic point.
- No point $\tau\in D_N$ is equivalent to a point $\beta\in\partial D_N\cap\mathbb H$. If $g\tau=\beta$, then $g$ maps a neighbourhood of $\tau$ inside $D_N$ onto a neighbourhood of $\beta$, which meets $D_N$. That contradicts the first consequence.

**(3b) The Hauptmodul.** Put $j=E_4^3/\eta^{24}$ for $N=1$, $t=(\eta(N\tau)/\eta(\tau))^{24/(N-1)}$ for $N=2,3$, and $c_N=N^{-12/(N-1)}$, that is $c_2=1/4096$ and $c_3=1/729$. For $N=4$ put $\lambda=\theta_2^4/\theta_3^4=16\eta(\tau)^8\eta(4\tau)^{16}/\eta(2\tau)^{24}$. By (2d), $X_1=1728/j$, $X_2=256t/(1+64t)^2$, $X_3=108t/(1+27t)^2$ and $X_4=4\lambda(1-\lambda)$.

- *Fricke.* By Apostol, *Modular Functions and Dirichlet Series*, Thm 3.1, p. 48, $\eta(-1/\tau)=(-i\tau)^{1/2}\eta(\tau)$, and every power that occurs is an integer, so the branch does not matter. This gives $t(W_N\tau)=c_N/t(\tau)$. For $N=4$ it gives $\lambda(W_4\tau)=\eta(\tau)^{16}\eta(4\tau)^8/\eta(2\tau)^{24}=\theta_4^4/\theta_3^4=1-\lambda(\tau)$, by [Ono] Thm 1.60 and Jacobi's identity from (2d).
- *Degree one on $X_0(N)$.* Let $f$ denote $j$, $t$ or $\lambda$ according to $N$. For $N=2,3,4$, $f$ is an eta quotient of weight $0$ and trivial character that passes [Ono] Thm 1.64. For $N=1$, $j$ is $E_4^3$ (in $M_{12}(\mathrm{SL}_2(\mathbb Z))$, [DS] §1.1) divided by $\eta^{24}$, which never vanishes on $\mathbb H$. So $f$ lies in $A_0(\Gamma_0(N))$, the field of meromorphic functions on the compact Riemann surface $X_0(N)$ ([DS] pp. 62, 72). None of them has a pole in $\mathbb H$, because $\eta$ never vanishes there. By [DS] §3.8, p. 103, $\Gamma_0(N)$ has $\varphi(\gcd(d,N/d))$ cusps for each $d\mid N$. Each count is $1$ here, giving the cusp $\infty$ for $N=1$, $\infty,0$ for $N=2,3$, and $\infty,0,\tfrac12$ for $N=4$; `sturm.py` checks the count. The orders at the cusps are as follows. $j=q^{-1}+744+\cdots$ has a simple pole at $\infty$ (`sturm.py`). By [Ono] Thm 1.65, $t$ has order $1$ at $\infty$ and $-1$ at $0$, and $\lambda$ has orders $1,0,-1$ at $\infty,0,\tfrac12$; `sturm.py` prints these. So in every case $f$ has exactly one pole on $X_0(N)$, and it is simple. For $w\in\mathbb C$, $\operatorname{div}(f-w)$ has degree $0$ ([DS] §3.4, p. 83), so $f-w$ has exactly one zero, counted with multiplicity. Hence $j$, $t$ and $\lambda$ are bijections $X_0(N)\to\mathbb P^1$. Since $\Gamma_0(4)$ has no elliptic points ([DS] Cor. 3.7.2, p. 96), the local coordinate of $X_0(4)$ at every point of $\mathbb H$ is $\tau$ itself, so $\lambda'\neq0$ on $\mathbb H$.
- *Injective on $\Gamma_0(N)^+\backslash\mathbb H$.* The differences factor as follows:
  - $X_2(t)-X_2(t')$ is a nonzero multiple of $(t-t')(1-4096\,tt')$;
  - $X_3(t)-X_3(t')$ is a nonzero multiple of $(t-t')(1-729\,tt')$;
  - $X_4(\lambda)-X_4(\lambda')=4(\lambda-\lambda')(1-\lambda-\lambda')$.

  So $X_N(\tau_1)=X_N(\tau_2)$ forces $f(\tau_2)=f(\tau_1)$ or $f(\tau_2)=f(W_N\tau_1)$, and hence $\tau_2\in\Gamma_0(N)\tau_1\cup\Gamma_0(N)W_N\tau_1$. For $N=1$, $W_1=S\in\mathrm{SL}_2(\mathbb Z)$, and $X_1=1728/j$ is injective because $j$ is.
- *Image.* $X_N(\mathbb H)\supseteq\mathbb C\setminus\{0\}$. For $N=2,3$, $t(\mathbb H)=\mathbb C\setminus\{0\}$, and for $w\neq0$ the quadratic $X_N(t)=w$ has roots with nonzero product $c_N$. For $N=4$, $\lambda(\mathbb H)=\mathbb C\setminus\{0,1\}$, and $\lambda(1-\lambda)=w/4\neq0$ avoids $\lambda=0,1$. For $N=1$, $j(\mathbb H)=\mathbb C$, since the cusp goes to $\infty$, so $X_1(\mathbb H)=\mathbb P^1\setminus\{0\}$. $X_N$ has rational $q$-coefficients, so $X_N(-\bar\tau)=\overline{X_N(\tau)}$.
- *Values.* At $i/\sqrt N$, which $W_N$ fixes, $t^2=c_N$ with $t>0$ (since $\eta(iy)>0$), so $t=N^{-6/(N-1)}$ and $X_N=1$. For $N=4$, $\lambda=1-\lambda$ gives $\lambda=\tfrac12$ and $X_4=1$. For $N=1$, $E_6(-1/\tau)=\tau^6E_6(\tau)$ at $\tau=i$ gives $E_6(i)=0$. Then $j(i)=1728$, because $j-1728=E_6^2/\eta^{24}$ by the discriminant identity of (2d), and so $X_1(i)=1$. At the corner $\tau_c$ for $N=2,3$, we have $|\tau_c|^2=1/N$, so $W_N\tau_c=-\bar\tau_c=\tau_c-1$ and $t(\tau_c)^2=c_N$. On $\operatorname{Re}\tau=\tfrac12$, $q<0$ and $t/q>0$, so $t(\tau_c)=-N^{-6/(N-1)}$, that is $-1/64$ or $-1/27$, which is exactly the pole of $X_N$. For $N=1$ the corner is $\rho+1$ with $\rho=e^{2\pi i/3}$. $E_4(-1/\tau)=\tau^4E_4(\tau)$ at $\tau=\rho$, together with $-1/\rho=\rho+1$, gives $E_4(\rho)=\rho E_4(\rho)$. So $E_4(\rho)=0$, and $X_1=1728\eta^{24}/E_4^3$ has its pole at the corner. By injectivity, $X_N=\infty$ only on the orbit of the corner (for $N=4$, nowhere in $\mathbb H$), and by (3a) that orbit misses $D_N$.
- *Derivative.* $X_N'(\sigma)\neq0$ at every $\sigma\in\mathbb H$ with $X_N(\sigma)\notin\{1,\infty\}$.
  - *Elliptic points of $\Gamma_0(N)$.* [DS] Cor. 3.7.2, p. 96, gives one orbit of period $2$ and one of period $3$ for $N=1$; one of period $2$ for $N=2$; one of period $3$ for $N=3$; and none for $N=4$. For $N=1$ these are the orbits of $i$ and $\rho$, where $X_1=1$ and $\infty$. For $N=2,3$ the corner is fixed by $\left[\begin{smallmatrix}1&-1\\2&-1\end{smallmatrix}\right]\in\Gamma_0(2)$, whose square is $-I$, so it acts on $\mathbb H$ with order $2$; and by $\left[\begin{smallmatrix}1&-1\\3&-2\end{smallmatrix}\right]\in\Gamma_0(3)$, whose cube is $I$, so it acts with order $3$. So the corner is the elliptic orbit, where $X_N=\infty$.
  - *At other points.* At a point of $\mathbb H$ that is not elliptic, the local coordinate of $X_0(N)$ is $\tau$ itself ([DS] p. 62). The degree-one function $f$ therefore has $f'\neq0$ there.
  - *The rational function.* $dX_N/df$ is a nonzero multiple of $1-64t$, $1-27t$ or $1-2\lambda$ for $N=2,3,4$, and equals $-1728/j^2$ for $N=1$. For $N=2,3,4$ it vanishes only where $X_N=1$ ($t=1/64$, $t=1/27$, $\lambda=\tfrac12$); for $N=1$ it never vanishes on $\mathbb H$. These derivatives are a direct computation: $dX_2/dt=256(1-64t)/(1+64t)^3$, $dX_3/dt=108(1-27t)/(1+27t)^3$, $dX_4/d\lambda=4(1-2\lambda)$.

**Lemma.** $X_N$ maps $D_N$ biholomorphically onto $P=\mathbb C\setminus\bigl((-\infty,0]\cup[1,\infty)\bigr)$. It maps $D_N^+$ onto the upper half-plane, $D_N^-$ onto the lower half-plane, and the ray $\{iv:v>1/\sqrt N\}$ onto $(0,1)$.

*Proof.* $X_N$ is holomorphic on $D_N$ because the corner orbit misses it. It is injective on $D_N$ by (3a) and (3b): equal values force equivalent points, and distinct points of $D_N$ are inequivalent.

$X_N$ is real on the imaginary axis, on $\operatorname{Re}\tau=\pm\tfrac12$ (where $X_N(-\bar\tau)=X_N(\tau-1)$), and on the arc $|\tau|=1/\sqrt N$ (where $W_N\tau=-\bar\tau$).

If $\tau\in D_N^+$ had $X_N(\tau)$ real, then $-\bar\tau\in D_N^-$ would be a different point of $D_N$ with the same value, contradicting injectivity. So $X_N(D_N^+)$ is a connected subset of $\mathbb C\setminus\mathbb R$. Near the cusp $X_N\approx C_Nq$, which has positive imaginary part for small $\operatorname{Re}\tau>0$, so $X_N(D_N^+)$ lies in the upper half-plane. By symmetry $X_N(D_N^-)$ lies in the lower half-plane.

Conversely, take $w$ with $\operatorname{Im}w>0$. By (3b) $w=X_N(\tau)$ for some $\tau\in\mathbb H$, and by (3a)(i) we may take $\tau\in\overline{D_N}$. $\tau$ is not on $\partial D_N$ or on the axis, where $X_N$ is real, so $\tau\in D_N^+$.

On the ray, $X_N$ is real, continuous and injective, so it is monotone. It tends to $1$ at $i/\sqrt N$ and to $0$ at $i\infty$, so its image is $(0,1)$. Hence $X_N(D_N)$ is the upper half-plane, the lower half-plane and $(0,1)$ together, which is $P$. $\square$

## Step 4. Continuation

Let $\Omega=(X_N|_{D_N})^{-1}(D')$. Since $D'\subset P$, the Lemma shows $\Omega$ is biholomorphic to $D'$. In particular it is connected, and it contains the whole ray $\{iv: v>1/\sqrt N\}$. The function

$$G(\tau)=Y_1(X_N(\tau))-2\pi i\,\tau\;y_0(X_N(\tau))$$

is holomorphic on $\Omega$, and by (M) it vanishes on the ray. By the identity theorem $G\equiv0$ on $\Omega$.

Now fix $\tau$ as in the Theorem, and put $z_0=X_N(\tau)=1/a$.

**Case (a).** $\tau$ is on the ray, so $z_0\in(0,1)$ by the Lemma, $a>1$, and $\tau\in\Omega$. Hence

$$\textstyle\sum f'(n)=Y_1(z_0)=2\pi i\tau\,y_0(z_0)=2\pi i\tau\sum f(n).$$

**Case (b), $-1<z_0<0$.** Here $a<-1$. Put $\sigma=\tau-1$, which lies on the left side of $D_N$, and note $X_N(\sigma)=z_0$ by periodicity. Since $|\sigma|^2=\tfrac14+v^2>\tfrac14+v_N^2=1/N$, the points $\tau_k=\sigma+1/k$ lie in $D_N^-$ for large $k$. Then $z_k=X_N(\tau_k)\to z_0$ with $\operatorname{Im}z_k<0$ by the Lemma, and $|z_k|<1$ for large $k$, so $\tau_k\in\Omega$. Using $G(\tau_k)=0$, letting $k\to\infty$, and applying (0):

$$\textstyle\sum f'(n)=\lim Y_1(z_k)=\lim 2\pi i\tau_k\,y_0(z_k)=2\pi i(\tau-1)\sum f(n).$$

**Case (b), $z_0=-1$.** Here $a=-1$ and $1/a$ is on the unit circle. By (3b), $X_N'(\sigma)\neq0$, since $X_N(\sigma)=-1\notin\{1,\infty\}$. So $X_N$ maps a small disk $U$ around $\sigma$ biholomorphically onto a neighbourhood of $-1$. On $U$, points with $\operatorname{Re}\tau>-\tfrac12$ lie in $D_N^-$ and go to the lower half-plane. Points with $\operatorname{Re}\tau<-\tfrac12$ go to the upper half-plane, by periodicity and the Lemma. Choose $w_k\to-1$ with $|w_k|<1$ and $\operatorname{Im}w_k<0$, and let $\tau_k\in U$ be the point with $X_N(\tau_k)=w_k$. Then $\tau_k\in D_N^-\cap\Omega$ and $\tau_k\to\sigma$, and the limit is as before, using the continuity of $y_0$ and $h$ on the closed disk from the Notation section. $\square$

## Step 5. The table rows

The printed $\tau$ satisfy the position hypotheses:

- *Rows 8–11, 18–23, 30–32, 35, 36* have $\tau=iv^*$ with $v^*>1/\sqrt N$. The smallest cases are $\sqrt2>1$, $1>1/\sqrt2$, $\sqrt6/3>1/\sqrt3$ and $\sqrt3/2>\tfrac12$.
- *Rows 1–7, 12–17, 24–29, 33, 34* have $\tau=\tfrac12+iv^*$ with $v^*>v_N$. The smallest cases are $\sqrt7/2>\sqrt3/2$, $\sqrt5/2>\tfrac12$, $\sqrt{27}/6>1/(2\sqrt3)$ and $\sqrt2/2>0$.

By **[C5]**, $J_N(\tau)=a$ for the printed $a$. For the rows of the second kind, $a\le-1$, so $X_N(\tau)=1/a\in[-1,0)$; row 33 is the case $a=-1$. So the Theorem applies to every row, with $f$ built from the printed $a$, and this proves the Corollary. $\square$

## What is cited and what is checked

- **[C1]** Abramowitz–Stegun, *Handbook of Mathematical Functions* (NBS 1964), 15.3.10, p. 559, read on the page.
- **[C2]** Not cited: proved in (2b). The symmetric-square property and Clausen's identity are proved in (2b) by an exact symbolic computation (`proof_checks.py`) and elementary linear ODE facts.
- **[C3]** Ramanujan's inversion theorem: Jacobi for $N=4$ (Abramowitz–Stegun §16.27, p. 576, and 16.38.5, 16.38.7, p. 579), and [BBG] Theorem 2.10 with Lemmas 2.6 and 2.9 for $N=3$, Theorem 9.3 for $N=2$, Theorem 11.3 for $N=1$. Read on the page in the journal (*Trans. AMS* 347 (1995), 4163–4244): (1.6) p. 4165 and (1.7)–(1.9) pp. 4165–4166, Lemma 2.6 p. 4171, Lemma 2.9 and Theorem 2.10 p. 4172, Theorem 9.3 p. 4214, (11.1) and Theorem 11.1 p. 4228 (the monotonicity of $\alpha,\beta$ is at the end of its proof, p. 4230), Theorem 11.3 p. 4230. The author-hosted PDF, https://qseries.org/fgarvan/papers/alternativebases.pdf, has the same numbering on different pages: 3, 9, 10, 54, 69 and 70. The links from these statements to $J_N$ are proved in (2d).
- **[C4]** The fundamental domain of $\Gamma_0(N)^+$ and the Hauptmodul property of $X_N$, used in Step 3. The domain is proved in (3a) by an elementary inequality. The Hauptmodul property is proved in (3b) from these sources, each read on the page cited:
  - Serre, *A Course in Arithmetic*, Ch. VII, Thm 1, pp. 77–78: the classical case $N=1$ of (3a), which (3a) also proves directly.
  - Apostol, *Modular Functions and Dirichlet Series*, Thm 3.1, p. 48 (the transformation of $\eta$). As a remark, Apostol's Thm 1.20, p. 21, gives $1728J=q^{-1}+744+\cdots$ for his modular invariant $J$, and `sturm.py` gives the same two leading terms for $j=E_4^3/\eta^{24}$. By his Thm 2.6, p. 39, $j-1728J$ is then a bounded modular function vanishing at $\infty$, so $j=1728J$ and $X_1=1/J$. The proof does not need this identification.
  - [DS] pp. 62 and 72 ($X_0(N)$ is a compact Riemann surface whose meromorphic functions are $A_0(\Gamma_0(N))$), §3.4, p. 83 (divisors of meromorphic functions have degree $0$), §3.8, p. 103 (cusps of $\Gamma_0(N)$), and Cor. 3.7.2, p. 96 ($\Gamma_0(4)$ has no elliptic points).
  - [Ono] Thms 1.60, 1.64 and 1.65.
- **[C5]** The pairing $J_N(\tau)=a$ for rows 1–36, printed in Cohen–Guillera, arXiv:2101.12592v1, Section 3, p. 7, read on the page. Used only in Step 5. The same reading shows that rows 37–44, on p. 8, match `table.py` in $N$, $\tau$, $a$, $P$, and $k$. The 45-digit agreement in `modular.py` checks that transcription against the Eisenstein series; it is not the citation.
- **[C6]** Standard facts on modular forms, used in (2d). [Ono] is *The Web of Modularity*, CBMS 102, AMS 2004. [DS] is Diamond–Shurman, *A First Course in Modular Forms*, GTM 228. Each item below was read on the page cited.
  - Sturm's theorem, [Ono] Thm 2.58, §2.9, p. 40 (Sturm 1987). Ono states it modulo an ideal $\mathfrak m$, for weight written $k/2$. For integer weight $w$ the condition reads $\operatorname{ord}_{\mathfrak m}(f)>w\,[\Gamma_0(1):\Gamma_0(N)]/12$. To get equality, clear denominators in the difference of the two sides and apply the theorem for every prime. Remark 2.59's extra condition $4\mid N$ is for half-integral weight only, so it does not affect the weight-1 identity on $\Gamma_0(9)$. The bound is applied only to holomorphic forms of equal weight, level and character, comparing coefficients $n\le w[\mathrm{SL}_2(\mathbb Z):\Gamma_0(N)]/12$, that index included.
  - The Gordon–Hughes–Newman criterion, [Ono] Thm 1.64, p. 18. If $\sum r_\delta$ is even and $\sum\delta r_\delta\equiv\sum(N/\delta)r_\delta\equiv0\pmod{24}$, the quotient is weakly holomorphic of weight $k=\tfrac12\sum r_\delta$ and character $\bigl((-1)^k\prod\delta^{r_\delta}/d\bigr)$ on $\Gamma_0(N)$.
  - Ligozat's cusp orders, [Ono] Thm 1.65, p. 18: the order at $c/d$ is $\tfrac N{24}\sum_\delta\gcd(d,\delta)^2r_\delta/(\gcd(d,N/d)\,d\,\delta)$. Holomorphy requires all of these to be $\ge0$.
  - `sturm.py` checks the hypotheses of Thms 1.64 and 1.65 for each eta quotient it uses.
  - The eta-quotient forms of $\theta_2,\theta_3,\theta_4$, [Ono] Thm 1.60, p. 17.
  - $E_4\in M_4(\mathrm{SL}_2(\mathbb Z))$: [DS] §1.1, pp. 5–6, with $E_k=G_k/(2\zeta(k))$. Also $E_2(\tau)-NE_2(N\tau)\in M_2(\Gamma_0(N))$: [DS] §1.2, p. 18, states $G_{2,N}=G_2(\tau)-NG_2(N\tau)\in M_2(\Gamma_0(N))$ (proof left as Exercise 1.2.8(e)), and $G_2=2\zeta(2)E_2$. The level-4 form $G_2=2\bigl(E_2-2E_2(2\tau)\bigr)-\bigl(E_2-4E_2(4\tau)\bigr)$ is a combination of two such forms, using $\Gamma_0(4)\subset\Gamma_0(2)$.
  - $f(M\tau)\in M_k(\Gamma_0(NM),\chi)$ for $f\in M_k(\Gamma_0(N),\chi)$. This is used for $E_4(2\tau)$, $E_4(3\tau)$ and $a(3\tau)$. [DS] Exercise 1.2.11(d), p. 24, states it without the character, so here is the short proof. Let $\delta=\left[\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right]\in\Gamma_0(NM)$, so $NM\mid c$, and put $\delta'=\left[\begin{smallmatrix}a&Mb\\c/M&d\end{smallmatrix}\right]$. Then $\det\delta'=1$, $N\mid c/M$, so $\delta'\in\Gamma_0(N)$, and $M\cdot\delta\tau=\delta'(M\tau)$. Since $ad\equiv1\pmod{NM}$, $d$ is coprime to $NM$, so $\chi(d)$ is defined and $\chi$ serves as a character modulo $NM$. Since $(c/M)(M\tau)+d=c\tau+d$, this gives $f(M\delta\tau)=\chi(d)(c\tau+d)^kf(M\tau)$. For the cusps, take $\gamma\in\mathrm{SL}_2(\mathbb Z)$ and write $\left[\begin{smallmatrix}M&0\\0&1\end{smallmatrix}\right]\gamma=\gamma'\left[\begin{smallmatrix}\alpha&\beta\\0&\delta_0\end{smallmatrix}\right]$ with $\gamma'\in\mathrm{SL}_2(\mathbb Z)$ and $\alpha\delta_0=M$ (Hermite normal form). With the weight-$k$ operator extended to $\mathrm{GL}_2^+(\mathbb Q)$ as in [DS] Exercise 1.2.11, $f(M\tau)$ is a nonzero constant times $f|_k\left[\begin{smallmatrix}M&0\\0&1\end{smallmatrix}\right]$. So $f(M\tau)|_k\gamma$ is a constant times $(f|_k\gamma')|_k\left[\begin{smallmatrix}\alpha&\beta\\0&\delta_0\end{smallmatrix}\right]$, which is $(f|_k\gamma')\bigl((\alpha\tau+\beta)/\delta_0\bigr)$ up to a constant. Since $f$ is holomorphic at the cusps, $f|_k\gamma'$ has an expansion in nonnegative powers of $e^{2\pi i\tau/h}$ for some $h$. Substituting $(\alpha\tau+\beta)/\delta_0$ keeps the powers nonnegative, so $f(M\tau)$ is holomorphic at every cusp.
  - The theta series $a(q)$ of $m^2+mn+n^2$ lies in $M_1(\Gamma_0(3),(\tfrac{-3}{\cdot}))$: [DS] Thm 4.11.3, p. 159, with $N=1$ and $\chi$ trivial. That theorem puts $\tfrac16\sum_{n\in\mathbb Z[\mu_3]}e^{2\pi i|n|^2\tau}$ in $M_1(\Gamma_0(3),\psi)$ with $\psi(d)=(d/3)=(\tfrac{-3}{d})$. The norm there is $|x|^2=x_1^2-x_1x_2+x_2^2$, which becomes $m^2+mn+n^2$ under $x_2\mapsto-n$, so this sum is $a(q)/6$.

`python -m ramanujan_harmonic.proof_checks` confirms the intermediate claims numerically. It covers (2a) and (2b). For (2c) and (2d) it checks the identification $4x(1-x)=1/J_N(\tau)$ with $J_N$ from `modular.py`, the [BBG] statements as transcribed above (Theorem 2.2 with the eta products, Lemma 2.9, Theorems 9.3 and 11.3), and the algebra from each to $1/J_N$, which it verifies symbolically where it is a rational identity. It also covers the values of $X_N$ at $i/\sqrt N$ and at the corner, and the monotonicity on the ray and the sides. It exits 1 if any check fails. The identities in (2d) are proved, not just checked, by `python -m ramanujan_harmonic.sturm`, which writes `results/sturm.txt`.

For $a<0$, the informal phrase "continuing from the cusp along $z\in[0,1/a]$" should be read as *approaching $1/a$ from $\operatorname{Im}z<0$*. The segment itself lies on the branch cut of $\operatorname{Log}$, and approaching from below is what selects $\tau-1$.
