# Ramanujan's 1/π series and harmonic numbers

A proof, with computer verification, of Conjecture 1 from [MathOverflow question 507693](https://mathoverflow.net/questions/507693/ramanujans-pi-series-and-harmonic-numbers).

## The result

For a rational hypergeometric Ramanujan series $\sum (an+b)f(n)$, the conjecture says that the index derivative of the summand satisfies

$$\mathrm{Re}\,\sum_{n\ge0} f'(n) = -2\pi\mathrm{Im}(\tau)\sum_{n\ge0} f(n),$$

where $\tau$ is the point attached to the series in the table of Cohen and Guillera ([arXiv:2101.12592](https://arxiv.org/abs/2101.12592), Section 3).

The proof here establishes a stronger, complex identity:

$$\sum_{n\ge0} f'(n) = 2\pi i\,\tilde\tau\sum_{n\ge0} f(n),\qquad \tilde\tau=\begin{cases}\tau, & a>0,\\ \tau-1, & a<0,\end{cases}$$

where $f(x) = A(x)\,e^{-x\mathrm{Log}\,a}$ and $a = J_N(\tau)$.

- The identity holds for every $\tau$ on the imaginary axis above $i/\sqrt N$. It also holds for every $\tau$ on the line $\mathrm{Re}\,\tau=\tfrac12$ above the corner of the fundamental domain with $-1\le 1/J_N(\tau)<0$.
- It doesn't use complex multiplication, so $\tau$ need not be a CM point.
- Conjecture 1, for all 36 convergent series in the table, is the special case at the printed CM points.
- Series (2) of the question (Q1) is the level-4 case $\tau=i\sqrt3/2$.

## Contents

| Path | Description |
| --- | --- |
| [`docs/proof.md`](docs/proof.md) | The full proof, with every cited input given by page |
| [`ramanujan_harmonic/`](ramanujan_harmonic) | Python package: the series table, numerical tests and exact verifications |
| [`results/`](results) | Reports produced by the scripts below |

## The package

| Module | Purpose |
| --- | --- |
| `table.py` | The 44 series of the Cohen–Guillera table (36 convergent, 8 divergent), with $N$, $\tau$, $a$ and $k$ |
| `terms.py` | The summands $f(n)$ and their index derivatives $f'(n)$, via the digamma function |
| `modular.py` | The Hauptmoduln $J_N$, computed independently from Eisenstein series |
| `q1_audit.py` | Numerical audit of the existing answer to Q1 |
| `conjecture.py` | Tests the identity on all 36 convergent series |
| `proof_checks.py` | Numerical and symbolic checks of each intermediate claim in the proof |
| `sturm.py` | Exact proofs, by Sturm's bound, of the modular-form identities used in the proof |

## Running

Requires Python 3.11 or later.

```bash
python -m pip install -r requirements.txt

python -m ramanujan_harmonic.q1_audit       # writes results/q1-audit.txt
python -m ramanujan_harmonic.conjecture     # writes results/conjecture.txt
python -m ramanujan_harmonic.proof_checks   # prints to the console
python -m ramanujan_harmonic.sturm          # writes results/sturm.txt
```

Each command exits with status 1 if any check fails.

To test only the three series discussed in the question, plus the complex example from its note, run:

```bash
python -m ramanujan_harmonic.conjecture --ids 35 23 11 1
```

## What is proved and what is cited

The proof is complete modulo published results. Each one is cited by page in [`docs/proof.md`](docs/proof.md):

- **[C1]** Abramowitz and Stegun, *Handbook of Mathematical Functions*, 15.3.10: the hypergeometric connection formula.
- **[C3]** Ramanujan's inversion theorems, in the form of Jacobi (level 4) and Berndt, Bhargava and Garvan, *Trans. AMS* 347 (1995), 4163–4244 (levels 1–3).
- **[C4]** Standard facts on modular curves: Diamond and Shurman, Serre, and Apostol. These are used to show that $1/J_N$ is a Hauptmodul. The fundamental domain itself is proved in the text.
- **[C5]** The pairing $J_N(\tau)=a$ printed in the Cohen–Guillera table. It is used only to identify the table rows.
- **[C6]** Sturm's theorem and the eta-quotient criteria (Ono, *The Web of Modularity*), and modularity of the Eisenstein and theta series (Diamond and Shurman).

Clausen's identity and the symmetric-square property, which are often cited, are proved in step (2b) by an exact symbolic computation. The modular identities linking the inversion theorems to $J_N$ are proved in step (2d) by exact coefficient comparison up to the Sturm bound (`sturm.py`).

## References

- M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions*, NBS, 1964.
- B. C. Berndt, S. Bhargava and F. G. Garvan, Ramanujan's theories of elliptic functions to alternative bases, *Trans. Amer. Math. Soc.* 347 (1995), 4163–4244.
- H. Cohen and J. Guillera, Rational hypergeometric Ramanujan identities for 1/π^c: survey and generalizations, [arXiv:2101.12592](https://arxiv.org/abs/2101.12592).
- F. Diamond and J. Shurman, *A First Course in Modular Forms*, GTM 228, Springer, 2005.
- K. Ono, *The Web of Modularity*, CBMS 102, AMS, 2004.
- J.-P. Serre, *A Course in Arithmetic*, GTM 7, Springer, 1973.
- T. M. Apostol, *Modular Functions and Dirichlet Series in Number Theory*, GTM 41, Springer.

## Citation and licensing

Use the author name **Open** when citing this work. GitHub citation metadata is in [CITATION.cff](CITATION.cff), and the release is archived on Zenodo. This is not a peer-reviewed publication.

Code and software documentation are licensed under [MIT](LICENSE). The proof in `docs/proof.md` and the numerical reports in `results/` are licensed under [CC BY 4.0](LICENSE-CONTENT.md). The series parameters in `table.py` are transcribed from the Cohen–Guillera table and credited to its authors; cited publications are outside these grants.
