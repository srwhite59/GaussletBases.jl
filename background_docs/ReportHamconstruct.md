# Radial boundary–gausslet basis with $Y_{\ell m}$: operators, angular algebra, and a **two-index** integral–diagonal plan

This note assumes you already have:

* A **boundary gausslet–like** radial basis ${\chi_a(r)}*{a=1}^{N_r}$ on $r\in[0,\infty)$ that is orthonormal in the radial metric,
  $$\int_0^\infty \chi_a(r),\chi_b(r),dr ;=; \delta*{ab}.$$
* A high-quality **radial grid** ${(r_i,w_i)}_{i=1}^{N_g}$ that integrates smooth functions on $[0,\infty)$ accurately.
* Samples $\chi_a(r_i)$ and derivatives $\chi_a'(r_i)$ available with good accuracy.

All one- and two-electron quantities below are built **numerically** on this grid. No oddization tricks and no analytic radial integrals are required.

---

## 1) Radial samplers and one-electron operators

Define the grid samplers
$$
\Phi_{ia} \equiv \chi_a(r_i), \qquad D_{ia} \equiv \chi_a'(r_i),
$$
and the weight–scaled versions
$$
\widetilde\Phi_{ia}\equiv \sqrt{w_i},\Phi_{ia}, \qquad \widetilde D_{ia}\equiv \sqrt{w_i},D_{ia}.
$$

Because the continuum basis is orthonormal, the discrete overlap is the identity to high accuracy:
$$
S_{ab} ;\equiv; \sum_i w_i,\chi_a(r_i)\chi_b(r_i)
;=; (\widetilde\Phi^\top \widetilde\Phi)*{ab};\approx;\delta*{ab}.
$$

### Kinetic energy (radial core, $\ell$-independent)

$$
T^{(0)}*{ab} ;=; -\tfrac12\int_0^\infty \chi_a'(r)\chi_b'(r),dr
;\approx; \tfrac12,(\widetilde D^\top \widetilde D)*{ab}.
$$

### Nuclear attraction ($\ell=0$)

$$
V^{\text{nuc}}_{ab} ;=; \int_0^\infty \chi_a(r),\Big(-\tfrac{Z}{r}\Big),\chi_b(r),dr
;\approx; \widetilde\Phi^\top !\Big[\mathrm{diag}!\big(-Z/r_i\big)\Big]\widetilde\Phi.
$$

### Centrifugal term (general $\ell$)

For an $(\ell,m)$ channel,
$$
V^{\text{cent}}_{ab}(\ell) ;=; \int_0^\infty \chi_a(r),\frac{\ell(\ell+1)}{2,r^2},\chi_b(r),dr
;\approx; \widetilde\Phi^\top !\Big[\mathrm{diag}!\big(\tfrac{\ell(\ell+1)}{2,r_i^2}\big)\Big]\widetilde\Phi.
$$

Hence the one-electron block for fixed $(\ell,m)$ is
$$
H^{(1)}(\ell) ;=; T^{(0)} + V^{\text{nuc}} + V^{\text{cent}}(\ell),
$$
and it is **block-diagonal in $(\ell,m)$**.

---

## 2) Extend the basis with $Y_{\ell m}$

Use the product basis
$$
\langle{a;\ell m} ;\equiv; \frac{\chi_a(r)}{r},Y_{\ell m}(\Omega),
$$
so a one-electron orbital reads
$$
\psi(\mathbf r)=\sum_{a,\ell,m} C_{a\ell m},\frac{\chi_a(r)}{r},Y_{\ell m}(\Omega).
$$

The angular part is orthonormal on the sphere, so the global overlap is $\delta_{ab}\delta_{\ell\ell'}\delta_{mm'}$, and one-electron operators remain diagonal in $(\ell,m)$.

---

## 3) **Two-index** integral–diagonal approximation (IDA) for Coulomb

For $s$-shell ($L=0$), the radial Coulomb kernel is
$$
K^{(0)}(r,r')=\frac{1}{\max(r,r')}.
$$

Let the **integrated weights of basis functions** be
$$
w_a^\chi ;\equiv; \int_0^\infty \chi_a(r),dr
;\approx; \sum_{i=1}^{N_g} w_i,\chi_a(r_i).
$$

Form
$$
M ;\equiv; \Phi^\top \mathrm{diag}(w) \in \mathbb R^{N_r\times N_g},
\qquad
K^{(0)}_{ij} ;\equiv; w_i w_j ,\frac{1}{\max(r_i,r_j)}.
$$

Then the **two-index** IDA Coulomb matrix **in the contracted basis** is
$$
\boxed{ ;
V^{(ee)} ;\approx; D^{-1},\Big( M,K^{(0)},M^\top \Big),D^{-1},
\qquad
D \equiv \mathrm{diag}(w^\chi_1,\dots,w^\chi_{N_r}).
;}
$$

This is exactly the “integral diagonal” rule you’ve been using:
$$
V^{(ee)}_{ab} ;\approx; \frac{ \displaystyle \iint \chi_a(r),\frac{dr,dr'}{\max(r,r')},\chi_b(r') }{ w_a^\chi,w_b^\chi }.
$$

### Higher multipoles for general $Y_{\ell m}$ couplings

With the spherical expansion,
$$
\frac{1}{|\mathbf r-\mathbf r'|}=\sum_{L=0}^\infty \frac{4\pi}{2L+1},\frac{r_<^L}{r_>^{L+1}}
\sum_{M=-L}^L Y_{LM}(\Omega)Y_{LM}^*(\Omega'),
$$
define on the grid
$$
K^{(L)}*{ij} ;\equiv; w_i w_j,\frac{r*<^L}{r_>^{L+1}}\Big|_{(r_i,r_j)}.
$$
The **radial** IDA block for that $L$ is the same two-index contraction:
$$
V^{(ee,L)} ;\approx; D^{-1},\big(M,K^{(L)},M^\top\big),D^{-1}.
$$

Angular coupling among $(\ell,m)$ channels is handled separately via Gaunt coefficients (next section). The radial part feeding those couplings is precisely $V^{(ee,L)}$.

---

## 4) Angular algebra (Gaunt tensors)

Define the Gaunt coefficients
$$
G^{L}*{\ell_1 m_1,\ell_2 m_2;M}
=\int Y*{\ell_1 m_1}(\Omega),Y_{LM}(\Omega),Y_{\ell_2 m_2}^*(\Omega),d\Omega
=(-1)^{m_2}\sqrt{\frac{(2\ell_1+1)(2L+1)}{4\pi(2\ell_2+1)}}
\begin{pmatrix}\ell_1 & L & \ell_2\ 0 & 0 & 0\end{pmatrix}
\begin{pmatrix}\ell_1 & L & \ell_2\ m_1 & M & -m_2\end{pmatrix}.
$$

The four-index Coulomb elements in the product basis factor as
$$
\begin{aligned}
\langle a\ell_1 m_1,, b\ell_2 m_2 | r_{12}^{-1} | c\ell_3 m_3,, d\ell_4 m_4\rangle
&=\sum_{L=0}^{L_{\max}} \frac{4\pi}{2L+1};
\underbrace{V^{(ee,L)}*{ab,cd}}*{\text{radial from IDA}};
\sum_{M=-L}^{L}G^{L}*{\ell_1 m_1,\ell_3 m_3;M},G^{L,*}*{\ell_2 m_2,\ell_4 m_4;M}.
\end{aligned}
$$

If you only need the **Hartree (J)** potential, you never have to materialize the full four-index tensor: contract the angular Gaunt couplers with the occupied-space density to form **density multipoles**, and combine them with the radial $V^{(ee,L)}$.

---

## 5) HF assembly with $Y_{\ell m}$ (practical outline)

Let the composite basis index be $\mu\equiv(a,\ell,m)$ and $D_{\alpha\beta}$ the one-particle density matrix. The Fock operator is
$$
F_{\mu\nu} ;=; H^{(1)}*{\mu\nu} + J*{\mu\nu} - K_{\mu\nu}.
$$

* **One-electron**: for each $\ell$, assemble $H^{(1)}(\ell)=T^{(0)}+V^{\text{nuc}}+V^{\text{cent}}(\ell)$ (reused across all $m$).
* **Coulomb $J$**: for each $L\le L_{\max}$,

  1. Precompute $V^{(ee,L)}=D^{-1}(M K^{(L)} M^\top)D^{-1}$.
  2. Precompute Gaunt couplers $G^{L}$ and selection masks.
  3. Contract angular density multipoles with $G^{L}$ and combine with the **radial** $V^{(ee,L)}$ to obtain $J_{\mu\nu}$.
* **Exchange $K$ (optional)**: same angular algebra and same radial $V^{(ee,L)}$, but with the usual exchange index pattern.

Because the gausslets are localized and $K^{(L)}$ is smooth, $M K^{(L)} M^\top$ is well conditioned and near-banded; angular selection rules make $G^{L}$ extremely sparse.

---

## 6) Diagnostics and shortcuts

* **Overlap check**: $S\approx I$ on the grid.
* **Virial**: for $\ell=0$, $2T+V^{\text{nuc}}\approx 0$ (one-electron) and the HF virial with $J,K$.
* **Precompute/cache**: $\widetilde\Phi$, $\widetilde D$, $\mathrm{diag}(-Z/r_i)$, $\mathrm{diag}(1/r_i^2)$, the per-$L$ kernels $K^{(L)}$, the row-sums $w^\chi_a$, and Gaunt couplers $G^{L}$.
* **Scaling**: one-electron $\mathcal O(N_r^2 N_Y)$; per-$L$ radial build $\mathcal O(N_r^2 N_g)$ with reuse; angular contractions are tiny due to selection rules.

---

## 7) Minimal equation set (ready to implement)

**Radial blocks**
$$
T^{(0)}=\tfrac12,\widetilde D^\top \widetilde D,\qquad
V^{\text{nuc}}=\widetilde\Phi^\top,\mathrm{diag}(-Z/r_i),\widetilde\Phi,\qquad
V^{\text{cent}}(\ell)=\widetilde\Phi^\top,\mathrm{diag}!\big(\tfrac{\ell(\ell+1)}{2r_i^2}\big),\widetilde\Phi.
$$

**Two-index IDA Coulomb (radial)**
$$
w_a^\chi = \sum_i w_i,\chi_a(r_i),\quad
M=\Phi^\top \mathrm{diag}(w),\quad
K^{(L)}*{ij}=w_i w_j,\frac{r*<^L}{r_>^{L+1}},
$$
$$
V^{(ee,L)}=D^{-1}\big(M K^{(L)} M^\top\big)D^{-1},\qquad D=\mathrm{diag}(w^\chi).
$$

**Angular coupling**
$$
(\mu\nu|\alpha\beta)
=\sum_{L,M}\frac{4\pi}{2L+1},V^{(ee,L)}*{ab,cd},
G^{L}*{\ell_\mu m_\mu,\ell_\alpha m_\alpha;M},
G^{L,*}*{\ell*\nu m_\nu,\ell_\beta m_\beta;M}.
$$

---

### Bottom line

* The Coulomb is assembled **directly in the contracted basis** via the two-index IDA
  $$V^{(ee)}=D^{-1}(M K^{(0)} M^\top)D^{-1},$$
  and generalized to higher $L$ by swapping $K^{(0)}!\to!K^{(L)}$.
* One-electron blocks are simple diagonal weightings on the same grid.
* Angular dependence is cleanly separated by Gaunt couplers.
* No odd gausslets and no analytic radial integrals are needed; accuracy is set by your boundary gausslets and your radial grid.

### Addendum

# Addendum: Linear-time ($O(N_\text{grid})$) build for the two-index IDA Coulomb

For the kernel $K^{(0)}(r,r')=1/\max(r,r')$, you can avoid forming the $N_\text{grid}\times N_\text{grid}$ kernel and replace the double sum by **two single sums** using prefix integrals.

## Continuous identity

For any $f,g$ on $[0,\infty)$,
$$
\iint \frac{f(r),g(r')}{\max(r,r')},dr,dr'
;=; \int_0^\infty \frac{f(r)}{r},G(r),dr ;+; \int_0^\infty \frac{g(r)}{r},F(r),dr,
$$
with prefix integrals
$$
F(r)=\int_0^{r} f(s),ds,
\qquad
G(r)=\int_0^{r} g(s),ds.
$$

This follows by splitting the $(r,r')$–domain into $r'\le r$ and $r'>r$.

## Discrete $O(N_\text{grid})$ rule

Given a radial grid ${(r_i,w_i)}_{i=1}^{N_g}$ and samples $f_i=f(r_i)$, $g_i=g(r_i)$:

1. Prefix sums (one pass)
   $$
   F_i=\sum_{k\le i} w_k f_k,
   \qquad
   G_i=\sum_{k\le i} w_k g_k.
   $$

2. Single-loop accumulation
   $$
   I(f,g);\approx;\sum_{i=1}^{N_g} w_i\left(\frac{f_i}{r_i},G_i ;+; \frac{g_i}{r_i},F_i\right).
   $$

This reproduces the quadrature of the split integral and costs $O(N_g)$ per $(f,g)$ pair.

## Drop-in for the **two-index** IDA Coulomb

Your Coulomb in the contracted basis ${\chi_a}$ is
$$
V^{(ee)}*{ab} ;\approx; \frac{I(\chi_a,\chi_b)}{w_a^\chi,w_b^\chi},
\qquad
w_a^\chi=\sum*{i=1}^{N_g} w_i,\chi_a(r_i).
$$

So for each pair $(a,b)$, set $f=\chi_a$, $g=\chi_b$, compute $F_i,G_i$ once, and evaluate $I(f,g)$ with the single loop above. Overall cost $O(N_b^2 N_g)$ (no $N_g^2$ kernel or matrix–matrix product).

## General multipoles $L>0$

For
$$
K^{(L)}(r,r')=\frac{r_<^{,L}}{r_>^{,L+1}},
$$
the same trick works with **weighted prefix sums**:
$$
\iint f(r)g(r'),K^{(L)}(r,r'),dr,dr'
= \int \frac{f(r)}{r^{L+1}}\Big[\int_0^r g(s),s^L ds\Big] dr
;+;
\int \frac{g(r)}{r^{L+1}}\Big[\int_0^r f(s),s^L ds\Big] dr.
$$

Discrete form:
$$
F^{(L)}*i=\sum*{k\le i} w_k, f_k, r_k^L,
\qquad
G^{(L)}*i=\sum*{k\le i} w_k, g_k, r_k^L,
$$
$$
I_L(f,g)\approx \sum_{i=1}^{N_g} w_i\left(\frac{f_i}{r_i^{L+1}},G^{(L)}_i;+;\frac{g_i}{r_i^{L+1}},F^{(L)}*i\right),
\qquad
V^{(ee,L)}*{ab}\approx \frac{I_L(\chi_a,\chi_b)}{w_a^\chi,w_b^\chi}.
$$

Again $O(N_g)$ per pair.

## Notes

* Handle the origin with your existing cutoff (e.g., ignore points with $r_i<r_\text{cut}$ or use a safe limit for $f_i/r_i$); your boundary gausslets are regular at $r=0$.
* Symmetry check: $I(f,g)=I(g,f)$ to machine precision.
* This replaces the $D^{-1}(M K M^\top)D^{-1}$ route without ever forming $K$ and preserves the exact **two-index** IDA formulation.

# Addendum²: Fast and Accurate Assembly of the Two-Index IDA Coulomb $V_{ee}$

This note restates the two-index integral–diagonal approximation (IDA) in a
form that avoids $N_{\text{grid}}^2$ work and records the discretization that
produced $\sim 10^{-9}$–level HF accuracy for He with modest grids.

---

## 1) Problem statement (two-index IDA)

With orthonormal, localized radial basis functions ${\chi_a(r)}$ and their integral weights
$$
w_a^\chi ;=; \int_0^\infty \chi_a(r),dr,
$$
the two-index IDA Coulomb is
$$
V^{(ee)}_{ab} ;\approx; \frac{1}{w_a^\chi,w_b^\chi};
\iint_0^\infty \frac{\chi_a(r),\chi_b(r')}{\max(r,r')},dr,dr'.
$$

---

## 2) Linear-time split identity

Split the square into $r'\le r$ and $r'>r$:
$$
\int \frac{\chi_a(r),\chi_b(r')}{\max(r,r')} ,dr,dr'
= \int_0^\infty \frac{\chi_a(r)}{r},G_b(r),dr
* \int_0^\infty \frac{\chi_b(r)}{r},F_a(r),dr,
  $$
  with prefix integrals
  $$
  F_a(r)=\int_0^{r} \chi_a(s),ds,
  \qquad
  G_b(r)=\int_0^{r} \chi_b(s),ds.
  $$
  Each $(a,b)$ build now costs $O(N_{\text{grid}})$.

---

## 3) Discretization that actually works

Use a mapped grid from a **uniform** coordinate $u$ (step $h$), with physical radii $r(u)$, Jacobian $J(u)=dr/du$, and accurate **outer** quadrature weights $w_k$ for integrals in $r$.

### 3.1 Prefixes (trapezoid in $u$, second order)

Let $u_k$ be uniform, $r_k=r(u_k)$, $J_k=J(u_k)$, and $\chi_a^k=\chi_a(r_k)$. For each column $a$,
$$
F_a[k] ;\approx; \sum_{m=1}^{k-1} \frac{h}{2},\Big(J_m,\chi_a^m + J_{m+1},\chi_a^{m+1}\Big),
\qquad
F_a[1]=0,
$$
and similarly $G_b$ (or reuse $F$ with $a\leftrightarrow b$).

### 3.2 Outer accumulation (keep your accurate $w_k$)

Define $r_k^{-1}=1/r_k$. Then
$$
I_{ab} ;\approx; \sum_{k} w_k;\Big(\chi_a^k,G_b[k] + F_a[k],\chi_b^k\Big),r_k^{-1}.
$$

### 3.3 Two-sided IDA normalization

Form the vector of basis weights
$$
w^\chi_a ;\approx; \sum_k w_k,\chi_a^k,
$$
and apply the two-sided scaling
$$
V^{(ee)}*{ab} ;\approx; \frac{I*{ab}}{w^\chi_a,w^\chi_b}.
$$
In matrix form: $V^{(ee)} \leftarrow D^{-1} I D^{-1}$ with $D=\mathrm{diag}(w^\chi)$, or elementwise $V^{(ee)} \mathrel{.=/} \big(w^\chi (w^\chi)^\top\big)$.

---

## 4) Extensions and small-$r$ hygiene

* **Multipoles $L>0$** (for spherical-harmonic coupling):
  $$
  K^{(L)}(r,r')=\frac{r_<^{,L}}{r_>^{,L+1}}
  \quad\Rightarrow\quad
  F_a^{(L)}[k] \approx \sum_{m\le k-1} \tfrac{h}{2}\Big(J_m,r_m^{L}\chi_a^m + J_{m+1},r_{m+1}^{L}\chi_a^{m+1}\Big),
  $$
  $$
  I^{(L)}*{ab} \approx \sum_k w_k\left(\frac{\chi_a^k}{r_k^{L+1}},G_b^{(L)}[k] + \frac{\chi_b^k}{r_k^{L+1}},F_a^{(L)}[k]\right),
  \qquad
  V^{(ee,L)}*{ab}\approx \frac{I^{(L)}_{ab}}{w^\chi_a,w^\chi_b}.
  $$

* **Origin:** boundary gausslets satisfy $\chi(0)=0$. Keep a tiny cutoff $r_\text{cut}$ when multiplying by $1/r_k$ (skip the very first point if needed).

* **Centrifugal term:** for $\ell>0$,
  $$
  T_{\text{centrif}}^{(\ell)} ;=; \frac{\ell(\ell+1)}{2};\chi^\top !\left(\frac{w}{r^2}\odot\chi\right),
  $$
  while it vanishes for $\ell=0$.

---

## 5) Validation and scaling

* **Symmetry:** $V^{(ee)}$ is symmetric to machine precision.
* **Convergence:** errors decay predictably with **corespacing** (improve further with Simpson prefixes in $u$ or Richardson in corespacing).
* **Complexity:** $O(N_b^2 N_{\text{grid}})$ without any $N_{\text{grid}}^2$ kernel.

---

## 6) What the numbers showed

With trapezoid prefixes in $u$ and two-sided IDA scaling, your HF for He reached

* $\sim 10^{-9}$ error for corespacing $=0.02$,
* $\sim 10^{-6}$ for corespacing $=0.1$,

consistent with a smooth mapped grid and second-order prefixes, and robust to parameter changes.

---

**TL;DR:** Use the split identity plus prefixes integrated in 
uniform $u$ (trapezoid with $J=dr/du$); keep the outer 
sum with your accurate $w_k$; finish with two-sided 
normalization by $w^\chi$. This 
yields linear-time $V_{ee}$ builds and the accuracy you observed.
