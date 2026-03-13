# Boundary Gausslets on the Half-Line, Erf-Grids, and Hydrogen Tests — A Consolidated Report

*(motivations, equations, logic, implementation notes, diagnostics; code names are given for reference only)*

---

## 1) Why boundary gausslets?

Uniform gausslets are short linear combinations of primitive Gaussians that behave like compact, nearly orthonormal, and high-order quadrature‐friendly basis functions on a uniform grid. On the radial half-line $r\in[0,\infty)$ we need basis functions that:

* enforce the $s$-wave boundary condition $u(0)=0$ (or, equivalently, make the mode values vanish at $r=0$),
* retain even-moment content near the origin (so we don’t cripple accuracy),
* remain localized and well-conditioned, and
* support sparse operators and fast quadrature on simple grids.

Two constructions meet these needs. Below we emphasize the **δ-surgery** route we actually used here (and briefly contrast the polynomial route for context), then show how we **localize by COMX**, assemble operators with **erf-type grids**, and verify everything with **hydrogen**.

*(For a complementary write-up—history, motivations, and variants—see the earlier report the user provided.)* 

---

## 2) δ-surgery: remove precisely the “value-at-0” direction

Let ${\phi_j(x)}$ be the full-line, (numerically) orthonormal gausslets on centers $x_j=jH$. Define the “value vector” at the origin
$$
v_j=\phi_j(0),\qquad |v|^2=\sum_j v_j^2.
$$

Project away the one coefficient-space direction that controls $\psi(0)$:
$$
P_0 ;=; I - \frac{v,v^\top}{v^\top v},
\qquad
\tilde\phi_j(x) ;=; \phi_j(x) - \frac{v_j}{|v|^2}\sum_k v_k,\phi_k(x).
$$
Then $\tilde\phi_j(0)=0$ for all $j$, yet all even derivatives at $0$ are **preserved** (we did *not* oddize).

Restrict to the half-line and keep a small number of negative centers (their positive tails sharpen near-origin even content):
$$
\phi_j^+(r)=\tilde\phi_j(r),\Theta(r)\quad (r\ge 0).
$$

**Half-line metric and position.** Build
$$
S^+*{ij}=\int_0^\infty \phi_i^+(r),\phi_j^+(r),dr,\qquad
X^+*{ij}=\int_0^\infty \phi_i^+(r),r,\phi_j^+(r),dr.
$$

**Symmetric orthonormalization on $[0,\infty)$.** With $S^{+^{-1/2}}$,
$$
\Phi_c(r)=\sum_j \phi_j^+(r)\big(S^{+^{-1/2}}\big)*{j\cdot}
\quad\Longrightarrow\quad
\langle \Phi*{c,m},\Phi_{c,n}\rangle_+=\delta_{mn}.
$$

**COMX localization.** Diagonalize $X$ in the orthonormal frame,
$$
X_c=S^{+^{-1/2}}{}^\top,X^+,S^{+^{-1/2}}=U,\Lambda,U^\top,
$$
and rotate:
$$
\boxed{\ \psi_m(r)=\sum_n \Phi_{c,n}(r),U_{n m}\ },\qquad
\Lambda_{mm}=\langle \psi_m|r|\psi_m\rangle_+.
$$
The $\psi_m$ are orthonormal, satisfy $\psi_m(0)=0$, and are sharply localized near their COMX centers.

**Seed pre-normalization.** Because we **ultimately use the half-line**, it is advantageous to rescale the seeds first:
$$
D_j=\frac{1}{\sqrt{\int_0^\infty \phi_j(r)^2,dr}},\qquad
\phi_j \leftarrow D_j,\phi_j.
$$
We fold this $D$ into the final seed→mode map (so the cache is “ready-to-use” on $[0,\infty)$).

> **What we cache (JLD2):** for each $j_{\rm neg}\in{2,3,4,5,6}$ we store one group with
> – trimmed index vector $j\in[-j_{\rm neg},\dots,J_{\rm pos}]$,
> – the **final** seed→mode matrix $C$ (already includes pre-normalization and all orthogonalizations), and
> – the list of COMX centers on the **uniform** axis (to be mapped later).
> File: `BoundaryGausslets.jld2` (built by `make_boundary_cache.jl`).

*(A polynomial augmentation pathway—oddize then add orthonormal residuals $r^{2},r^{4},\dots$—is also viable; we skip details since the δ-surgery route is what we used operationally here.)* 

---

## 3) Erf-type grids on $[0,\infty)$ and a composed, mapped grid

Simple, **uniform** quadrature works very well for smooth, same-scale integrands. On the half-line we want very fine resolution near $r=0$ (to accommodate $1/r,1/r^2$ when multiplied by functions that vanish at $0$) while transitioning smoothly to **constant spacing** for large $r$.

We therefore use a **uniform grid in an auxiliary coordinate $s$** with a smooth mapping $r=r(s)$ of “erf-type,”
$$
r(s)=h,s,g(s),\qquad g(s)=\tfrac{1}{2}\big[1+\mathrm{erf}\big(\tfrac{s-s_0}{\sigma}\big)\big],
$$
so that $r'(s)!\to!0$ near the left endpoint and $r'(s)!\to!h$ for large $s$. We build weights with the **trapezoidal rule** in $s$:
$$
\int_0^\infty f(r),dr ;\approx; \sum_{i=1}^N f\big(r(s_i)\big),r'(s_i),\Delta s,
\quad s_i=i,\Delta s,\ \Delta s=h.
$$
This yields **exponential‐accuracy behavior** for smooth $f$, while keeping a **uniform** grid in practice.

For radial Schrödinger work we then **compose** this erf-grid with the standard **inverse-square-root mapping** used for gausslet distortions:
$$
u = \frac{1}{\omega},r + F^{-1}(r;,c,s),
\qquad \text{(schematically; actual algebra in `CoordinateMapping.jl`),}
$$
and use the standard orthonormal **distortion factor** so that a uniform $u$-grid gausslet $\phi_j(u)$ becomes the mapped, orthonormal
$$
\chi_j(r)=\frac{\phi_j!\big(u(r)-jH\big)}{\sqrt{J(r)}},
\qquad J(r)=\frac{du}{dr}.
$$
Derivatives follow by the chain rule,
$$
\chi_j'(r)
= \frac{1}{\sqrt{J(r)}},\phi'_j!\big(u(r)-jH\big),u'(r) ;-; \frac{J'(r)}{2,J(r)},\chi_j(r).
$$
We evaluate $\phi_j$ and $\phi'_j$ by the gausslet library’s primitive routines (unit spacing), and we supply $u'(r),J(r),J'(r)$ from `CoordinateMapping.jl`.

> **Two grids used here**
>
> * `Erfgrid.jl`: plain erf-grid for constructing/computing half-line seed matrices and checks.
> * `ErfMappedGrid.jl`: erf-grid **composed** with the coordinate mapping (used in `solveHatom.jl`).

**Accuracy checks.** With $f(r)=(r e^{-r})^2$ we verified rapid convergence of
$$
\int_0^\infty f(r),r^p,dr \quad (p=-2,-1,0,1),
$$
e.g., with $\sigma=3$ we observed $\lesssim10^{-12}$ relative errors at moderate steps $h$ (see the development logs). Near $r=0$ the basis satisfies $\chi_j(0)=0$, so integrals of $1/r$ and $1/r^2$ stay **absolutely convergent** and numerically stable on the erf-grid.

---

## 4) Assembling the boundary gausslets and trimming their support

Given the cached $(C,,j\text{s})$ for a chosen $j_{\rm neg}$:

1. **Extend** the $j$ range if needed so that modes cover the requested $R_{\max}$ after mapping; beyond the cached right edge we append **pure mapped gausslets** (identity coefficients).
2. **Map** the cached COMX centers to physical $r$ with the coordinate mapping (centers are stored on the uniform axis).
3. **Trim locality (optional)** by inspecting coefficient vectors of each mode $m$ and defining the *tail norm*
   $$
   |C_{\cdot,m}|^2_{\text{tail}}=\sum_{j\neq j^\ast(m)}C_{j,m}^2,
   \qquad j^\ast(m)=\arg\max_j |C_{j,m}|,
   $$
   then cutting the support where the tail falls below a user tolerance (we used $10^{-12}$).
   In practice we found that choosing a **fixed** $J_{\rm pos}=24$ across $j_{\rm neg}=2..6$ is simple and safe.

**How large should $j_{\rm neg}$ be?** From energy errors for H at fixed $R_{\max}=16$, $h_{\rm grid}=5\times 10^{-4}$, we observed

```
jneg=2 → 0.5 + E1 ≈ 1.55e-5
jneg=3 → 0.5 + E1 ≈ 4.81e-6
jneg=4 → 0.5 + E1 ≈ 6.78e-7
jneg=5 → 0.5 + E1 ≈ 3.17e-8
jneg=6 → 0.5 + E1 ≈ 1.84e-9
```

Moving from $j_{\rm neg}=4$ to $6$ costs only a few more modes but materially improves near-origin completeness; we recommend **$j_{\rm neg}=6$** for production.

---

## 5) Hydrogen (ℓ=0 and ℓ>0): operators, energies, virial

**Basis vectors.** We combine

* **boundary modes** $\psi_m$ (columns of $C$ applied to seed gausslets, then mapped with the distortion factor) up to the cached right edge,
* followed by **pure mapped gausslets** beyond the edge so overall coverage reaches the requested $R_{\max}$.

**Grid assembly.** On the erf-mapped grid ${(r_i,w_i)}$,

* Overlap: $S_{mn}\approx\sum_i w_i,\chi_m(r_i)\chi_n(r_i)$ (checked $|S-I|\sim10^{-10}$–$10^{-12}$ for well-resolved grids).
* Kinetic:
  $$
  T_{mn}\approx \tfrac12\sum_i w_i,\chi_m'(r_i),\chi_n'(r_i).
  $$
* Coulomb:
  $$
  V^{(C)}_{mn}\approx -Z\sum_i w_i,\frac{\chi_m(r_i)\chi_n(r_i)}{r_i}.
  $$
* Centrifugal ($\ell>0$):
  $$
  V^{(\ell)}_{mn}\approx \frac{\ell(\ell+1)}{2}\sum_i w_i,\frac{\chi_m(r_i)\chi_n(r_i)}{r_i^2}.
  $$

We **solve a standard eigenproblem** (since $S\approx I$ numerically):
$$
H=T+V^{(C)}+V^{(\ell)}\quad\Rightarrow\quad H,c=E,c.
$$

**Results (typical, well-resolved run).** With $(\text{corespacing},\omega,s,R_{\max},h_{\rm grid},j_{\rm neg})=(0.2,10,0.2,60,5\times10^{-4},6)$ and $N!=!36$:

* $|S-I|\approx1.8\times10^{-10}$,
* $E_{1s}\approx-0.49999999818$ (variational error $\approx 1.82\times10^{-9}$),
* $E_{2s}\approx-0.125000000$,
* higher $s$ levels spot-on to $\sim10^{-9}$.

**Virial checks.**
For $\ell=0$ (no centrifugal), the usual virial $2\langle T\rangle+\langle V^{(C)}\rangle=0$ holds to $\sim10^{-9}$.
For $\ell>0$, split $T=T_{\rm rad}+T_{\rm cent}$ and verify either equivalent form:
$$
2\big\langle T_{\rm rad}+T_{\rm cent}\big\rangle + \langle V^{(C)}\rangle \approx 0,
\qquad
2\langle T_{\rm rad}\rangle \approx -\langle V^{(C)}\rangle -2\langle V_{\rm cent}\rangle,
$$
again satisfied to $\sim10^{-10}$–$10^{-11}$ across $\ell=1..4$ in our tests.

---

## 6) Practical recipe (what we actually run)

* **Cache once** (JLD2): build `BoundaryGausslets.jld2` for $j_{\rm neg}=2..6$ with a fixed $J_{\rm pos}=24$. The stored $C$ already includes **seed pre-normalization** and **half-line orthonormalization** and **COMX rotation**.
* **Load for a run**: choose $j_{\rm neg}$ (we default to 6), map centers through the **current coordinate mapping**, append pure mapped gausslets until the mapped center list exceeds your $R_{\max}$, and truncate modes by $R_{\max}$ if desired.
* **Grid**: use the **erf-mapped** grid with parameters $(h_{\rm grid},\sigma,s_0)$; we typically use $h_{\rm grid}=10^{-2}$ for diagnostics, $5\times10^{-4}$ or $10^{-4}$ for production; $\sigma=2$–$3$, $s_0\approx6.5$ are excellent defaults.
* **Quick health checks**:
  (i) overlap $|S-I|$,
  (ii) exponential moments $\int_0^\infty e^{-r/a}r^k,dr$ sampled on the grid,
  (iii) one‐function plots of $\chi_m(r)$ and $\chi_m'(r)$ near $r=0$ (should vanish linearly and stay smooth).
* **Solve** $H$ and (optionally) print **virial diagnostics** using either of the $\ell>0$ forms above.

---

## 7) What we learned (and pitfalls we fixed)

* **Integrals on $[0,\infty)$ need care.** Naïve `quadgk` over a vast interval can silently miss tails. Splitting integrals or, better, **moving to an explicit grid** made issues transparent and reproducible. With the erf-grid we achieved machine-precision overlaps for well-resolved runs.
* **Seed pre-normalization matters.** Pre-scaling seeds by their half-line norms tightens conditioning and simplifies caching (the stored $C$ can be used directly on $[0,\infty)$).
* **COMX centers are mapping-agnostic.** Cache COMX centers on the **uniform** axis; always map them at use time. (We had one misstep early by treating cached centers as physical $r$—fixed.)
* **Sign fixing must use the correct measure.** We determine the sign of each mode from $\int_0^\infty \psi_m(r),dr$ **on the same grid/measure** used for operator assembly (not from uniform-axis coefficients). Once we aligned this, near-origin behavior/plots were clean.
* **Where to stop on the left.** Increasing $j_{\rm neg}$ from 4→6 costs little (a handful of modes) and clearly reduces the $1s$ energy error by ~3 orders of magnitude in our standard test.

---

## 8) Parameter defaults that worked well

* **Cache build:** uniform spacing $H=1$, $J_{\rm pos}^{\rm build}=40$ (safe headroom), $j_{\rm neg}\in[2,6]$, seed pre-norm on $[0,\infty)$, δ-surgery projector with a **single** pivot removal, half-line S-ortho, COMX rotation, stored $C$ includes all factors.
* **Run time:** $(\text{corespacing},\omega,s)=(0.2,10,0.2)$ unless physics dictates otherwise; erf-mapped grid with $h_{\rm grid}=5\times10^{-4}$ (or $10^{-4}$), $\sigma=2$–$3$, $s_0=6.5$, and $R_{\max}$ chosen so tails are negligible for the target $n,\ell$.
* **Left/right support:** $j_{\rm neg}=6$ and a fixed $J_{\rm pos}=24$ for the cached part; append mapped gausslets to reach $R_{\max}$.

---

## 9) Summary

* **Boundary gausslets via δ-surgery**: remove only the “value-at-0” direction in coefficient space, preserve even moments, restrict to $[0,\infty)$, S-orthonormalize, then **COMX** to localize.
* **Erf-mapped grids**: a simple, uniform $s$-grid plus a smooth $r(s)$ delivers **exponential accuracy** for smooth integrands, resolves $1/r$ and $1/r^2$ safely (since $\chi(0)=0$), and keeps assembly straightforward.
* **Hydrogen**: with $j_{\rm neg}=6$, modest $N$, and sensible grid parameters we obtained machine-level agreement for $s$ and excellent results for $\ell>0$; the **virial** holds to $\sim10^{-10}$–$10^{-11}$.
* **Caching**: JLD2 groups per $j_{\rm neg}$ with a consistent right edge ($J_{\rm pos}=24$) make subsequent radial runs trivial—just load, map, append, assemble, and solve.

---

### References to the earlier write-up

A prior document the user shared covers motivations, alternatives (e.g., polynomial augmentation), and complementary details. We have integrated its main ideas here while focusing on the δ-surgery + erf-grid path actually implemented and tested in this work. 

---

*Notes on style and conventions used here:* we avoid Unicode Greek in names and equations that need to be typed verbatim in code; when we refer to code, we use ASCII identifiers (`chi`, `chip`, `rho`, `rhop`, etc.), and we keep matrix multiplications explicit (prefer `.*` for diagonals) to match the project’s coding preferences.

