// TODO: mention multi index implementation? Simplices and multivectors are multi indices.
// TODO: better images


#import "@preview/touying:0.7.3": *
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

#import "math.typ": *
#show: math-template

#let fgcolor = white
#let bgcolor = black

#let lwirth-theme(
  fgcolor,
  bgcolor,
  ..args,
  body,
) = {
  set text(font: "New Computer Modern Sans")
  set text(size: 20pt)
  set text(fill: fgcolor)

  show: touying-slides.with(
    config-page(
      paper: "presentation-16-9",
      fill: bgcolor,
      margin: (left: 1.5cm, right: 1.5cm, top: 1.0cm, bottom: 0.1cm),
    ),
    config-common(
      slide-fn: slide,
    ),
    ..args,
  )

  body
}
#show: lwirth-theme.with(white, black)

#show heading.where(level: 1): set text(35pt)
#show heading.where(level: 1): set block(spacing: 10pt)


#let weblink(..args) = text(
    fill: blue,
    link(..args)
  )

#let snote(content) = [
  #set text(20pt)
  #content
]

#slide[
  #set page(background: image("res/bg.png"))

  #set align(center + horizon)
  #block(
    fill: black.transparentize(50%),
    outset: 20pt,
    radius: 0pt,
  )[
    #set align(center)

    #[
      #set text(size: 33pt, weight: "bold")
      #set par(spacing: 5mm)

      Finite Element Exterior Calculus \
      #v(1cm)
      Discretization of Differential Forms \
      for the Numerical Solution of PDEs
    ]

    #v(1cm)
    #block(
      fill: black.transparentize(20%),
      outset: 4pt,
      radius: 5pt,
    )[
      #smallcaps[Luis Wirth]\
      #weblink("http://ethz.lwirth.com")[ethz.lwirth.com] \
      #weblink("mailto:luwirth@ethz.ch")[luwirth\@ethz.ch]
    ]
  ]
]

#slide[
  = What is FEEC?
  #snote[Finite Element Exterior Calculus.]

  A unified framework for PDEs and FEM formulated in differential forms.

  #pause
  *Core idea:* Discretize the whole de Rham complex, not just the individual function spaces. \
  Respect topology, geometry of the domain $->$ obtain correct solutions.

  #pause
  *Big theme:* Cohomology

  #pause
  *Implementation:* Formoniq, a Rust library for FEEC on coordinate-free simplicial complexes in arbitrary dimension (BSc thesis).
]

#slide[
  = The PDE at the heart of FEEC
  #snote[The Hodge–Laplace equation.]

  We want to solve
  #text(50pt)[$
    Delta u = f
  $]
  on a Riemannian manifold $Omega$, where $u, f in Lambda^k (Omega)$ are differential $k$-forms.

  - Generalises scalar ($k=0$) Laplacian.
  - Underlies more sophisticated PDEs like Maxwell's eqs, elasticity eqs, and Stokes flow.
]


#slide[
  = Codifferential $delta: Lambda^k -> Lambda^(k-1)$
  #v(0.5cm)

  Defined as formal $L^2$-adjoint of $dif$ (on closed manifolds):
  $
    delta := dif^*
    quad <==> quad
    inner(dif omega, eta)_(L^2) = inner(omega, delta eta)_(L^2)
  $
  
  #pause
  With boundary we get generalized integration by parts for DF.
  $
    inner(dif omega, eta)_(L^2)
    = inner(omega, delta eta)_(L^2)
    + integral_(partial M) i^* omega wedge i^* (hodge eta)
  $
  
  #pause
  In $RR^3$ (no boundary) these are the adjoint relationships:
  $
    grad^* &= -div
    \
    curl^* &= curl
    \
    div^* &= -grad
  $

  //Using Hodge star we can also write
  //$
  //  delta = plus.minus hodge dif hodge .
  //$
]


#slide[
  = Scalar Laplacian
  #snote[The familiar case, $k = 0$.]

  On $Omega subset.eq RR^n$:
  $
    -Delta_s u = -div grad u
  $

  In the language of forms:
  $
    grad^flat u = dif u
    quad quad quad
    -div X = delta X^flat
  $

  So
  $
    -Delta_s = -div grad quad arrow.squiggly quad Delta^0 := delta dif
  $
]

#slide[
  = The Hodge–Laplace Operator
  #snote[Both compositions, permuted.]

  $
    Delta^k : Lambda^k (Omega) -> Lambda^k (Omega), quad
    Delta^k := dif delta + delta dif
  $

  Symmetric, positive semi-definite, second-order elliptic.

  $
    k &= 0:            quad &&Delta^0 = delta dif             quad &&arrow.squiggly quad -div grad \
    k &in {1,...,n-1}: quad &&Delta^k = dif delta + delta dif quad &&arrow.squiggly quad "Vector Laplacian" \
    k &= n:            quad &&Delta^n = dif delta             quad &&arrow.squiggly quad -div grad \
  $

  All vector-calculus Laplacians are special cases of one operator.
]

#slide[
  = Hodge-Laplace in $RR^3$
  #v(0.5cm)

  //$
  //  0 -> H(grad) ->^grad Hvec(curl) ->^curl Hvec(div) ->^div L^2 -> 0
  //  \
  //  0 <- L^2 <-^(-div) H0vec(div) <-^curl H0vec(curl) <-^(-grad) H0(grad) <- 0
  //$

  #[
    #set align(center)

    #diagram(
      spacing: (2.2em, 1.6em),
      //cell-size: (5em, 1em),
      node-stroke: none,
      edge-stroke: fgcolor,
  
      // Primal (top): left to right
      node((0, 0), $0$),
      node((1, 0), $H(grad)$),
      node((2, 0), $Hvec(curl)$),
      node((3, 0), $Hvec(div)$),
      node((4, 0), $L^2$),
      node((5, 0), $0$),
      edge((0, 0), (1, 0), "->"),
      edge((1, 0), (2, 0), $grad$, "->"),
      edge((2, 0), (3, 0), $curl$, "->"),
      edge((3, 0), (4, 0), $div$, "->"),
      edge((4, 0), (5, 0), "->"),
  
      // Dual (bottom): right to left, but operators horizontally aligned with primal
      node((0, 1), $0$),
      node((1, 1), $L^2$),
      node((2, 1), $H0vec(div)$),
      node((3, 1), $H0vec(curl)$),
      node((4, 1), $H0(grad)$),
      node((5, 1), $0$),
      edge((5, 1), (4, 1), "->"),
      edge((4, 1), (3, 1), $-grad$, "->"),
      edge((3, 1), (2, 1), $curl$, "->"),
      edge((2, 1), (1, 1), $-div$, "->"),
      edge((1, 1), (0, 1), "->"),
    )

    #v(1cm)
  
    #set text(17pt)
    #table(
      columns: 6,
      align: center,
      inset: 6pt,
      //stroke: (x, y) => if y == 0 {(bottom: fgcolor)},
      stroke: fgcolor,
      table.header($k$, $Delta^k = delta dif + dif delta$, $inner(dif, dif) + inner(delta, delta)$, [natural BC], [essential BC], [Name]),
      $0$, $-div grad + 0$, $inner(grad, grad) + 0$, $partial u\/partial n$, [-], [Scalar Neumann],
      $1$, $curl curl - grad div$, $inner(curl, curl) + inner(div, div)$, $curl u times n$, $u dot n$, [Vector],
      $2$, $-grad div + curl curl$, $inner(div, div) + inner(curl, curl)$, $div u$, $u times n$, [Vector],
      $3$, $0 -div grad$, $0 + inner(grad, grad)$, [-], $u$, [Scalar Dirichlet]
    )
  ]

]

#slide[
  = Primal Weak Formulation
  #snote[Integrate against test forms.]

  Multiply by $v in Lambda^k$ and integrate by parts:
  $
    inner((dif delta + delta dif) u, v)
    = inner(dif u, dif v) + inner(delta u, delta v)
  $

  #pause
  We get the following variational equation:
  $
    u in H Lambda^k inter H^* Lambda^k: quad
    inner(dif u, dif v)_(L^2 Lambda^k ) + inner(delta u, delta v)_(L^2 Lambda^k) = inner(f, v)_(L^2 Lambda^k)
    quad forall v in H Lambda^k inter H^* Lambda^k
  $
]

#slide[
  = The Sobolev Spaces
  #snote[Needed for variational formulation]

  $
    H Lambda^k inter H^* Lambda^k
    \
    H Lambda^k (Omega) = { omega in L^2 Lambda^k mid(|) dif omega in L^2 Lambda^(k+1) }
    \
    H^* Lambda^k (Omega) = { omega in L^2 Lambda^k mid(|) delta omega in L^2 Lambda^(k-1) }
  $

  #pause
  In $RR^3$:
  $
    H Lambda^0 &= H (grad; Omega) \
    H Lambda^1 &= H (curl; Omega) \
    H Lambda^2 &= H (div; Omega) \
    H Lambda^3 &= L^2 (Omega) \
  $
]

#slide[
  = Complication: Non-Trivial Kernel
  #snote[Hodge-Laplace is singular.]

  In general $Delta^k$ has nontrivial kernel $ker Delta^k != {0}$ \
  The operator is singular: not invertible, not unconditional existence and uniqueness.

  #pause
  #v(1cm)
  In other words our bilinear form $a(u, v)$, is *not coercive* on $H Lambda^k inter H^* Lambda^k$. \
  It vanishes on the whole kernel. \
  $->$ Lax-Milgram fails.
]

#slide[
  = Harmonic Forms
  #snote[The kernel of $Delta^k$.]

  $
    frak(H)^k (Omega) := ker Delta^k
    = { omega in Lambda^k mid(|) Delta^k omega = 0}
  $

  #pause
  #v(1cm)
  *Hodge's theorem*:
  $
    frak(H)^k tilde.equiv H^k_"dR" (Omega)
  $

  Homology (the study of holes) enters the picture:\
  Each hole gives rise to a unique harmonic representative.

  Number of $k$-dim holes is Betti number
  $
    beta_k = dim H^k_"dR" (Omega)
  $

  $=> $Topology determines kernel of PDE.
]

#slide[
  = de Rham Cohomology

  De Rham Complex
  $
    0 -> Lambda^0 (Omega) limits(->)^dif dots.c limits(->)^dif Lambda^n (Omega) -> 0
    \
    dif compose dif = 0
  $

  Exact ($omega = dif eta$) $=>$ Closed ($dif omega = 0$) \
  Closed ($dif omega = 0$) $arrow.r.double.not$ Exact ($omega = dif eta$)

  #pause
  De Rham cohomology
  $
    H^k_"dR" (Omega) = (ker dif^k)/(im dif^(k-1)) = {omega | dif omega = 0}/{omega | omega = dif eta} = "closed" / "exact"
  $

  *De Rham's theorem*
  $
    H^k_"dR" (Omega) tilde.equiv (H_k (Omega))^*
  $
]

#slide[
  = Fixing the Kernel Issue

  Operator not invertible, but solutions exist *conditionally*!

  *Fredholm alternative* (self-adjoint case):\
  - Existence of solution iff $f perp ker Delta = frak(H)^k$.
  - Uniqueness of solution modulo $ker Delta = frak(H)^k$

  #v(1cm)

  #pause
  Possible fix:
  - Uniqueness: quotient out by harmonics
  - Existance: project out harmonic part of $f$
  $
    u in H Lambda^k inter H^* Lambda^k slash frak(H)^k text(" with ") a(u, v) = inner(f - P_frak(H) f, v).
  $

  #pause
  But there are actually two better fixes.
]

#slide[
  = Uniqueness: Gauge Condition

  Instead of using a quotient space and having to deal with equivalence classes.

  We get uniqueness by enforcing a Gauge
  $
    u perp frak(H)^k
  $
  This picks out a unique representative.

  #v(1cm)

  In weak form
  $
    inner(u, q) &= 0
    quad &&forall q in frak(H)^k
  $
]

#slide[
  = Existence: Hodge Decomposition
  #snote[Every form splits into three $L^2$-orthogonal pieces.]
  
  $
    Lambda^k (Omega) =
    underbrace(dif Lambda^(k-1), "exact"\ omega = dif alpha)
    plus.o
    underbrace(delta Lambda^(k+1), "coexact"\ omega = delta beta)
    plus.o
    underbrace(frak(H)^k, "harmonic"\ Delta omega = 0)
  $

  #pause
  Concretely, every $omega in Lambda^k$ splits as
  $
    omega = dif alpha + delta beta + h, quad h in frak(H)^k.
  $

  #pause
  Compare to Hodge-Laplace
  $
    dif delta u + delta dif u = f
  $
  #pause
  we're missing the harmonic part. Let's introduce it!
  $
    dif delta u + delta dif u + p = f
    quad => quad
    dif delta u + delta dif u = f - p
  $

  Existence!
]

#slide[
  = Well-posed Problem
  #v(1cm)

  We arrive at\
  Find $(u, p) in (H Lambda^k inter H^* Lambda^k) times cal(H)^k$, s.t.
  $
    inner(dif u, dif v) + inner(delta u, delta v) + inner(p, v) &= inner(f, v)
    quad &&forall v in H Lambda^k inter H^* Lambda^k
    \
    inner(u, q) &= 0
    quad &&forall q in frak(H)^k
  $

  #pause
  But:\
  Hard to find FE spaces conforming to both $H Lambda^k$ and $H^* Lambda^k$ simultaneously... \
  #pause
  $->$ We need to get rid of $H^* Lambda^k$! \
  How? Eliminate codifferential!
]

#slide[
  = Mixed Weak Formulation
  #snote[Let's eliminate the codifferential]

  Introduce $sigma = delta u$ as an auxiliary variable.\
  #pause
  We enforce this weakly using the adjoint property.
  $
    inner(sigma, tau) = inner(delta u, tau)
    quad <==> quad
   inner(sigma, tau) = inner(u, dif tau)
  $

  #pause
  #v(1cm)
  And then the codifferential bilinear form becomes
  $
    inner(delta u, delta v) = inner(sigma, delta v) = inner(dif sigma, v)
  $

  No more codifferential! Only exterior derivative.
]

#slide[
  = Final Weak Formulation
  #snote[Built-in Hodge Decomposition]

  Find $(sigma, u, p) in H Lambda^(k-1) times H Lambda^k times frak(H)^k$ such that
  $
    inner(sigma, tau) - inner(u, dif tau) &= 0
    quad &&forall tau in H Lambda^(k-1)
    \
    inner(dif sigma, v) + inner(dif u, dif v) + inner(p, v) &= inner(f, v)
    quad &&forall v in H Lambda^k
    \
    inner(u, q) &= 0
    quad &&forall q in frak(H)^k
  $

  #v(0.3cm)
  - Eq. 1: $sigma$ is the *weak* codifferential of $u$.
  - Eq. 2: Weak Hodge decomposition of $f$: exact $dif sigma$ + coexact $delta dif u$ + harmonic $p$
  - Eq. 3: Gauge constraint $u perp frak(H)^k$

  Saddle-point system. Well-posed by Brezzi inf-sup (Arnold–Falk–Winther).
]

#slide[
  = Discretization
  #snote[How do we make all of this finite-dimensional?]

  Let's start with the manifold. We need a mesh.
]

#slide[
  = Simplices
  #snote[The atoms of the mesh.]

  A Euclidean $k$-simplex is the convex hull of $k+1$ affinely independent points:
  $
    sigma = "conv"{avec(v)_0, ..., avec(v)_k} =
    {
      sum_(i=0)^k lambda^i avec(v)_i
      mid(|)
      quad lambda^i >= 0,
      quad sum_(i=0)^k lambda^i = 1
    }.
  $

  - $0$-simplex: vertex.
  - $1$-simplex: edge.
  - $2$-simplex: triangle.
  - $3$-simplex: tetrahedron.

  #pause
  A abstract combinatorial $k$-simplex is coordinate-free and doesn't need an embedding:
  $
    sigma = [v_0, ..., v_k] in NN^(k+1), quad v_i in NN.
  $

  Two orientations: even and odd permutations of vertex order.
]

#slide[
  = Boundary of a Simplex

  $
    partial [v_0, ..., v_k] = sum_(i=0)^k (-1)^i [v_0, ..., hat(v)_i, ..., v_k]
  $

  Drop one vertex at a time, with alternating sign.

  $
    "edge: " &partial [v_0, v_1] = [v_1] - [v_0] \
    "triangle: " &partial [v_0, v_1, v_2] = [v_1, v_2] - [v_0, v_2] + [v_0, v_1] = [v_0, v_1] + [v_1, v_2] + [v_2, v_0] \
  $

  #pause
  #v(1cm)
  *Key property:*
  #[
  #set text(40pt)
  #set block(spacing: 1em)
  $
    partial compose partial = 0.
  $
  ]
  Reminds you of exterior derivative $dif compose dif = 0$?
]

#slide[
  = Simplicial Complex
  #snote[All subsimplicies included.]

  A simplicial complex $mesh$ is a collection of simplices such that
  - every subsimplex of a simplex in $mesh$ is also in $mesh$,
  - any two simplices intersect in a common subsimplex.
  
  #v(0.3cm)
  #[
    #set align(center)
    #[
      #set block(below: 1pt)
      #image("res/simplices.png", width: 80%)
    ]
    $
      Delta_0 (mesh) #h(1.3cm) limits(<--)^partial #h(1.3cm) Delta_1 (mesh) #h(1.3cm) limits(<--)^partial #h(1.3cm) Delta_2 (mesh) #h(1.3cm) limits(<--)^partial #h(1.3cm) Delta_3 (mesh)
    $
  ]
]

#slide[
  = Simplicial Chain Complex
  #snote[Linear combinations of simplices.]

  $
    C_k (mesh) = { sum_i c_i sigma_i mid(|) sigma_i in Delta_k (mesh), c_i in RR }
  $

  #pause
  This is a complex!
  $
    0 limits(<-)^partial C_0 (mesh) limits(<-)^partial C_1 (mesh) limits(<-)^partial C_2 (mesh) limits(<-)^partial C_3 (mesh) limits(<-)^partial 0
    \
    partial^2 = partial compose partial = 0
  $

  #pause
  *Simplicial homology*\
  Boundary ($c = partial e$) $=>$ cycle ($partial c = 0$)\
  Cycle ($partial c = 0$) $arrow.r.double.not$ boundary ($c= partial e$). Hole!
  $
    H_k (mesh) = (ker partial_k)/(im partial_(k+1)) = {c | partial c = 0}/{c | c = partial d} = "cycles"/"boundaries"
  $
]

#slide[
  = Computational Boundary Operator
  #snote[Just a signed incidence matrix.]

  Order the simplices in each $Delta_k (mesh)$. \
  The signed incidence matrix
  $
    amat(D)_k in {-1, 0, +1}^(N_(k-1) times N_k)
  $
  has entries
  $
    (amat(D)_k)_(i j) = cases(
      +1 quad &"if " sigma_i subset.sq.eq +sigma_j,
      -1 quad &"if " sigma_i subset.sq.eq -sigma_j,
      0  quad &"otherwise."
    )
  $
]

#slide[
  = Mesh Geometry
  #v(0.5cm)

  Two routes:
  - *Extrinsic:* Embedding $Omega arrow.r.hook RR^N$, inherit Euclidean geometry (standard FEM)
  - *Intrinsic:* Riemannian metric directly on each simplex.

  #pause
  Intrinsic is better:
  - Mesh abstract manifolds (flat torus, Klein bottle) without ambient $RR^N$,
  - Extends naturally to Pseudo-Riemannian geometry, e.g.
    - FEEC on 4D spacetimes with Lorentzian metric \
]

#slide[
  = Metric on a Simplex
  #snote[Gram matrix from a basis of the tangent space.]

  A simplex is intrinsically flat $=>$ metric is constant on it. \
  Globally metric is piecewise constant.

  By picking a basis we can represent the metric as a matrix.\
  $=>$ Spanning vectors $avec(e)_i = avec(v)_i - avec(v)_0$
  $
    amat(G)_(i j) = g(avec(e)_i, avec(e)_j) in RR.
  $

  One Gram matrix per top-level simplex.
]

#slide[
  = Intrinsic Geometry via Edge Lengths
  #snote[Some Regge calculus]

  No extrinsic vertex coordinates but intrinsic edge lengths!

  Given edge lengths
  $
    d : Delta_1 (mesh) -> RR^+, quad d_(i j) = d([v_i, v_j]),
  $
  the Gram matrix is determined:
  $
    amat(G)_(i j) = 1/2 (d_(0 i)^2 + d_(0 j)^2 - d_(i j)^2)
  $
  and conversely
  $
    d_(i j) = sqrt(G_(i i) + G_(j j) - 2 G_(i j))
  $

]

#slide[
  = Mesh finished.
  #v(0.5cm)

  - Topology as simplicial complex.
  - Geometry as edge lengths.
]

#slide[
  = Discrete Differential Forms
  #snote[Defined on the mesh]

  - Simplicial Cochains:
    - Discretized Differential Forms
    - DOF coefficients
    - Discrete Exterior Derivative
  - Whitney Forms:
    - Reconstructed Continuous Differential Forms
    - FE basis functions
]

#slide[
  = Simplicial Cochains
  #snote[Discrete differential forms.]

  These are the actual discrete objects.\
  There exists a discrete exterior derivative $dif_h$.

  #pause
  A discrete $k$-form is a $k$-cochain: a real-valued function on $k$-simplices:
  $
    omega_h: Delta_k (mesh) -> RR, quad omega_h in C^k (mesh).
  $

  #pause
  Discretization of continuous forms by the *de Rham map* (integration):
  $
    cal(I): cases(
      Lambda^k (Omega) -> C^k (mesh),
      omega |-> omega_h
    )
    quad "with" quad
    omega_h (sigma) = integral_sigma omega.
  $
]

#slide[
  = Chains and Cochains
  #v(0.5cm)

  $C_k (mesh)$ and $C^k (mesh)$ are *dual* finite-dimensional vector spaces.

  #[
    #set align(center)
    #set text(17pt)
    #table(
      columns: 2,
      align: left,
      stroke: fgcolor,
      inset: 8pt,
      table.header([*Chain* $C_k$ --- subscript], [*Cochain* $C^k$ --- superscript]),
      [Linear combination of $k$-simplices], [Linear functional on $k$-chains],
      [Geometric --- what you integrate over], [Analytic --- what you integrate],
      [Boundary $partial_k: C_k -> C_(k-1)$], [Differential $dif^k: C^k -> C^(k+1)$],
      [Homology $H_k$ --- the holes themselves], [Cohomology $H^k$ --- dual to holes],
      //[Covariant (push-forward)], [Contravariant (pull-back)],
    )
  ]

  #pause
  The *duality pairing* represents discrete integration.
  $
    inner(dot, dot): C_k (mesh) times C^k (mesh) -> RR
    \
    inner(omega, c) = sum_((c_i, sigma_i) in c) c_i thin omega(sigma_i) = integral_c omega
  $
]

#slide[
  = Discrete Exterior Derivative
  #snote[$dif$ is the dual map of $partial$.]

  Stokes' theorem links integration and boundary:
  $
    integral_(partial c) omega = integral_c dif omega.
  $

  #pause
  Under the chain-cochain pairing, $dif$ is the *dual map* of $partial$:
  $
    inner(dif omega, c) = inner(omega, partial c)
    quad quad
    dif^k = partial_(k+1)^*.
  $

  #pause
  *Computationally:* transpose of the incidence matrix,
  $
    dif^k = amat(D)_(k+1)^transp.
  $

  Purely topological. No metric. \
]

#slide[
  = Whitney Forms
  #v(1cm)

  For Galerkin we need a finite-dimensional function subspace. \
  These are the spaces of *Whitney forms*.
  $
    cal(W) Lambda^k (Omega) subset.eq H Lambda^k (Omega)
  $

  #pause
  Unify standard FE families:
  $
    cal(W) Lambda^0 &= cal(P)_1 quad &&"Lagrange (nodal)" \
    cal(W) Lambda^1 &= cal(N) quad &&"Nédélec edge elements" \
    cal(W) Lambda^(n-1) &= cal(R T) quad &&"Raviart–Thomas" \
    cal(W) Lambda^n &= cal(P)_0^"DG" quad &&"piecewise constant"
  $
]


#slide[
  = Whitney Basis
  #v(0.5cm)

  Basis functions for Whitney $k$-forms live on $k$-simplices.
  $
    cal(W) Lambda^k (mesh) = "span" {phi_sigma : sigma in Delta_k (mesh)}
  $
  
  - $cal(W) Lambda^0 (mesh)$ on 0-simplices $Delta_0 (mesh)$
  - $cal(W) Lambda^1 (mesh)$ on 1-simplices $Delta_1 (mesh)$
  - $cal(W) Lambda^2 (mesh)$ on 2-simplices $Delta_2 (mesh)$

  #pause
  *Whitney basis property* (generalization of Lagrange nodal basis property):
  $
    integral_tau phi_sigma = cases(
      +1 quad &"if " tau = +sigma,
      -1 quad &"if " tau = -sigma,
      0  quad &"otherwise."
    )
  $
]

#slide[
  = Whitney LSF
  #v(0.5cm)

  Local shape functions defined using barycentric coordinate functions $lambda_i in Lambda^0 (K)$,
  of top-level simplex $K$. Each associated with a vertex $v_i in K$.

  Whitney $k$-form of a subsimplex $sigma = [v_(i_0), ..., v_(i_k)] subset.eq K$:
  $
    restr(phi_sigma)_K =
    lambda_(i_0 dots i_k) =
    k! sum_(l=0)^k (-1)^l lambda_i_l
    (dif lambda_i_0 wedge dots.c wedge hat(dif lambda)_i_l wedge dots.c wedge dif lambda_i_k)
  $

]

#slide[
  = Whitney 1-Form LSF on Reference 2-Simplex
  #v(0.5cm)

  #figure(
    grid(
      columns: (1fr, 1fr, 1fr),
      rows: 1,
      gutter: 0pt,
      image("res/ref_lsf01.cochain.svg", width: 100%),
      image("res/ref_lsf02.cochain.svg", width: 100%),
      image("res/ref_lsf12.cochain.svg", width: 100%),
    ),
  ) 
  $
    lambda_01 &= (1-y) dif x + x dif y
    quad quad quad
    lambda_02 &= y dif x + (1-x) dif y
    quad quad quad
    lambda_12 &= -y dif x + x dif y
  $
]

#slide[
  = Whitney 1-Form GSF on Triforce Mesh
  #v(1cm)
  
  #figure(
    grid(
      columns: (1fr, 1fr, 1fr),
      rows: 1,
      gutter: 0pt,
      image("res/triforce_gsf01.cochain.svg", width: 100%),
      image("res/triforce_gsf02.cochain.svg", width: 100%),
      image("res/triforce_gsf12.cochain.svg", width: 100%),
    ),
  ) 
]

#slide[
  = Example Whitney 1-Forms on Triforce Mesh
  #snote[Linear Combination]

  #figure(
    grid(
      columns: (1fr, 1fr, 1fr),
      rows: 1,
      gutter: 0pt,
      image("res/triforce_constant.cochain.svg", width: 100%),
      image("res/triforce_div.cochain.svg", width: 100%),
      image("res/triforce_rot.cochain.svg", width: 100%),
    ),
  )

]

#slide[
  = Whitney Reconstruction
  #v(1cm)
  
  Whitney forms are linear combinations of Whitney basis functions. Coefficents are cochain values (DoFs) on simplices.
  $
    u = sum_(sigma in Delta_k (mesh)) u_sigma thin phi_sigma
  $

  #pause
  #v(1cm)
  Reconstruction of discrete cochains by the *Whitney map* (interpolation):
  $
    cal(W) : C^k (mesh) -> cal(W) Lambda^k (Omega)
    \
    cal(W)(omega) = sum_(sigma in Delta_k (mesh)) u_sigma phi_sigma
  $
]

#slide[
  = The de Rham and Whitney Maps
  #v(1cm)

  $
    cal(I) : Lambda^k (Omega) -> C^k (mesh),
    quad quad
    cal(W) : C^k (mesh) -> cal(W) Lambda^k subset.eq H Lambda^k (Omega).
  $

  - Integration: $cal(I)(omega)(sigma) = integral_sigma omega$.
  - Reconstruction: $cal(W)(c) = sum_sigma c(sigma) thin phi^k_sigma$.

  #pause
  #v(1cm)
  *Cochain map property* (maps commute with $dif$):
  $
    cal(I) compose dif = dif compose cal(I),
    quad quad
    cal(W) compose dif = dif compose cal(W).
  $
    
  Because $cal(I)$ and $cal(W)$ are cochain maps, they induce isomorphisms on cohomology.
]

#slide[
  = The Structure Preservation
  #snote[The commutative diagram]

  #diagram(
    spacing: (3em, 2.4em),
    node-stroke: none,
    edge-stroke: fgcolor,

    node((0, 0), $H Lambda^(k) (Omega)$),
    node((1, 0), $H Lambda^(k+1) (Omega)$),
    node((0, 1), $C^k (mesh)$),
    node((1, 1), $C^(k+1) (mesh)$),
    node((0, 2), $cal(W) Lambda^k (mesh)$),
    node((1, 2), $cal(W) Lambda^(k+1) (mesh)$),

    edge((0, 0), (1, 0), $dif$, "->"),
    edge((0, 1), (1, 1), $dif$, "->"),
    edge((0, 2), (1, 2), $dif$, "->"),

    edge((0, 0), (0, 1), $cal(I)$, "->", label-side: left),
    edge((1, 0), (1, 1), $cal(I)$, "->", label-side: left),
    edge((0, 1), (0, 2), $cal(W)$, "->", label-side: left),
    edge((1, 1), (1, 2), $cal(W)$, "->", label-side: left),
  )

]

#slide[
  = The Whitney Complex
  #v(0.5cm)

  FEEC projection operator onto Whitney subspace
  $
    Pi_h := cal(W) compose cal(I) : H Lambda^k -> cal(W) Lambda^k
  $
  
  #v(0.5cm)
  #set align(center)

  #diagram(
    spacing: (3em, 3em),
    node-stroke: none,
    edge-stroke: fgcolor,

    // Top row: de Rham complex
    node((0, 0), $0$),
    node((1, 0), $H Lambda^0 (Omega)$),
    node((2, 0), $H Lambda^1 (Omega)$),
    node((3, 0), $dots.c$),
    node((4, 0), $H Lambda^n (Omega)$),
    node((5, 0), $0$),
    edge((0, 0), (1, 0), "->"),
    edge((1, 0), (2, 0), $dif$, "->"),
    edge((2, 0), (3, 0), $dif$, "->"),
    edge((3, 0), (4, 0), $dif$, "->"),
    edge((4, 0), (5, 0), "->"),

    // Bottom row: Whitney complex
    node((0, 1), $0$),
    node((1, 1), $cal(W) Lambda^0 (mesh)$),
    node((2, 1), $cal(W) Lambda^1 (mesh)$),
    node((3, 1), $dots.c$),
    node((4, 1), $cal(W) Lambda^n (mesh)$),
    node((5, 1), $0$),
    edge((0, 1), (1, 1), "->"),
    edge((1, 1), (2, 1), $dif$, "->"),
    edge((2, 1), (3, 1), $dif$, "->"),
    edge((3, 1), (4, 1), $dif$, "->"),
    edge((4, 1), (5, 1), "->"),

    // Vertical projections
    edge((1, 0), (1, 1), $Pi_h$, "->", label-side: left),
    edge((2, 0), (2, 1), $Pi_h$, "->", label-side: left),
    edge((4, 0), (4, 1), $Pi_h$, "->", label-side: left),
  )

  #v(1.0cm)
  Each square commutes: $quad Pi_h compose dif = dif compose Pi_h$. \
  $=>$ The Whitney complex is a *subcomplex* of the de Rham complex.
  $
    cal(W) Lambda^bullet subset.eq H Lambda^bullet
  $
 
]


#slide[
  = Galerkin Mixed Formulation
  #snote[Replace $H Lambda^k$ by $cal(W) Lambda^k (mesh)$.]

  #v(1cm)
  Find $(sigma_h, u_h, p_h) in cal(W) Lambda^(k-1) times cal(W) Lambda^k times frak(H)_h^k$ such that
  $
    inner(sigma_h, tau) - inner(u_h, dif tau) &= 0
    quad &&forall tau in cal(W) Lambda^(k-1)
    \
    inner(dif sigma_h, v) + inner(dif u_h, dif v) + inner(p_h, v) &= inner(f, v)
    quad &&forall v in cal(W) Lambda^k
    \
    inner(u_h, q) &= 0
    quad &&forall q in frak(H)_h^k
  $

  #v(1cm)
  Convergence guaranteed by FEEC projection $Pi_h = cal(W) compose cal(I)$. \
  Convergence rate $O(h)$ for first-order Whitney.
]

#slide[
  = Matrix Galerkin Mixed Formulation
  #v(0.5cm)

  $
    sum_j sigma_j inner(phi^(k-1)_j,phi^(k-1)_i) - sum_j u_j inner(phi^k_j,dif phi^(k-1)_i) &= 0
    \
    sum_j sigma_j inner(dif phi^(k-1)_j,phi^k_i) + sum_j u_j inner(dif phi^k_j,dif phi^k_i) + sum_j p_j inner(eta^k_j,phi^k_i) &= inner(f,phi^k_i)
    \
    sum_j u_j inner(phi^k_j,eta^k_i) &= 0
  $

  #v(1cm)
  #pause
  
  Given $avec(f) in RR^(N_k)$, find $(avec(sigma),avec(u),avec(p)) in (RR^(N_(k-1)) times RR^(N_k) times RR^(N_k))$ s.t.
  $
    amat(M)^(k-1) avec(sigma) - (amat(dif)^(k-1))^transp amat(M) avec(u) &= 0
    \
    amat(M) amat(dif) avec(sigma) + amat(dif)^transp amat(M)^(k+1) amat(dif) avec(u) + amat(M) amat(H) avec(p) &= avec(f)
    \
    amat(H)^transp amat(M) avec(u) &= 0
  $
]

#slide[
  = Mass Element Matrix
  #snote[Galerkin Geometry]

  Local mass matrix on a top-simplex $K$:
  $
    amat(M)^k_K = [inner(phi^k_i, phi^k_j)_(L^2 Lambda^k (K))]_(i, j)
  $

  #pause
  Expanding the Whitney LSF:
  $
    (M^k_K)_(i,j) &= inner(lambda_(i_0 dots i_k), lambda_(j_0 dots j_k))_(L^2 Lambda^k (K)) \
    &= k!^2 sum_(l,m) (-)^(l+m) innerlines(
      dif lambda_i_0 wedge dots.c wedge hat(dif lambda)_i_l wedge dots.c wedge dif lambda_i_k,
      dif lambda_j_0 wedge dots.c wedge hat(dif lambda)_j_m wedge dots.c wedge dif lambda_j_k,
    )_(Lambda^k)
    integral_K lambda_i_l lambda_j_m vol_g \
  $
]

#slide[
  = This Completes the Construction
  #v(1cm)

  My library can solve
  - Hodge-Laplace
    - Eigenvalue problem
    - Source problem
  - On any mesh
    - any dimension
    - any topology

  #v(1cm)
  Let's look at some results.
]


#slide[
  = Harmonic 1-Forms on $TT^2$, with $beta_1 = 2$
  #v(1cm)

  #figure(
    grid(
      columns: (auto, auto),
      rows: 1,
      gutter: 0pt,
      image("res/torus0.png", height: 80%),
      image("res/torus1.png", height: 80%),
    ),
  )
]

#slide[
  = A Harmonic 1-Forms On a Football Mesh
  #v(1cm)

  #image("res/football.png", height: 80%)
]

#slide[
  = Summary
  #v(0.4cm)

  *FEEC combines three threads:*
  - Differential geometry: manifolds, differential forms, de Rham.
  - Algebraic topology — simplicial complexes, boundaries, homology.
  - Functional analysis — Sobolev spaces, variational formulations.

  #pause
  *Implementation choices in Formoniq:*
  - Arbitary dimensional manifolds with non-trivial topology
  - Coordinate-free geometry
  - Arbitrary form degrees
  - 1st order Whitney forms
  - Hodge–Laplace source and eigenvalue problems.
]

#slide[
  = Extensions
  #v(0.5cm)

  - The Whitney space is the lowest-order *trimmed polynomial* space $cal(W) Lambda^k = cal(P)_1^- Lambda^k$.\
    - Extend to higher-degree polynomial FEEC spaces.
  - Lorentzian geometry for 4D spacetime.
  - Solve other PDEs
    - Maxwell!!
]

#slide[
  #set align(center + horizon)
  #set text(50pt)

  Thank you. \
  #v(0.2cm)
  #weblink("http://ethz.lwirth.com")[ethz.lwirth.com] \
  #weblink("mailto:luwirth@ethz.ch")[luwirth\@ethz.ch]
]

//#slide[
//  = $L^2$-Inner Product
//  #snote[Depends on Riemannian metric $g$ / geometry]
//
//  Pointwise inner product on $Lambda^k_p$ via determinant formula on wedges:
//  $
//    inner(alpha_1 wedge dots wedge alpha_k, beta_1 wedge dots wedge beta_k)_(Lambda^k_p)
//    = det[g^*_p (alpha_i, beta_j)]_(i,j)
//  $
//
//  $L^2$-inner product on $k$-forms (using metric volume form):
//  $
//    inner(omega, eta)_(L^2 Lambda^k) = integral_Omega inner(omega_x, eta_x) thin vol_g
//    quad quad quad
//    vol_g = sqrt(det g) thin dif x^1 wedge dots wedge dif x^n
//  $
//]

//#slide[
//  = Maxwell as Hodge–Laplace
//
//  Cavity resonance problem (source-free, time-harmonic):
//  $
//    E in Lambda^k: quad quad
//    curl curl E = omega^2 E, quad
//    div E = 0, quad
//    restr(E times n)_(partial Omega) = 0.
//  $
//
//  On source-free $E$ we have $div E = 0$ and the full Hodge–Laplace
//  $
//    Delta^1 = dif delta + delta dif quad <-> quad
//    -grad div + curl curl
//  $
//  reduces to $curl curl <-> delta dif$.
//  $==>$ Hodge–Laplace eigenvalue problem on 1-forms.
//
//  Magnetic problems ($B in Lambda^2$), elasticity, Stokes flow all follow the same template at other degrees. FEEC discretises them uniformly.
//]


//#slide[
//  = Simplicial Manifold
//
//  $mesh$ is a *simplicial $n$-manifold* if every $(n-1)$-simplex is the face of exactly one (on boundary) or two (in interior) $n$-simplices.
//
//  For any triangulation of the manifold the simplicial homology of the mesh
//  is isomorphic to the singular homology of the manifold.
//  $
//    H_k^"simp" (mesh) = H_k^"sing" (Omega)
//  $
//
//  We basically assume our manifold from now on is piecewise-linear.
//  Or in other words we replace the smooth manifold with a PL manifold.
//]
