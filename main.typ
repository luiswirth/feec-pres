#import "@preview/touying:0.5.3": *
#import "@preview/fletcher:0.5.2" as fletcher: diagram, node, edge
#import "@preview/tiaoma:0.2.1"

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
  set text(size: 23pt)
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

//#show raw.where(block: false): box.with(
//  //fill: black.lighten(10%),
//  //stroke: black.lighten(10%),
//  inset: (x: 3pt, y: 0pt),
//  outset: (y: 5pt),
//  radius: 4pt,
//)
//#show raw.where(block: true): block.with(
//  fill: black.lighten(5%),
//  inset: (x: 3pt, y: 0pt),
//  outset: 5pt,
//  radius: 4pt,
//)


#page(
  background: image("res/bg-vibrant.jpg", fit: "cover"),
  margin: 2cm,
)[
  #set align(center)
  #block(
    fill: black.transparentize(20%),
    outset: 20pt,
    radius: 0pt,
  )[
    #set align(center)
    
    #[
      #set text(size: 35pt, weight: "bold")
      #set par(spacing: 5mm)

      Rust Implementation of \
      Finite Element Exterior Calculus on \
      Coordinate-Free Simplicial Complexes
    ]

    #v(1cm)
    #smallcaps[Luis Wirth]\
    #weblink("luwirth@ethz.ch")\
    #weblink("ethz.lwirth.com")
  ]
]

#slide[
  = Finite Element Exterior Calculus
  #snote[FEEC]

  - Finite Element Method formulated using Differential Forms
  - Used to solve PDEs by tackling weak variational form.
]



#slide[
  = My Implementation: Formoniq
  #v(0.5cm)

  Goals:
  - Embrace Differential Geometry
  - Arbitrary dimension
  - Arbitrary rank differential forms
  - Non-Trivial Topologies (Cohomology)
  - Intrinsic geometry via Riemannian Metric
]

#slide[
  = PDE Domain as Riemannian Manifold
  #v(1cm)

  - PDE Domain is Riemannian Manifold $Omega$
  - Treating non-trivial domains (arbitrary Betti numbers)
  - Betti number $beta_k =$ number of $k$-dim holes

  #grid(
    columns: (50%, 50%), 
    align: center + horizon,
    image("res/embedding.png"),
    image("res/torus.png"),
  )
]

#slide[
  = Mesh as Simplical Manifold
  #v(1cm)

  #grid(
    columns: (60%, 40%), 
    [
      - Discretize PDE Domain into Mesh
      - Obtain Simplicial Manifold $mesh$ by triangulation of manifold $Omega$

    ], [
      #set align(center + horizon)
      #image("res/moebius.png")
    ]
  )
]


#slide[
  = Coordinate Simplex

  $
    sigma =
    "convex" {avec(v)_0,...,avec(v)_n} =
    {
      sum_(i=0)^n lambda^i avec(v)_i
      mid(|)
      quad lambda^i >= 0,
      quad sum_(i=0)^n lambda^i = 1
    }
  $

  - A 0-simplex is a point.
  - A 1-simplex is a line segment.
  - A 2-simplex is a triangle.
  - A 3-simplex is a tetrahedron.
]

#slide[
  = Coordiante-Free Abstract Simplex
  #snote[Drop the Coordinates!]

  - Only Combinatorics
  - Defines Topology
  - No Geometry

  $
    sigma = [v_0,...,v_n] in NN^(n+1)
    quad quad
    v_i in NN
  $

  Two orientations: Even and Odd Permutations

  Boundary operator
  $
    diff sigma = sum_(i=0)^n (-1)^i [v_0,...,hat(v)_i,...,v_n]
  $
]

#slide[
  = Simplicial Complex
  #snote[All subsimplicies]

  Simplicial Complex $cal(K)$ contains all $k$-dim subsimplicies $Delta_k (mesh), quad k <= n$

  
  #set align(center)
  #set block(below: 1pt)
  #image("res/simplices.png", width: 80%)
  $
    Delta_0 (mesh) #h(1cm) limits(<--)^diff #h(1cm) Delta_1 (mesh) #h(1cm) limits(<--)^diff #h(1cm) Delta_2 (mesh) #h(1cm) limits(<--)^diff #h(1cm) Delta_3 (mesh)
  $
]

#slide[
  = Simplicial Chain Complex

  Chain $c in C_k$ is linear combination of simplicies.

  $
    0 limits(<-)^diff C_0 (mesh) limits(<-)^diff C_1 (mesh) limits(<-)^diff C_2 (mesh) limits(<-)^diff C_3 (mesh) limits(<-)^diff 0
    \
    diff^2 = diff compose diff = 0
  $

  Boundary Operator
  $
    diff_k: C_k (mesh) -> C_(k-1) (mesh)
  $

  Simplicial Homology $cal(H)_k = (ker diff_k)/(im diff_k)$

  $
    "Betti number" beta_k = dim cal(H)_k
  $
]

#slide[
  = Computational Boundary Operator

  Signed Incidence Matrix
  $
    amat(D)_k in {-1,0,+1}^(N_(k-1) times N_k)
  $

  $
    &sigma_i in Delta_(k-1) (mesh) \
    &sigma_j in Delta_k (mesh) \
  $

  $
    (amat(D)_k)_(i j) = cases(
      +1 quad &"if" sigma_i subset.sq.eq +sigma_j,
      -1 quad &"if" sigma_i subset.sq.eq -sigma_j,
      0  quad &"if" sigma_i subset.sq.not plus.minus sigma_j, 
    )
  $
]

#slide[
  = Now Geometry!

  Topology is now defined by the simplicial complex.

  Now let's reintroduce the geometry using a metric.
]

#slide[
  = Metric Tensor as Gram Matrix
  #v(0.5cm)

  $
    g_p: T_p M times T_p M -> RR^+
  $

  Gives Geometry
  - lengths $norm(v)_g = sqrt(g(v, v))$
  - angles $phi(v, w) = arccos((g_p (v, w))/(norm(v)_g norm(w)_g))$.
  
  Discretized as Gram matrix, using basis, induces by chart map.
  $phi: p in Omega subset.eq M |-> (x_1,...,x_n)$
  $
    amat(G)_(i j) = g_p (restr(diff/(diff x^i))_p,restr(diff/(diff x^j))_p)
  $
]

#slide[
  = Metric on Simplex

  Simplex is flat $=>$ Metric is constant on simplex.
  $restr(g)_sigma = "const"$
]

#slide[
  = Derive Metric from Vertex Coordinates
  #snote[Needs Embedding]

  Inherits Geometry from Euclidean Ambient space.

  Spanning vectors $avec(e)_i = avec(v)_i - avec(v)_0 in RR^N$

  $
    amat(G)_(i j) = e_i dot e_j
  $
]

#slide[
  = Derive Metric from Edge Lengths
  #snote[Using Regge Calculus.]

  Edge lengths
  $
    d: Delta_1 (mesh) -> RR^+
    \
    d_(i j) = d([v_i, v_j])
  $


  Metric from Edge Lengths
  $
    amat(G)_(i j) = 1/2 (d_(0 i)^2 + d_(0 j)^2 - d_(i j)^2)
  $

  Edge Lengths from Metric
  $
    d_(i j) = sqrt(amat(G)_(i i) + amat(G)_(j j) - 2 amat(G)_(i j))
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
  Differential Forms on the mesh

  - Simplicial Cochains:
    - Discrete Combinatorial Objects
    - DOF coefficents
  - Whitney Forms:
    - Reconstructed Continuous Objects
    - FE basis functions
]



#slide[
  = Simplicial Cochains
  #snote[Just like in Discrete Exterior Calculus (DEC)]

  Discrete Differential $k$-Form is $k$-Cochain.

  $
    Lambda^k (Omega) arrow.squiggly C^k (mesh)
  $

  Real valued function on all $k$-simplicies $omega: Delta_k (mesh) -> RR$

  Dual to simplicial chains $C^k (mesh) = C_k (mesh)'$
]

#slide[
  = Discretization via Integration

  Discretization of continuous Differential Form via Integration map.
  $
    I: Lambda^k (Omega) -> C^k (mesh)
  $

  $
    I(omega) = (sigma |-> c_sigma) quad "where" quad c_sigma = integral_sigma omega quad forall sigma in Delta_k (mesh)
  $
]

#slide[
  = Discrete Exterior Derivative
  #snote[A little bit of Cochain calculus.]

  $
    dif: Lambda^k (Omega) -> Lambda^(k+1) (Omega)
    quad arrow.squiggly quad
    dif_h: C^k (mesh) -> C^(k+1) (mesh)
  $

  Introduce duality pairing:
  $
    inner(omega, c) := integral_c omega
  $

  Stokes' Theorem. \
  Exterior derivative is adjoint of boundary operator. \
  Coboundary operator.
  $
    integral_c dif omega = integral_(diff c) omega
    quad <==> quad
    inner(dif omega, c) = inner(omega, diff c)
    quad <==> quad
    dif = diff^*
  $
]

#slide[
  = Computational Discrete Exterior Derivative

  Transpose of signed incidence matrix.
  $
    amat(dif)^k = amat(D)_(k+1)^transp
  $

  Purely topological, no metric, no geometry.
]

#slide[
  = Whitney FE Space of Differential Forms
  #v(1cm)
  
  Finite dimensional subspaces of infinite-dimensional function space.
  
  Space of Whitney $k$-forms:\
  Piecewise-linear #text(blue)[coefficents] over cells $Delta_n (mesh)$

  
  $
    cal(W) Lambda^0 (mesh) &=^~ cal(S)^0_1 (mesh) \
    cal(W) Lambda^1 (mesh) &=^~ bold(cal(N)) (mesh) \
    cal(W) Lambda^2 (mesh) &=^~ bold(cal(R T)) (mesh) \
    cal(W) Lambda^3 (mesh) &=^~ cal(S)^(-1)_0 (mesh) \
  $  
]

#slide[
  Whitney Subcomplex
  $
    cal(W) Lambda^0 (mesh) limits(->)^dif dots.c limits(->)^dif cal(W) Lambda^n (mesh)
  $

  De Rham Complex
  $
    H Lambda^0 (Omega) limits(->)^dif dots.c limits(->)^dif H Lambda^n (Omega)
  $

  De Rham Cohomology $cal(H)^k = (ker dif^k)/(im dif^k)$
]

#slide[
  = Whitney Basis
  #v(1cm)

  $
    cal(W) Lambda^k (mesh) = "span" {phi_sigma : sigma in Delta_k (mesh)}
  $
  
  - $cal(W) Lambda^0 (mesh)$ on 0-simplices  $Delta_0 (mesh)$
  - $cal(W) Lambda^1 (mesh)$ on 1-simplicies $Delta_1 (mesh)$
  - $cal(W) Lambda^2 (mesh)$ on 2-simplicies $Delta_2 (mesh)$

  $
    restr(phi_sigma)_K =
    lambda_(i_0 dots i_k) =
    k! sum_(l=0)^k (-1)^l lambda_i_l
    (dif lambda_i_0 wedge dots.c wedge hat(dif lambda)_i_l wedge dots.c wedge dif lambda_i_k)
  $

  Whitney Basis property!
  $
    integral_sigma lambda_tau = cases(
      +&1 quad &"if" sigma = +tau,
      -&1 quad &"if" sigma = -tau,
       &0 quad &"if" sigma != plus.minus tau,
    )
  $
]

#slide[
  = Whitney 1-Forms on 2-Simplex
  #v(0.5cm)

  We have the following formula for Whitney 1-forms.
  $
    lambda_(i j) = lambda_i dif lambda_j - lambda_j dif lambda_i
  $

  For the reference 2-simplex, we get the following Whitney basis 1-forms.
  $
    lambda_01 &= (1-y) dif x + x dif y
    \
    lambda_02 &= y dif x + (1-x) dif y
    \
    lambda_12 &= -y dif x + x dif y
  $
]

#slide[
  = Local Shape Functions
  #figure(
    grid(
      columns: (1fr, 1fr, 1fr),
      rows: 1,
      gutter: 0pt,
      image("res/ref_lambda01.png", width: 100%),
      image("res/ref_lambda02.png", width: 100%),
      image("res/ref_lambda12.png", width: 100%),
    ),
  ) 
]

#slide[
  = Global Shape Functions
  #figure(
    grid(
      columns: (1fr, 1fr, 1fr),
      rows: 1,
      gutter: 0pt,
      image("res/eq_phi01.png", width: 100%),
      image("res/eq_phi02.png", width: 100%),
      image("res/eq_phi12.png", width: 100%),
    ),
  ) 
]


#slide[
  = Whitney Forms
  #figure(
    grid(
      columns: (1fr, 1fr, 1fr),
      rows: 1,
      gutter: 0pt,
      image("res/triforce_constant.cochain.png", width: 100%),
      image("res/triforce_div.cochain.png", width: 100%),
      image("res/triforce_rot.cochain.png", width: 100%),
    ),
  ) 

  $
    u = sum_j u_j phi^k_j
  $
]


#slide[
  = Hodge-Laplace Operator
  #snote[Generalization of scalar Laplacian]

  $
    Delta^k u = f
  $

  Now $u$ and $f$ are Differential $k$-forms.
  $
    u in Lambda^k (Omega), f in Lambda^k (Omega)
  $

  Involves exterior derivative and codifferential.
  $
    Delta^k: Lambda^k (Omega) -> Lambda^k (Omega)
    \
    Delta^k := delta dif + dif delta =  delta^(k+1) dif^k + dif^(k-1) delta^k
  $
]


#slide[
  = Hodge-Laplace in $(RR^3; times, dot)$
  #v(0.5cm)

  $
    0 -> H(grad) ->^grad Hvec(curl) ->^curl Hvec(div) ->^div L^2 -> 0
    \
    0 <- L^2 <-^(-div) H0vec(div) <-^curl H0vec(curl) <-^(-grad) H0(grad) <- 0
  $
  
  #[
    #set align(center)
    #set text(17pt)
    #table(
      columns: 6,
      align: center,
      //stroke: (x, y) => if y == 0 {(bottom: fgcolor)},
      stroke: fgcolor,
      table.header($k$, $Delta^k = delta dif + dif delta$, $tilde(Delta)^k = inner(dif, dif) + inner(delta, delta)$, [natural BC], [essential BC], $V^(k-1) times V^k$),
      $0$, $-div grad + 0$, $inner(grad, grad) + 0$, $diff u\/diff n$, [-], $H(grad)$,
      $1$, $curl curl - grad div$, $inner(curl, curl) + inner(div, div)$, $curl u times n$, $u dot n$, $H(grad) times Hvec(curl)$,
      $2$, $-grad div + curl curl$, $inner(div, div) + inner(curl, curl)$, $div u$, $u times n$, $Hvec(curl) times Hvec(div)$,
      $3$, $0 -div grad$, $0 + inner(grad, grad)$, [-], $u$, $H(div) times L^2$,
    )
  ]

  $k=0$: Scalar Laplacian w/ Neumann B.C. \
  $k=1$: Vector Laplacian w/ Magnetic B.C. \
  $k=2$: Vector Laplacian w/ Electric B.C. \
  $k=3$: Scalar Laplacian w/ Dirichlet B.C. \
]


#slide[
  = Weak Variational Form & $L^2$ Inner Product
  #v(1cm)

  The inner product is defined as
  $
    inner(omega, eta)_(L^2 Lambda^k (Omega))
    := integral_Omega inner(omega_x, eta_x)_(Lambda^k) vol_g
    = integral_Omega omega wedge hodge eta
  $
  #pause

  Sobolev Space
  $
    H Lambda^k (Omega) = { omega in L^2 Lambda^k (Omega) mid(|) dif omega in L^2 Lambda^(k+1) (Omega) }
  $

  $
    cal(W) Lambda^k (Omega) subset.eq H Lambda^k
  $
]

#slide[
  = Integrate against Test function

  Take strong form and form $L^2$-inner product with test function $v$
  $
    Delta u = f
  $

  We obtain the variational equation
  $
    u in H Lambda^k (Omega): quad quad
    inner(Delta u, v)_(L^2 Lambda^k (Omega)) = inner(f, v)_(L^2 Lambda^k (Omega))
    quad quad forall v in H Lambda^k (Omega)
  $
]

#slide[
  = Mixed Hodge-Laplace Source Problem
  #v(0.5cm)

  $
    Delta^k u = dif delta u &+ delta dif u = f
    \
    &arrow.b
    \
    sigma &= delta u
    \
    dif sigma + delta dif u &= f - p
    \
    u &perp frak(H)^k \
  $

  #pause

  Given $f in L^2 Lambda^k$, find $(sigma,u,p) in (H Lambda^(k-1) times H Lambda^k times frak(H)^k)$ s.t.
  $
    inner(sigma,tau) - inner(u,dif tau) &= 0
    quad &&forall tau in H Lambda^(k-1)
    \
    inner(dif sigma,v) + inner(dif u,dif v) + inner(p,v) &= inner(f,v)
    quad &&forall v in H Lambda^k
    \
    inner(u,q) &= 0
    quad &&forall q in frak(H)^k
  $
]

#slide[
  = Galerkin Mixed Hodge-Laplace Source Problem
  #v(0.5cm)

  Given $f in L^2 Lambda^k$, find $(sigma,u,p) in (H Lambda^(k-1) times H Lambda^k times frak(H)^k)$ s.t.
  $
    inner(sigma,tau) - inner(u,dif tau) &= 0
    quad &&forall tau in H Lambda^(k-1)
    \
    inner(dif sigma,v) + inner(dif u,dif v) + inner(p,v) &= inner(f,v)
    quad &&forall v in H Lambda^k
    \
    inner(u,q) &= 0
    quad &&forall q in frak(H)^k
  $

  $
    sum_j sigma_j inner(phi^(k-1)_j,phi^(k-1)_i) - sum_j u_j inner(phi^k_j,dif phi^(k-1)_i) &= 0
    \
    sum_j sigma_j inner(dif phi^(k-1)_j,phi^k_i) + sum_j u_j inner(dif phi^k_j,dif phi^k_i) + sum_j p_j inner(eta^k_j,phi^k_i) &= inner(f,phi^k_i)
    \
    sum_j u_j inner(phi^k_j,eta^k_i) &= 0
  $
]

#slide[
  = Galerkin Hodge-Laplace Source Problem
  #v(0.5cm)

  Given $f in L^2 Lambda^k$, find $(sigma,u,p) in (H Lambda^(k-1) times H Lambda^k times frak(H)^k)$ s.t.
  $
    inner(sigma,tau) - inner(u,dif tau) &= 0
    quad &&forall tau in H Lambda^(k-1)
    \
    inner(dif sigma,v) + inner(dif u,dif v) + inner(p,v) &= inner(f,v)
    quad &&forall v in H Lambda^k
    \
    inner(u,q) &= 0
    quad &&forall q in frak(H)^k
  $

  
  Given $avec(b) in RR^(N_k)$, find $(avec(sigma),avec(u),avec(p)) in (RR^(N_(k-1)) times RR^(N_k) times RR^(N_k))$ s.t.
  $
    amat(M)^(k-1) avec(sigma) - (amat(dif)^(k-1))^transp amat(M) avec(u) &= 0
    \
    amat(M) amat(dif) avec(sigma) + amat(dif)^transp amat(M)^(k+1) amat(dif) avec(u) + amat(M) amat(H) avec(p) &= avec(b)
    \
    amat(H)^transp amat(M) avec(u) &= 0
  $

  $
    amat(M)^k
    = inner(phi^k_i, phi^k_j)_(L^2 Lambda^k (Omega))
  $
]

#slide[
  = Mass Element Matrix
  #v(0.5cm)

  $
    M^k_K &= inner(lambda_(i_0 dots i_k), lambda_(j_0 dots j_k))_(L^2 Lambda^k (K)) \
    &=
    k!^2 sum_(l=0)^k sum_(m=0)^k (-)^(l+m) innerlines(
      lambda_i_l (dif lambda_i_0 wedge dots.c wedge hat(dif lambda)_i_l wedge dots.c wedge dif lambda_i_k),
      lambda_j_m (dif lambda_j_0 wedge dots.c wedge hat(dif lambda)_j_m wedge dots.c wedge dif lambda_j_k),
    )_(L^2 Lambda^k (K)) \
    &= k!^2 sum_(l,m) (-)^(l+m) innerlines(
      dif lambda_i_0 wedge dots.c wedge hat(dif lambda)_i_l wedge dots.c wedge dif lambda_i_k,
      dif lambda_j_0 wedge dots.c wedge hat(dif lambda)_j_m wedge dots.c wedge dif lambda_j_k,
    )_(Lambda^k)
    integral_K lambda_i_l lambda_j_m vol_g \
  $
]

#slide[
  = Harmonic Forms 
  #v(0.5cm)

  $
    frak(H)^k = ker Delta = { u in Lambda^k mid(|) Delta u = 0}
  $

  Concrete Representatives of Cohomology group
  $
    frak(H)^k =^~ H^k
    \
    dim frak(H)^k = beta_k
  $

  Compute Harmonic Forms via Eigenvalue problem.
  $
    Delta u = 0
    quad <==> quad
    Delta u = lambda u "with" lambda = 0
  $
]


//#slide[
//  = 1-Form EVP on Annulus
//  #v(0.5cm)
//
//  #figure(
//    grid(
//      columns: (1fr, 1fr, 1fr),
//      rows: 1,
//      gutter: 0pt,
//      image("res/evp0.png", width: 100%),
//      image("res/evp5.png", width: 100%),
//      image("res/evp6.png", width: 100%),
//    ),
//  ) 
//]


#slide[
  = Harmonic 1-forms on Torus
  #v(0.5cm)

  #figure(
    grid(
      columns: (1fr, 1fr),
      rows: 1,
      gutter: 0pt,
      image("res/torus_eigen0_full.png", width: 100%),
      image("res/torus_eigen1_full.png", width: 100%),
    ),
  ) 
]

//#slide[
//  //#set page(background: image("res/bg-vibrant.jpg", width: 100%))
//
//  = Thank you for listening!
//
//  #set align(center + horizon)
//  #block()[
//    #set align(center)
//    #set par(spacing: 10pt)
//
//    Presentation Slides
//    #tiaoma.qrcode("https://github.com/luiswirth/feec-pres",
//      options: (
//        scale: 4.0,
//        fg-color: fgcolor,
//        bg-color: bgcolor,
//      )
//    )
//    #weblink("https://github.com/luiswirth/feec-pres", "github:luiswirth/feec-pres")
//  ]
//]
