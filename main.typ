#import "@preview/touying:0.5.3": *


#let contentcolor = white
#let bgcolor = black
#set text(fill: contentcolor)
#set page(fill: bgcolor)

#let lwirth-theme(
  content_color,
  bgcolor,
  ..args,
  body,
) = {
  set text(font: "New Computer Modern Sans")
  set text(size: 23pt)

  show: touying-slides.with(
    config-page(paper: "presentation-16-9", margin: 1.5cm, fill: bgcolor),
    config-common(
      slide-fn: slide,
    ),
    ..args,
  )

  body
}

#show: lwirth-theme.with(white, black)



//#set par(justify: true)
//#show math.equation: set text(font: "Noto Sans Math")


#let weblink(a) = text(fill: blue, link(a))

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

#show heading.where(level: 1): set text(35pt)

#set math.mat(delim: "[")
#set math.vec(delim: "[")

#let wedge = math.and

#let avec(a) = math.upright(math.bold(a))
#let vvec(a) = math.accent(math.bold(a), math.arrow)
#let nvec(a) = math.accent(avec(a), math.hat)
#let amat(a) = math.upright(math.bold(a))

#let xv = $avec(x)$
#let ii = $dotless.i$

#let conj(u) = math.overline(u)
#let transp = math.tack.b
#let hert = math.upright(math.sans("H"))

#let clos(a) = math.overline(a)
#let restr(a) = $lr(#a|)$
#let openint(a,b) = $lr(\] #a, #b \[)$

#let inner(a, b) = $lr(angle.l #a, #b angle.r)$


#let grad = $avec("grad")$
#let curl = $avec("curl")$
#let div = $"div"$

#let Hvec = $avec(H)$

#let mesh = $cal(M)$

// Hello and Welcome. Paper Title and Authors
#page(
  background: image("res/bg-vibrant.jpg", width: 100%),
  margin: 2cm,
)[
  #box(
    fill: black.transparentize(20%),
    outset: 20pt,
    radius: 0pt,
  )[
  #set align(center)
    
    #[
      #set text(size: 35pt, weight: "bold")
      #set par(spacing: 5mm)

      Finite Element Exterior Calculus,\
      Homological Techniques, and Applications
    ]

    #[
      #set text(size: 25pt)
      2006
    ]
  
    #[
      #set text(size: 20pt)
      #grid(
        columns: (1fr, 1fr, 1fr),
        align: center,
        [
          #smallcaps[Douglas N. Arnold]\
          University of Minnesota\
          #weblink("mailto:arnold@ima.umn.edu")
        ],
        [
          #smallcaps[Richard S. Falk]\
          Rutgers University\
          #weblink("mailto:falk@math.rutgers.edu")
        ],
        [
          #smallcaps[Ragnar Winther]\
          University of Oslo\
          #weblink("mailto:ragnar.winther@cma.uio.no")
        ]
      )
    ]

    #v(1cm)
    Presented by #smallcaps[Luis Wirth]\
    #weblink("luwirth@ethz.ch")\
    #weblink("lwirth.com")
  ]
]


#slide[
  = Last time: Discrete Exterior Calculus
  #snote[Previous case studies presentation with Felicia Scharitzer]

  - Using theory of Exterior Calculus to solve PDE numerically.
  - PDEs formulated using Differential Forms.
  - DEC tackles *strong form* of PDE.
]


#slide[
  = This time: Finite Element Exterior Calculus
  #snote[Connects to the last talk.]

  - Marriage of Finite Element Method and Exterior Calculus.
  - FEM formulated using Differential Forms
  - FEEC tackles *weak variational form* of PDE.
]

#slide[
  // What was Exterior Calculus again?
  = Exterior Calculus of Differential Forms
  #snote[What was that again?]

  - modern formulation of traditional *vector calculus*
  - *unification* of all kinds of integrals and derivatives
  - *generalization* to *arbitrary dimensions*
  - Calculus on arbitrary Riemannian Manifolds (Differential Geometry)
]


#slide[
  = Differential $k$-Form
  #snote[The central object of study.]
  #v(0.5cm)

  Differential $k$-form $omega in Lambda^k (Omega)$ is $k$-dimensional *integrand*!

  #let this(content) = text(fill: blue.lighten(20%), content)
  #let that(content) = text(fill: red, content)

  $
    // line integral
    omega^1 &=
    quad
    this(3 sin(x)) that(dif x) + this(e^(x y)) that(dif y) + this(log(z)) that(dif z)
    quad
    &in Lambda^1 (Omega)
    \
    // area integral
    omega^2 &=
    quad
    this(5) that(dif x wedge dif y) + this(2) that(dif x wedge dif z) + this(2) that(dif y wedge dif z)
    quad
    &in Lambda^2 (Omega)
    \
    // volume integral
    omega^3 &=
    quad
    this(cos(x y)) that(dif x wedge dif y wedge dif z)
    quad
    &in Lambda^3 (Omega)
  $

  Unifies integration over $k$-dimensional submanifold $M subset.eq Omega$.\
  $
    integral_M omega
  $
]

#slide[
  = Differential $k$-Form

  - $k$-dimensional measurement device $omega in Lambda^k (Omega)$
  - measures $integral_M omega$ $k$-dimensional submanifold $M subset.eq Omega$
  - coefficents determine ruler $omega_p$ at each point $p$
  - tangential $k$-vectors measured at each point

  $omega: p in Omega|-> omega_p quad quad omega_p: (T_p M)^k -> RR$

  
  $
    M subset.eq Omega
    quad quad
    M = "Image" phi
    quad quad
    phi: [0,1]^k -> Omega
    \
    integral_M omega =
    limits(integral dots.c integral)_([0,1]^k) quad
    omega_(vvec(phi)(u))
    ((diff vvec(phi))/(diff t_1) wedge dots wedge (diff vvec(phi))/(diff t_k))
    dif t_1 dots.c dif t_k
  $

]

#slide[ 
  #grid(
    columns: (1fr, 1fr),
    align: center + horizon,
    [
      #set align(left)
      = Bivectors
  
      *Exterior Product* of two vectors $avec(u), avec(v) in RR^3$
      #text(50pt)[$ avec(u) and avec(v) $]
      
      Results in *bivector* or *2-vector*. \
      Oriented plane segment.
    ],
    image("res/bivector-in-3d.svg")
  )
]

#slide[ 
  #grid(
    columns: (1fr, 1fr),
    align: center + horizon,
    [
      #set align(left)
      = Trivectors
  
      Take three vectors $avec(u), avec(v), avec(w) in RR^3$ and get a \
      *trivector* or *3-vector*.
      #text(50pt)[$ avec(u) and avec(v) and avec(w) $]

      It's an oriented plane segment.
    ],
    [
      #image("res/trivector-in-3d.svg")
    ]
  )
]

#slide[
  // Unification of derivatives from vector calculus
  = Derivative Unification
  #snote[One derivative to rule them all.]

  #set align(horizon + center)
  #grid(
    columns: (50%, 50%),
    align: horizon,
    [
      $
        grad f = (dif f)^transp
        \
        "curl" vvec(F) = star d star vvec(F)^transp
        \
        "div" vvec(F) = (star dif vvec(F)^transp)^transp
      $
    ],
    [
      Exterior derivative
      $
        dif: Lambda^k (Omega) -> Lambda^(k+1) (Omega)
      $
    ],
  )
]

#slide[
  = Analysis 2 Theorem Unification
  #snote[One theorem to rule them all.]
  
  #set align(horizon + center)
  #grid(
    columns: (50%, 50%),
    align: horizon,
    [
      // Gradient Theorem
      $
        integral_C grad f dot dif avec(s) =
        phi(avec(b)) - phi(avec(a))
      $
      // Curl Theorem
      $
        integral.double_S curl avec(F) dot dif avec(S) =
        integral.cont_(diff A) avec(F) dot dif avec(s)
      $
      // Divergence Theorem
      $
        integral.triple_V "div" avec(F) dif V =
        integral.surf_(diff V) avec(F) dot nvec(n) dif A
      $
    ],
    [
      Generalized Stokes Theorem
      $
        integral_D dif omega = integral_(diff D) omega
      $
    ]
  )

]

#slide[
  = How to FEEC?

  Unify and generalize:
  - PDE Domain
  - Mesh
  - Sobolev Spaces
  - FE spaces
  - Basis functions
]



#slide[
  = Domain and Mesh

  - PDE Domain $Omega$ is Riemannian Manifold
  - Need to discretize it!
  - Triangulation
  - Simplex: Vertex, Edge, Triangle, Tetrahedron, ...
  - Simplicial Complex
  - Captures Topology, Geometry and Homology of continuous Manifold
  - Structure-Preserving Discretization -> Correct PDE Solution
]

#slide[
  = Beyond scalar-valued FEM

  - NumPDE: Only scalar-valued PDEs $u: Omega -> RR$
  - Only Lagrangian FEM.
  - But there's more!
  - Maxwells Equations!

  $
    &div avec(E) = rho/epsilon_0
    quad quad
    &&curl E = -(diff avec(B))/(diff t)
    \
    &div avec(B) = 0
    quad quad
    &&curl avec(B) = mu_0 (avec(J) + epsilon_0 (diff avec(E))/(diff t))
  $
  
  - Electric Field $avec(E): Omega -> RR^3$ and Magnetic Field $avec(B): Omega -> RR^3$
  - We need function spaces for these vector fields!
  - Weak formulation -> Integrals over $div$ and $curl$!
  - More Sobolov Spaces:
  $
    H    (grad; Omega) = { u: Omega -> RR : integral_Omega norm(grad u)^2 < oo}
    \
    Hvec (curl; Omega) = { avec(u): Omega -> RR^3 : integral_Omega norm(curl u)^2 < oo}
    \
    Hvec (div ; Omega) = { avec(u): Omega -> RR^3 : integral_Omega abs(div u)^2 < oo}
  $
  - We need more FE spaces: Finite dimensiona subspaces of infite-dimensional Function space.

  $
    &H    (grad) &&arrow.squiggly "Lagrangian basis (vertices)"\
    &Hvec (curl) &&arrow.squiggly "Nedelec basis (edges)"\
    &Hvec (div)  &&arrow.squiggly "Raviart-Thomas basis (faces)"\
  $
]

#slide[
  = Unification through FEEC
  - Write Maxwells equation using Differential Forms. Neat!

  $
    &dif E = rho/epsilon_0
    quad quad
    &&dif E = -(diff B)/(diff t)
    \
    &dif B = 0
    quad quad
    &&dif B = mu_0 (J + epsilon_0 (diff E)/(diff t))
  $

  $
    E in Lambda^1(Omega)
    quad
    J in Lambda^2(Omega)
    \
    B in Lambda^2(Omega)
    quad
    rho in Lambda^3(Omega)
  $

  $
    // Faraday 2-form
    F = E wedge dif t + B
    quad in Lambda^2(Omega)
    \
    // current 3-form
    J = rho + J wedge dif t
    quad in Lambda^3(Omega)
  $

  $
    dif F = 0 \
    dif (star F) = J \
  $
  
  - Weak form -> Only integrals over exterior derivative!
  - Unification of Sobolov spaces.
  $
    &H    (grad) &&= H Lambda^0 \
    &Hvec (curl) &&= H Lambda^1 \
    &Hvec (div)  &&= H Lambda^2 \
  $
  - Only have Space of Differential k-Forms.
  - Need $H Lambda^k$-conforming FE spaces!
  - Also only one type of piecewise-linear FE space: Whitney Basis!
  - Whitney 0-forms = Lagrangian Basis (Barycentric)
  - Whitney 1-forms = Nedelec
  - Whitney 2-forms = Raviart-Thomas
  - But in FEM: Of course higher order spaces possible!
  - Space of polynomial differential-forms!!! Compare: DEC only linear.
]


#slide[
  = Example: Hodge-Laplace Problem

  Derive weak form of Hodge-Laplace problem.
  $
    Delta u = f
    \
    (dif delta + delta dif) u = f
    \
    integral_Omega ((dif delta + delta dif) u) wedge star v = integral_Omega f wedge star v
    \
    integral_Omega (dif delta u) wedge star v + integral_Omega (delta dif u) wedge star v = integral_Omega f wedge star v
    \
    integral_Omega (delta u) wedge star (delta v) + integral_Omega (dif u) wedge star (dif v) = integral_Omega f wedge star v
    \
    integral_Omega (delta u) wedge star (delta v) + integral_Omega (dif u) wedge star (dif v) = integral_Omega f wedge star v
    \
    u in H Lambda^k (Omega): quad quad
    inner(delta u, delta v) + inner(dif u, dif v) = inner(f, v)
    quad quad forall v in H Lambda^k (Omega)
  $
]

#slide[
  $
    u_h = sum_(i=1)^N mu_i phi_i
    quad quad
    v_h = phi_j
  $

  $
    u in H Lambda^k (Omega): quad quad
    inner(delta u, delta v) + inner(dif u, dif v) = inner(f, v)
    quad quad forall v in H Lambda^k (Omega)
    \
    sum_(i=1)^N mu_i (integral_Omega (delta phi_i) wedge star (delta phi_j) + integral_Omega (dif phi_i) wedge star (dif phi_j)) = 0
    \
    amat(A) vvec(mu) = 0
  $

  $
    A = C + D
    \
    C = [integral_Omega (delta phi_i) wedge star (delta phi_j)]_(i,j=1)^N
    \
    D = [integral_Omega (dif phi_i) wedge star (dif phi_j)]_(i,j=1)^N
    \
    vvec(phi) = [integral_Omega f wedge star phi_j]_(j=1)^N
  $

]


#slide[
  $L^2$-Inner product on Differential $n$-forms
  $
    inner(alpha, beta) = integral_Omega alpha wedge star beta
  $

  Define coderivative $delta: Lambda^k -> Lambda^(k-1)$ as
  $
    star delta omega = (-1)^k dif star omega
  $
  From this we get integration by parts formulas
  $
    inner(dif omega, eta) = inner(omega, delta eta) + integral_(diff Omega) "Tr" omega wedge "Tr" star eta
  $

  Therefore the coderivative $delta$ is the formal adjoint of the exterior derivative $dif$ w.r.t.
  $L^2$-inner product if $omega$ or $eta$ vanish near the boundary $diff Omega$:
  $
    inner(dif omega, eta) = inner(omega, delta eta)
  $
]

#slide[

= Chain-Complex and Cochain-Complex
De Rham Complex.
Homology and Cohomology
Simplicial Complex vs Differential Complex
Preserve structure!
Actually not only subspace of function space, but actually subcomplex! Preserve all homological structure! => Correct PDE solution.

Develop discretizations that are compatible with geometric, topological and algebraic structures. -> Well-posedness of PDE problem
Discretization is compatible with respect to elliptic complex! Our discret structures reproduce the properties of the continuous structure.

  $
    0 -> C^oo (Omega) limits(->)^grad [C^oo (Omega)]^3 limits(->)^curl [C^oo (Omega)]^3 limits(->)^div C^oo (Omega) -> 0
    \
    0 -> H(grad; Omega) limits(->)^grad Hvec (curl; Omega) limits(->)^curl Hvec (div; Omega) limits(->)^div L^2(Omega) -> 0
    \
    0 -> cal(S)^0_r (mesh) limits(->)^grad bold(cal(N))_(r-1) (mesh) limits(->)^curl bold(cal(R T))_(r-1) (mesh) limits(->)^div cal(S)^(-1)_(r-1) (mesh) -> 0
  $
  $
    0 -> Lambda^0 (Omega) limits(->)^dif Lambda^1 (Omega) limits(->)^dif Lambda^2 (Omega) limits(->)^dif Lambda^3 (Omega) -> 0
    \  
    0 -> H Lambda^0 (Omega) limits(->)^dif H Lambda^1 (Omega) limits(->)^dif H Lambda^2 limits(->)^dif H Lambda^3 (Omega) -> 0
    \
    0 -> cal(P)^-_r Lambda^0 (mesh) limits(->)^dif cal(P)^-_r Lambda^1 (mesh) limits(->)^dif cal(P)^-_r Lambda^2 (mesh) limits(->)^dif cal(P)^-_r Lambda^3 (mesh) -> 0
  $
  $
    0 -> Lambda^0 (Omega) limits(->)^dif Lambda^1 (Omega) limits(->)^dif dots.c limits(->)^dif Lambda^n (Omega) -> 0
    \
    0 -> H Lambda^0 (Omega) limits(->)^dif H Lambda^1 (Omega) limits(->)^dif dots.c limits(->)^dif H Lambda^n (Omega) -> 0
    \
    0 -> cal(P)^-_r Lambda^0 (mesh) limits(->)^dif cal(P)^-_r Lambda^1 (mesh) limits(->)^dif dots.c limits(->)^dif cal(P)^-_r Lambda^n (mesh) -> 0
  $
]

#slide[
  = What is FEEC?

  - Method for creating FEM codes of dimensional generality.
  - Theoretical framework for establishing well-posedness of PDE problems.
  
  - Unification of conforming FEM spaces
  - Generalization to arbitrary dimensions!
  - Respecting Topological and Homological properties of Domain -> Good discretization.


  - My bachelor's thesis: Rust Implementation of Finite Element Exterior Calculus on Coordinate-Free Simplicial Manifolds
]


#slide[
  = Thank you for listening!

  
  - Any questions?
  - Slides on GitHub!
]
