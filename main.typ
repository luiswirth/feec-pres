#import "@preview/touying:0.5.3": *
#import "@preview/tiaoma:0.2.1"


#let fgcolor = white
#let bgcolor = black
#set text(fill: fgcolor)
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
#let H0 = $limits(H)^circle.stroked.small$

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
  = Finite Element Exterior Calculus
  #snote[Connects to last case study with Felicia Scharitzer.]

  - Marriage of Finite Element Method and Exterior Calculus
  - FEM formulated using Differential Forms
  - FEEC tackles weak variational form of PDE
]

#slide[
  // What was Exterior Calculus again?
  = Exterior Calculus of Differential Forms
  #snote[What was that again?]

  - Modern formulation of traditional vector calculus
  - Calculus on Riemannian Manifolds (Differential Geometry)
  - Unification of all kinds of integrals and derivatives
  - Generalization to arbitrary dimensions
]


#slide[
  = Differential $k$-Form
  #snote[The central object of study.]
  #v(0.5cm)

  Differential $k$-form $omega in Lambda^k (Omega)$ is $k$-dimensional *integrand*!

  #let this(content) = text(fill: blue.lighten(20%), content)
  #let that(content) = text(fill: red, content)

  // basis k-forms in red
  // coefficent function in blue
  // coefficents vary over space
  // compare vector field also with space varying components
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
  #snote[What is it?]

  - $k$-dimensional ruler $omega in Lambda^k (Omega)$
  - ruler continuously varies across space $omega: p in Omega |-> omega_p$ according to #text(blue.lighten(20%))[coefficent functions]
  - locally measures tangential $k$-vectors $omega_p: (T_p M)^k -> RR$
  - globally measures $k$-dimensional submanifold $integral_M omega in RR$

  #v(1cm)
  $
    phi: [0,1]^k -> Omega
    quad quad
    M = "Image" phi
    \
    integral_M omega =
    limits(integral dots.c integral)_([0,1]^k) quad
    omega_(avec(phi)(t))
    ((diff avec(phi))/(diff t_1) wedge dots.c wedge (diff avec(phi))/(diff t_k))
    dif t_1 dots dif t_k
  $

]

//#slide[ 
//  = $k$-Vectors?
//  #snote[What's that now?]
//
//  - Generalization of vectors (1-dimensional oriented line segment)
//  - Oriented $k$-dimensional segments
//
//  #v(0.5cm)
//
//  #set par(spacing: 6pt)
//  #grid(
//    //stroke: 1pt + white,
//    columns: (1fr, 1fr),
//    align: center + horizon,
//    [
//      Bivector / 2-vector
//      #image("res/bivector-in-3d.svg", height: 50%)
//    ], [
//      Trivector / 3-vector
//      #image("res/trivector-in-3d.svg", height: 50%)
//    ]
//  )
//]

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
  = Theorem Unification
  #snote[One theorem to rule them all.]
  
  // Analysis 2
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
      Stokes' Theorem
      $
        integral_D dif omega = integral_(diff D) omega
      $
    ]
  )

]

#slide[
  = How to FEEC?
  #v(1cm)

  Unify and generalize FEM using Differential Forms!
]


#slide[
  = Domain
  #v(1cm)

  - PDE Domain is Riemannian Manifold $Omega$
  - Treating Domains of full topological generality (arbitrary Betti numbers)
  - $k$-th Betti numbers = number of $k$-dim holes

  #grid(
    columns: (50%, 50%), 
    align: center + horizon,
    image("res/embedding.png"),
    image("res/torus.png"),
  )
]

#slide[
  = Mesh
  #v(1cm)

  #grid(
    columns: (60%, 40%), 
    [
      - Discretize PDE Domain
      - Triangulation $mesh$ of manifold $Omega$
      - Simplicial Complex
      - Structure-preserving discretization of continuous manifold
      - Topology, Geometry and Homology
    ], [
      #set align(center + horizon)
      #image("res/moebius.png")
    ]
  )

  
  #grid(
    columns: (50%, 50%), 
    align: center + horizon,
    image("res/simplices.png", width: 100%),
    image("res/simplicialcomplex.png"),
  )

]

#slide[
  = Vector-valued FEM
  #snote[Beyond scalar-valued FEM]

  - NumPDE: Only scalar-valued PDEs $u: Omega -> RR$
  - Only Lagrangian FEM.
  - But there's more!
  - Maxwells Equations!
  - Electric Field $avec(E): Omega -> RR^3$ and Magnetic Field $avec(B): Omega -> RR^3$

  $
    &div avec(E) = rho/epsilon_0
    quad quad
    &&curl E = -(diff avec(B))/(diff t)
    \
    &div avec(B) = 0
    quad quad
    &&curl avec(B) = mu_0 (avec(J) + epsilon_0 (diff avec(E))/(diff t))
  $
]

#slide[
  = Vector-valued Sobolev Spaces
  #snote[It's not so bad, don't worry.]
  
  - We need function spaces for these vector fields
  - Weak formulation: Integrals over $curl$ and $div$

  #set align(center)
  #alternatives(
    position: top,
    $
      H^1(Omega) &= { u: Omega -> RR : integral_Omega norm(grad u)^2 < oo}
    $,
    $
      H (grad; Omega) &= { u: Omega -> RR : integral_Omega norm(grad u)^2 < oo}
    $,
    $
      H (grad; Omega) &= { u: Omega -> RR : integral_Omega norm(grad u)^2 < oo}
      \
      Hvec (curl; Omega) &= { avec(u): Omega -> RR^3 : integral_Omega norm(curl u)^2 < oo}
      \
      Hvec (div ; Omega) &= { avec(u): Omega -> RR^3 : integral_Omega abs(div u)^2 < oo}
    $
  )
]

#slide[
  = Vector-valued FE Spaces
  #v(1cm)

  - Finite dimensional subspaces of infite-dimensional function space.
  - Need to meticulously construct $H(circle.filled.small; Omega)$-conforming FE space for each.

  #grid(
    columns: (50%, 50%),
    align: horizon,
    $
      &H    (grad; Omega) &&supset.eq cal(S)^0_1   (mesh) \
      &Hvec (curl; Omega) &&supset.eq bold(cal(N))   (mesh) \
      &Hvec (div ; Omega) &&supset.eq bold(cal(R T)) (mesh) \
    $,
    [
      - Lagrangian basis on vertices $mesh_0$
      - Nédélec basis on edges $mesh_1$
      - Raviart-Thomas basis on faces $mesh_2$
    ]
  )


]

#slide[
  = Classical FEM vs FEEC
  #v(1cm)

  - These spaces seem very seperate...
  - This is because of Vector Calculus #emoji.face.angry
  - Can we unify them using Exterior Calculus?
  - Can we extend them to more dimensions?
  - Yes with FEEC!
]

#slide[
  = Maxwell's Equations with Differential Forms
  #v(1cm)

  - Write Maxwells equation using Differential Forms.
  - Electric Field is 1-form $E in Lambda^1(Omega)$
  - Magnetic Field is 2-form $B in Lambda^2(Omega)$
  - Current Density is 2-form $J in Lambda^2(Omega)$
  - Electric Charge Density is 3-form $rho in Lambda^3(Omega)$

  $
    &dif E = rho/epsilon_0
    quad quad
    &&dif E = -(diff B)/(diff t)
    \
    &dif B = 0
    quad quad
    &&dif B = mu_0 (J + epsilon_0 (diff E)/(diff t))
  $
]

#slide[
  = Relativistic Electrodynamics
  #snote[Einstein #emoji.hands.shake Maxwell]

  - We can do even better!
  - Make use of Relativistic Electrodynamics!
  - Maxwell's Equations on 4D Spacetime Manifold!
  - Faraday 2-form $F = E wedge dif t + B$
  - Current 3-form $J = rho + J wedge dif t$
  - Only two equations!

  #v(1cm)
  $
    dif F = 0 \
    dif (star F) = J \
  $
]


#slide[
  = Sobolev Space of Differential Forms
  #v(1cm)

  - Weak form only involves one kind of derivative!
  - The exterior derivative!
  - Unification: Only one kind of Sobolev space!
  $
    //H Lambda^k (Omega) = { omega in L^2 Lambda^k (Omega) : dif omega in L^2 Lambda^(k+1) (Omega) }
    H Lambda^k (Omega) = { omega in Lambda^k (Omega) : integral_Omega dif omega < oo }
  $

  $
    &H    (grad; Omega) &&=^~ H Lambda^0 (Omega) \
    &Hvec (curl; Omega) &&=^~ H Lambda^1 (Omega) \
    &Hvec (div ; Omega) &&=^~ H Lambda^2 (Omega) \
  $
]

#slide[
  = Whitney FE Space of Differential Forms
  #v(1cm)
  
  - There is only one type of Sobolev space: $H Lambda^k (Omega)$
  - Can we find a simple $H Lambda^k$-conforming FE-Space?
  - Piecewise-linear over cells $mesh_n$ of triangulation $mesh$?
  - Space of Whitney $k$-forms:
  $
    cal(W) Lambda^k (mesh) = "span" {lambda_sigma : sigma in mesh_k}
  $

  $
    cal(W) Lambda^0 (mesh) &=^~ cal(S)^0_1 (mesh) \
    cal(W) Lambda^1 (mesh) &=^~ bold(cal(N)) (mesh) \
    cal(W) Lambda^2 (mesh) &=^~ bold(cal(R T)) (mesh) \
  $
]

#slide[
  = The Structure behind FEEC
  #v(1cm)

  - There is a rich algebraic structure behind FEEC, called the de Rham complex,
    which connects the different spaces and derivatives.
  - It looks like this in vector calculus:

  #[
    #set align(center)
    #alternatives(
      $
        0 -> C^oo (Omega) limits(->)^grad [C^oo (Omega)]^3 limits(->)^curl [C^oo (Omega)]^3 limits(->)^div C^oo (Omega) -> 0
      $,
      $
        0 -> H(grad; Omega) limits(->)^grad Hvec (curl; Omega) limits(->)^curl Hvec (div; Omega) limits(->)^div L^2(Omega) -> 0
      $
    )
  ]

  #pause
  #pause

  - It's a cochain-complex and is intimitaly connected to the dual
    Chain-Complex, which is here the Simplical Complex (mesh)
  - Connects calculus structure (cohomology) to topology of mesh (homology)
  - Necessary to treat domains of full topology generality.
]

#slide[
  = Discrete Subcomplexes of de Rham complex
  #v(1cm)

  - In order to obtain a good discretization of a PDE, we need to preserve the structure of the continuous problem.
  - We need to preserve this complex!
  - So instead of only finding a discrete subspace of the Sobolev spaces, we want to find a elliptic subcomplex!
  $
    0 -> cal(S)^0_1 (mesh) limits(->)^grad bold(cal(N)) (mesh) limits(->)^curl bold(cal(R T)) (mesh) limits(->)^div cal(S)^(-1)_0 (mesh) -> 0
  $
]
  
#slide[
  = Differential Form Complex
  #v(1cm)

  Exterior Calculus streamlines the de Rham complex.\
  #only("2-")[And it extends the de Rham complex to arbitrary $n$ dimensions.\ ]
  #only("3-")[And there are piecewise polynomial differential forms $cal(P)_r Lambda^k$ for any degree $r$.]
  
  #set align(center)
  #alternatives(
    $
      0 -> H Lambda^0 (Omega) limits(->)^dif H Lambda^1 (Omega) limits(->)^dif H Lambda^2 limits(->)^dif H Lambda^3 (Omega) -> 0
      \
      0 -> cal(W) Lambda^0 (mesh) limits(->)^dif cal(W) Lambda^1 (mesh) limits(->)^dif cal(W) Lambda^2 (mesh) limits(->)^dif cal(W) Lambda^3 (mesh) -> 0
    $,
    $
      0 -> H Lambda^0 (Omega) limits(->)^dif dots.c limits(->)^dif H Lambda^n (Omega) -> 0
      \
      0 -> cal(W) Lambda^0 (mesh) limits(->)^dif dots.c limits(->)^dif cal(W) Lambda^n (mesh) -> 0
    $,
    $
      0 -> H Lambda^0 (Omega)  limits(->)^dif dots.c limits(->)^dif H Lambda^n (Omega) -> 0
      \
      0 -> cal(P)_r Lambda^0 (mesh) limits(->)^dif dots.c limits(->)^dif cal(P)_r Lambda^n (mesh) -> 0
    $
  )
]

#slide[
  = What is FEEC?
  #v(1cm)

  - It is extremly general:
    - Arbitrary dimensions
    - Arbitrary topological manifolds
    - Arbitrary $k$-forms
    - Arbitrary polynomial degree FE solutions
  - It is a theoretical framework for establishing well-posedness of PDE problems, by respecting co-/homology.
  - It is a manual for creating a FEM library of extreme generality.
    - My bachelor's thesis:\ Rust Implementation of Finite Element Exterior Calculus on Coordinate-Free Simplicial Manifolds
]


#slide[
  = Thank you for listening!

  //#set page(background: image("res/bg-vibrant.jpg", width: 100%))
  #set align(center + horizon)

  #block()[
    #set align(center)
    #set par(spacing: 10pt)

    #tiaoma.qrcode("https://github.com/luiswirth/feec-pres",
      options: (
        scale: 4.0,
        fg-color: fgcolor,
        bg-color: bgcolor,
      )
    )
    #weblink("https://github.com/luiswirth/feec-pres", "github:luiswirth/feec-pres")
  ]
]


#slide[
  = Hodge-Laplace Problem
  #snote[Generalization of prototypical Poisson equation]

  $
    Delta u = f
  $

  Now $u$ and $f$ are Differential $k$-forms.
  $
    u in Lambda^k (Omega), f in Lambda^k (Omega)
  $

  And the Laplacian becomes the Hodge-Laplace operator.
  $
    Delta: Lambda^k (Omega) -> Lambda^k (Omega)
    \
    Delta = dif delta + delta dif
  $
]

#slide[
  = Coderivative Operator

  Coderivative operator $delta: Lambda^k (Omega) -> Lambda^(k-1) (Omega)$
  defined such that
  $
    star delta omega = (-1)^k dif star omega
  $
]

#slide[
  = Weak Variational Form
  #v(1cm)

  In order to do FEM, we need to change from the strong PDE form into the weark variational form.

  We need to form the $L^2$-inner product with a test "function" $v in Lambda^k (Omega)$.

  The inner product is defined as
  $
    inner(omega, eta)_(L^2 Lambda^k)
    =
    integral_Omega inner(omega_x, eta_x) "vol"
    =
    integral omega wedge star eta
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
    inner(Delta u, v) = inner(f, v)
    quad quad forall v in H Lambda^k (Omega)
  $

  Or in integral form
  $
    integral_Omega ((dif delta + delta dif) u) wedge star v = integral_Omega f wedge star v
  $
]

#slide[
  = Integration by Parts

  $
    integral_Omega dif omega wedge eta
    =
    (-1)^(k-1)
    integral_Omega omega wedge dif eta
    +
    integral_(diff Omega) "Tr" omega wedge "Tr" eta
  $

  $
    inner(dif omega, eta) = inner(omega, delta eta) + integral_(diff Omega) "Tr" omega wedge "Tr" star eta
  $

  If $omega$ or $eta$ vanishes on the boundary, then
  $delta$ is the formal adjoint of $dif$ w.r.t. the $L^2$-inner product.
  $
    inner(dif omega, eta) = inner(omega, delta eta)
  $

  $
    inner(Delta u, v) = inner(f, v)
    \
    inner((dif delta + delta dif) u, v) = inner(f, v)
    \
    inner((dif delta + delta dif) u, v) = inner(f, v)
    \
    inner(dif delta u, v) + inner(delta dif u, v) = inner(f, v)
    \
    inner(delta u, delta v) + inner(dif u, dif v) = inner(f, v)
  $

  $
    u in H Lambda^k (Omega): quad quad
    inner(delta u, delta v) + inner(dif u, dif v) = inner(f, v)
    quad
    forall v in H Lambda^k (Omega)
  $

  $
    u in H Lambda^k (Omega): quad
    integral_Omega (delta u) wedge star (delta v) + integral_Omega (dif u) wedge star (dif v) = integral_Omega f wedge star v
    quad
    forall v in H Lambda^k (Omega)
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




