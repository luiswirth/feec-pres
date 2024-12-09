#import "@preview/touying:0.5.3": *
#import "@preview/fletcher:0.5.2" as fletcher: diagram, node, edge
#import "@preview/tiaoma:0.2.1"


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
  //fill: blue,
  link(..args)
)

#let snote(content) = [
  #set text(20pt)
  #content
]

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


#slide[
  #set page(background: image("res/bg-vibrant.jpg", fit: "cover"))
  #set align(center)

  #v(1cm)
  
  #block(
    fill: black.transparentize(20%),
    outset: 20pt,
    radius: 0pt,
  )[
    
    #[
      #set text(size: 35pt, weight: "bold")

      Finite Element Exterior Calculus
    ]

    #v(1cm)
    #smallcaps[Luis Wirth]\
    ETH Zürich\
    #weblink("luwirth@ethz.ch")\
    #weblink("lwirth.com")
  ]
]

#slide[
  = Finite Element Exterior Calculus
  #v(0.5cm)

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
  - ruler $omega: p in Omega |-> omega_p$ varies continuously  across manifold according to #text(blue.lighten(20%))[coefficent functions]
  - locally measures tangential $k$-vectors $omega_p: (T_p M)^k -> RR$
  - globally measures $k$-dimensional submanifold $integral_M omega in RR$
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
        grad f &= (dif f)^transp
        \
        "curl" vvec(F) &= star d star vvec(F)^transp
        \
        "div" vvec(F) &= (star dif vvec(F)^transp)^transp
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

  $  
    0 limits(<-)^diff C_0 (Omega) limits(<-)^diff C_1 (Omega) limits(<-)^diff C_2 (Omega) limits(<-)^diff C_3 (Omega) limits(<-)^diff 0
  $

]

#slide[
  = Mesh
  #v(1cm)

  #grid(
    columns: (60%, 40%), 
    [
      - Discretize PDE Domain
      - Obtain Simplicial Complex $mesh$ by triangulating manifold $Omega$
      - Contains all $k$-dim simplicies $Delta_k (mesh)$
      - Preserve topology and geometry of continuous manifold
    ], [
      #set align(center + horizon)
      #image("res/moebius.png")
    ]
  )

  #set align(center)
  #set block(below: 1pt)
  #image("res/simplices.png", width: 80%)
  $
    0 limits(<--)^diff #h(1cm) Delta_0 (mesh) #h(1cm) limits(<--)^diff #h(1cm) Delta_1 (mesh) #h(1cm) limits(<--)^diff #h(1cm) Delta_2 (mesh) #h(1cm) limits(<--)^diff #h(1cm) Delta_3 (mesh) #h(1cm) limits(<--)^diff_0
  $
]

#slide[
  = FEM in Vector-Calculus
  #snote[To motivate FEEC]

  NumPDE course: Mostly scalar-valued PDEs $u: Omega -> RR$ \
  But we could also have a vector-valued $avec(u): Omega -> RR^3$.
  #pause

  Maxwells Equations \
  Electric Field $avec(E): Omega -> RR^3$ and Magnetic Field $avec(B): Omega -> RR^3$

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
  #v(1cm)
  
  Need function spaces for vector fields.\
  #only("3-")[Weak formulation: Integrals over $curl$ and $div$]

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
  = The de Rham Complex
  #snote[The underlying differential structure]

  - There is a rich algebraic structure connecting these Sobolev spaces
    through their respective derivatives.

  $
    0 -> H (grad; Omega) limits(->)^grad Hvec (curl; Omega) limits(->)^curl Hvec (div; Omega) limits(->)^div L^2(Omega) -> 0
  $
  #pause
  
  $
    curl compose grad = 0
    quad quad
    div compose curl = 0
  $
]

#slide[
  = Vector-valued FE Spaces
  #v(1cm)

  - Finite dimensional subspaces of infinite-dimensional function space.
  - Need to meticulously construct $H(circle.filled.small; Omega)$-conforming FE space for each.

  #grid(
    columns: (40%, 60%),
    align: horizon,
    $
      &H    (grad; Omega) &&supset.eq cal(S)^0_1   (mesh) \
      &Hvec (curl; Omega) &&supset.eq bold(cal(N))   (mesh) \
      &Hvec (div ; Omega) &&supset.eq bold(cal(R T)) (mesh) \
    $,
    [
      - Lagrangian basis on vertices $Delta_0 (mesh)$
      - Nédélec basis on edges $Delta_1 (mesh)$
      - Raviart-Thomas basis on faces $Delta_2 (mesh)$
    ]
  )
]


#slide[
  = Discrete Subcomplexes of de Rham complex
  #v(1cm)

  - Good discretization of PDE?
  - Preserve the structure of the continuous problem!
  - The de Rham complex!
  $
    0 -> H (grad; Omega) limits(->)^grad Hvec (curl; Omega) limits(->)^curl Hvec (div; Omega) limits(->)^div L^2(Omega) -> 0
  $
  #pause

  Find discrete subcomplex!
  $
    0 -> cal(S)^0_1 (mesh) limits(->)^grad bold(cal(N)) (mesh) limits(->)^curl bold(cal(R T)) (mesh) limits(->)^div cal(S)^(-1)_0 (mesh) -> 0
  $
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

  - Maxwell's Equations on 4D Spacetime Manifold!
  - Faraday 2-form $F = E wedge dif t + B$
  - Current 3-form $J = rho + J wedge dif t$

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
  #pause

  $
    &H    (grad; Omega) &&=^~ H Lambda^0 (Omega) \
    &Hvec (curl; Omega) &&=^~ H Lambda^1 (Omega) \
    &Hvec (div ; Omega) &&=^~ H Lambda^2 (Omega) \
  $
]


#slide[
  = de Rham Complex of Differential Forms
  #v(1cm)

  Exterior Calculus streamlines the de Rham complex.\
  
  $
    0 -> H Lambda^0 (Omega) limits(->)^dif H Lambda^1 (Omega) limits(->)^dif H Lambda^2 limits(->)^dif H Lambda^3 (Omega) -> 0
    \
    dif^2 = dif compose dif = 0
  $
  #pause

  Connected to continuous chain complex\
  $
    0 limits(<-)^diff C_0 (Omega) limits(<-)^diff C_1 (Omega) limits(<-)^diff C_2 (Omega) limits(<-)^diff C_3 (Omega) limits(<-)^diff 0
    \
    diff^2 = diff compose diff = 0
  $
  #pause

  Relating calculus structure (cohomology) to topology of mesh (homology) \
  Necessary to treat domains of full topology generality.
]

#slide[
  = Whitney FE Space of Differential Forms
  #v(1cm)
  
  - Only one type of Sobolev space: $H Lambda^k (Omega)$
    - Only one kind of FE Space!
  - Space of Whitney $k$-forms:
  - Piecewise-linear #text(blue)[coefficents] over cells $Delta_n (mesh)$
  #pause
  
  $
    cal(W) Lambda^k (mesh) = "span" {lambda_sigma : sigma in Delta_k (mesh)}
  $

  #grid(
    columns: (50%, 50%),
    align: center + horizon,
    only("3-")[
      - $cal(W) Lambda^0 (mesh)$ on 0-simplices  $Delta_0 (mesh)$
      - $cal(W) Lambda^1 (mesh)$ on 1-simplicies $Delta_1 (mesh)$
      - $cal(W) Lambda^2 (mesh)$ on 2-simplicies $Delta_2 (mesh)$
    ],
    only("4-")[$
      cal(W) Lambda^0 (mesh) &=^~ cal(S)^0_1 (mesh) \
      cal(W) Lambda^1 (mesh) &=^~ bold(cal(N)) (mesh) \
      cal(W) Lambda^2 (mesh) &=^~ bold(cal(R T)) (mesh) \
    $],
  )
]




#slide[
  = Subcomplex of Differential Forms
  #v(1cm)

  $
    0 -> H Lambda^0 (Omega) limits(->)^dif H Lambda^1 (Omega) limits(->)^dif H Lambda^2 limits(->)^dif H Lambda^3 (Omega) -> 0
  $

  $
    0 -> cal(W) Lambda^0 (mesh) limits(->)^dif cal(W) Lambda^1 (mesh) limits(->)^dif cal(W) Lambda^2 (mesh) limits(->)^dif cal(W) Lambda^3 (mesh) -> 0
  $

  
  $
    0 limits(<-)^diff Delta_0 (mesh) limits(<-)^diff Delta_1 (mesh) limits(<-)^diff Delta_2 (mesh) limits(<-)^diff Delta_3 (mesh) limits(<-)^diff 0
    \
    diff^2 = diff compose diff = 0
  $
]

#slide[
  = Generalization in FEEC
  #v(1cm)

  Extend the de Rham complex to arbitrary $n$ dimensions.\
  #only("2-")[Piecewise polynomial differential forms $cal(P)_r Lambda^k$ for any degree $r$.]

  #set align(center)

  #alternatives(
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
  #only("2-")[
  - A theoretical framework for establishing well-posedness of PDE problems, by respecting co-/homology.
  ]
  #only("3-")[
  - It is a manual for creating a FEM library of extreme generality.
    - My bachelor's thesis:\ Rust Implementation of FEEC on Coordinate-Free Simplicial Manifolds
  ]
]


#slide[
  = Thank you for listening!

  #set align(center + horizon)
  #block()[
    #set align(center)
    #set par(spacing: 10pt)


    #grid(
      columns: (50%, 50%),
      [
      Presentation Slides
      #tiaoma.qrcode("https://github.com/luiswirth/feec-pres",
        options: (
          scale: 4.0,
          fg-color: fgcolor,
          bg-color: bgcolor,
        )
    )
    #weblink(
      "https://github.com/luiswirth/feec-pres",
      "github:luiswirth/feec-pres"
    )
    ], [
      Rust FEEC Library
      #tiaoma.qrcode("https://github.com/luiswirth/formoniq",
        options: (
          scale: 4.0,
          fg-color: fgcolor,
          bg-color: bgcolor,
        )
      )
      #weblink(
        "https://github.com/luiswirth/formoniq",
        "github:luiswirth/formoniq"
      )
    ])
  ]
]

