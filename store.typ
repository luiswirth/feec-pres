#slide[
= Simplicial complexes
#v(10pt)
#image("res/simplices.svg")
\
#text(weight: "bold")[k-simplex] $sigma_k$ = convex hull of $k+1$ points $v_0,v...v_k in RR^n$ ($n >= k$)
\
#only(1)[
$sigma_k={v_0v_1...v_k}$
]
#only(2)[
  #text(fill: red)[$sigma_k=-{v_0v_1..v_k}$]
]
]


#slide[
== Simplex boundary 
#v(20pt)
#only((1,2))[#image(width: 50%,"res/boundary-operator.svg")]
\
\
#only(1)[
#text(size: 25pt)[$diff(sigma_k)=diff{v_0v_1...v_k}=sum_(j=0)^k (-1)^j {v_0,...,hat(v_j),...,v_k}$]
]
#only((2,3))[
#text(size: 25pt)[$diff(sigma_k)=diff{v_0v_1...v_k}=sum_(j=0)^k #text(red)[$(-1)^j$] {v_0,...,hat(v_j),...,v_k}$]
]\
\
#only(3)[$diff{v_0,v_1,v_2} = {v_0,v_1} #text(fill: red)[$-$] {v_1,v_2} + {v_2,v_0}$]
]

#polylux-slide[
  #image("res/cochain-complex.svg")
  #image("res/chain-complex.svg")
]
