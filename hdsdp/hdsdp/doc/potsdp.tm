<TeXmacs|2.1>

<style|<tuple|generic|triangle-list>>

<\body>
  <\hide-preamble>
    \;

    <assign|fa|<macro|\<cal-A\>>>
  </hide-preamble>

  <doc-data|<doc-title|On the design of a dual potential reduction
  solver>|<doc-author|<author-data|<\author-affiliation>
    \;

    HDSDP Developers' Group

    \;

    \;

    <date|>
  </author-affiliation>>>>

  In this note we describe the implementation of a dual potential reduction
  solver that exploits either the embedding or the big-<math|M> potential
  reduction method to solve\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|x>>|<cell|<around*|\<langle\>|c,x|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|\<cal-A\> x=b>|<cell|>>|<row|<cell|>|<cell|x\<succeq\><rsub|\<cal-K\>>0>|<cell|>>>>
  </eqnarray*>

  via its dual

  <\eqnarray*>
    <tformat|<table|<row|<cell|max<rsub|y,s>>|<cell|b<rsup|\<top\>>y>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|\<cal-A\><rsup|\<ast\>>y+s=c>|<cell|>>|<row|<cell|>|<cell|s\<succeq\><rsub|\<cal-K\><rsup|\<ast\>>>0.>|<cell|>>>>
  </eqnarray*>

  The main solver contains the following components

  <\itemize>
    <item>Solver/data

    Interface: set data, set parameter, set solution, optimize, get solution

    Sparse, dense, rank-one

    <item>Algorithm

    HSD/infeasible start/dual potential reduction

    Presolve, phase A, phase B

    Schur complement setup*, potential line-search, barrier line-search

    <item>Linear algebra

    sparse, dense, low rank; eigen, trace, decomposition (Cholesky)

    Lanczos, conjugate gradient, block buffer computation*

    <item>Other utilities

    Parameter tuner, IO and things like that
  </itemize>

  <section|Algorithm>

  <subsection|Notations>

  <math|<around*|\<langle\>|\<cdummy\>,\<cdummy\>|\<rangle\>>> denotes inner
  product; <math|\<succeq\><rsub|\<cal-K\>>> denotes partial order defined
  over cone <math|\<cal-K\>>; the operator
  <math|\<cal-A\>:<big|otimes><rsub|i=1><rsup|k>\<bbb-R\><rsup|n<rsub|i>>\<rightarrow\>\<bbb-R\><rsup|m>>
  and its adjoint <math|\<cal-A\><rsup|\<ast\>>:\<bbb-R\><rsup|m>\<rightarrow\><big|otimes><rsub|i=1><rsup|k>\<bbb-R\><rsup|n<rsub|i>>>
  are respectively defined. They have different expressions for different
  cones. The dual barrier is defined by <math|log det s>. <math|e> is
  generalized to express the unit vector in the conic space.

  <subsection|Basic dual algorithm from a KKT view>\ 

  In the dual method, we are interested in applying the dual algorithm to
  solve the conic problem, whose KKT conditions give

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\>x>|<cell|=>|<cell|b>>|<row|<cell|\<cal-A\><rsup|\<ast\>>y+s>|<cell|=>|<cell|c>>|<row|<cell|x
    s>|<cell|=>|<cell|0<around*|(|\<mu\>e|)>>>|<row|<cell|x,s>|<cell|\<succeq\>>|<cell|0,>>>>
  </eqnarray*>

  where <math|x s=\<mu\>e> denotes the perturbed KKT system often used in the
  analysis of the interior point method.

  From the Newton's perspective towards the condition, given
  <math|<around*|(|x,y,s|)>> such that <math|<around*|(|y,s|)>> satisfies
  <math|\<cal-A\><rsup|\<ast\>>y+s=c>, we wish to find
  <math|\<Delta\>x,\<Delta\>s> such that

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\><around*|(|x+\<Delta\>x|)>>|<cell|=>|<cell|b>>|<row|<cell|\<cal-A\><rsup|\<ast\>>\<Delta\>y+\<Delta\>s>|<cell|=>|<cell|<with|color|red|0>>>|<row|<cell|<around*|(|x+\<Delta\>x|)><around*|(|s+\<Delta\>s|)>>|<cell|=>|<cell|\<mu\>e.>>>>
  </eqnarray*>

  Since we assume a dual feasible solution, the second term has no related
  residual. Dual method modifies the third condition and

  <\eqnarray*>
    <tformat|<table|<row|<cell|x+\<Delta\>x>|<cell|=>|<cell|\<mu\><around*|(|s+\<Delta\>s|)><rsup|-1>.>>>>
  </eqnarray*>

  Linearizing <math|<around*|(|s+\<Delta\>s|)><rsup|-1>\<approx\>s<rsup|-1>-s<rsup|-1>\<Delta\>s
  s<rsup|-1>> gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<mu\>s<rsup|-1>\<Delta\> s
    s<rsup|-1>+\<Delta\>x>|<cell|=>|<cell|\<mu\>s<rsup|-1>-x.>>>>
  </eqnarray*>

  Combined with the above relations, we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\>\<Delta\>x>|<cell|=>|<cell|\<mu\>\<cal-A\>s<rsup|-1>-\<cal-A\>x-\<mu\>\<cal-A\>s<rsup|-1>\<Delta\>
    s s<rsup|-1>>>|<row|<cell|>|<cell|=>|<cell|\<mu\>\<cal-A\>s<rsup|-1>-\<cal-A\>x+\<mu\>\<cal-A\>s<rsup|-1><around*|(|\<cal-A\><rsup|\<ast\>>\<Delta\>y|)>s<rsup|-1>>>|<row|<cell|>|<cell|=>|<cell|\<mu\>\<cal-A\>s<rsup|-1><with|color|red|-\<cal-A\>x>+\<mu\>\<cal-A\>s<rsup|-2>\<cal-A\><rsup|\<ast\>>\<Delta\>y>>|<row|<cell|>|<cell|=>|<cell|b<with|color|red|-\<cal-A\>x>>>>>
  </eqnarray*>

  and we arrive at

  <\equation*>
    \<cal-A\>s<rsup|-2>\<cal-A\><rsup|\<ast\>>\<Delta\>y=<frac|1|\<mu\>>b-\<cal-A\>s<rsup|-1>.
  </equation*>

  Solving the system provides <math|<around*|(|\<Delta\>y,\<Delta\>s=-\<cal-A\><rsup|\<ast\>>\<Delta\>y|)>>
  we need.

  <subsection|Homogeneous dual method>

  Homogeneous dual method works almost the same as dual method but applies
  simplified embedding trick

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\>x-b\<tau\>>|<cell|=>|<cell|0>>|<row|<cell|-\<cal-A\><rsup|\<ast\>>y-s+c\<tau\>>|<cell|=>|<cell|0>>|<row|<cell|<around*|\<langle\>|b,y|\<rangle\>>-<around*|\<langle\>|c,x|\<rangle\>>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|x
    s>|<cell|=>|<cell|\<mu\>e>>|<row|<cell|\<kappa\>\<tau\>>|<cell|=>|<cell|\<mu\>>>>>
  </eqnarray*>

  By treating <math|\<tau\>> as dual variable, we have, similarly

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\><around*|(|x+\<Delta\>x|)>-b<around*|(|\<tau\>+\<Delta\>\<tau\>|)>>|<cell|=>|<cell|0>>|<row|<cell|-\<cal-A\><rsup|\<ast\>><around*|(|y+\<Delta\>y|)>-<around*|(|s+\<Delta\>s|)>+c<around*|(|\<tau\>+\<Delta\>\<tau\>|)>>|<cell|=>|<cell|0>>|<row|<cell|<around*|\<langle\>|b,y+\<Delta\>y|\<rangle\>>-<around*|\<langle\>|c,x+\<Delta\>x|\<rangle\>>-<around*|(|\<kappa\>+\<Delta\>\<kappa\>|)>>|<cell|=>|<cell|0>>|<row|<cell|<around*|(|x+\<Delta\>x|)><around*|(|s+\<Delta\>s|)>>|<cell|=>|<cell|\<mu\>e>>|<row|<cell|<around*|(|\<kappa\>+\<Delta\>\<kappa\>|)><around*|(|\<tau\>+\<Delta\>\<tau\>|)>>|<cell|=>|<cell|\<mu\>.>>>>
  </eqnarray*>

  Since we do not assume a dual feasible solution this time,
  <math|r<rsup|d>\<assign\>-\<cal-A\><rsup|\<ast\>>y-s+c\<tau\>> may be
  nonzero and we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|<fa>\<Delta\>x-b\<Delta\>\<tau\>>|<cell|=>|<cell|b\<tau\>-<fa>x>>|<row|<cell|-<fa><rsup|\<ast\>>\<Delta\>y-\<Delta\>s+c\<Delta\>\<tau\>>|<cell|=>|<cell|\<cal-A\><rsup|\<ast\>>y+s-c\<tau\>>>|<row|<cell|<around*|\<langle\>|b,\<Delta\>y|\<rangle\>>-<around*|\<langle\>|c,\<Delta\>x|\<rangle\>>-\<Delta\>\<kappa\>>|<cell|=>|<cell|-<around*|\<langle\>|b,y|\<rangle\>>+<around*|\<langle\>|c,x|\<rangle\>>+\<kappa\>>>|<row|<cell|\<mu\>s<rsup|-1>\<Delta\>
    s s<rsup|-1>+\<Delta\>x>|<cell|=>|<cell|\<mu\>s<rsup|-1>-x>>|<row|<cell|\<mu\>\<tau\><rsup|-2>\<Delta\>\<tau\>+\<Delta\>\<kappa\>>|<cell|=>|<cell|\<mu\>\<tau\><rsup|-1>-\<kappa\>.>>>>
  </eqnarray*>

  After some tedius simplification we arrive at\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<fa>s<rsup|-2><fa><rsup|*\<ast\>>\<Delta\>y-<around*|(|b+\<mu\><fa>s<rsup|-1>c
    s<rsup|-1>|)>\<Delta\>\<tau\>>|<cell|=>|<cell|b\<tau\>-\<mu\><fa>s<rsup|-1>+\<mu\><fa>s<rsup|-1>r<rsup|d>s<rsup|-1>>>|<row|<cell|<around*|(|b-\<mu\><fa>s<rsup|-1>c
    s<rsup|-1>|)>\<Delta\>y+<around*|(|\<mu\><around*|\<langle\>|c,s<rsup|-1>c
    s<rsup|-1>|\<rangle\>>+\<kappa\>\<tau\><rsup|-1>|)>\<Delta\>\<tau\>>|<cell|=>|<cell|-<around*|\<langle\>|b,y|\<rangle\>>+\<mu\>\<tau\><rsup|-1>+\<mu\><around*|\<langle\>|c,s<rsup|-1>|\<rangle\>>-\<mu\><around*|\<langle\>|c,s<rsup|-1>r<rsup|d>s<rsup|-1>|\<rangle\>>>>>>
  </eqnarray*>

  and we can update <math|<around*|(|\<Delta\>y,\<Delta\>s,\<Delta\>\<tau\>|)>>
  in the same fashion.

  <subsection|Conic operations and mixed cones>

  It is clear that we need <math|<around*|\<langle\>|c,x|\<rangle\>>,<fa>s<rsup|-1>c
  s<rsup|-1>,<around*|\<langle\>|c,s<rsup|-1>c s<rsup|-1>|\<rangle\>>> and
  nutorious <math|\<cal-A\>s<rsup|-2>\<cal-A\><rsup|\<ast\>>> in our
  computation. In different conic contexts we have

  <big-table|<block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<cwith|1|-1|1|-1|cell-valign|c>|<cwith|5|5|1|2|cell-halign|c>|<cwith|5|5|1|2|cell-valign|c>|<cwith|1|-1|1|-1|cell-tborder|0ln>|<cwith|1|-1|1|-1|cell-bborder|0ln>|<cwith|1|-1|1|-1|cell-lborder|1ln>|<cwith|1|-1|1|-1|cell-rborder|1ln>|<cwith|2|-1|1|1|cell-lborder|0ln>|<cwith|2|-1|4|4|cell-rborder|0ln>|<cwith|1|1|1|-1|cell-tborder|1ln>|<cwith|1|1|1|-1|cell-bborder|1ln>|<cwith|1|1|1|-1|cell-lborder|1ln>|<cwith|1|1|1|-1|cell-rborder|1ln>|<cwith|2|2|1|-1|cell-tborder|1ln>|<cwith|1|1|1|1|cell-tborder|1ln>|<cwith|1|1|1|1|cell-bborder|1ln>|<cwith|2|2|1|1|cell-tborder|1ln>|<cwith|1|1|1|1|cell-lborder|0ln>|<cwith|1|1|1|1|cell-rborder|1ln>|<cwith|1|1|2|2|cell-lborder|1ln>|<cwith|10|10|1|-1|cell-tborder|0ln>|<cwith|9|9|1|-1|cell-bborder|0ln>|<cwith|10|10|1|-1|cell-bborder|1ln>|<cwith|10|10|1|1|cell-lborder|0ln>|<cwith|10|10|4|4|cell-rborder|0ln>|<cwith|1|1|4|4|cell-tborder|1ln>|<cwith|1|1|4|4|cell-bborder|1ln>|<cwith|2|2|4|4|cell-tborder|1ln>|<cwith|1|1|4|4|cell-lborder|1ln>|<cwith|1|1|3|3|cell-rborder|1ln>|<cwith|1|1|4|4|cell-rborder|0ln>|<table|<row|<cell|Operation/Cone>|<cell|LP
  <math|<around*|(|x\<geq\>0|)>>>|<cell|SDP
  <math|<around*|(|X\<succeq\>0|)>>>|<cell|<text-dots>>>|<row|<cell|<math|e>>|<cell|<math|<around*|(|1,\<ldots\>,1|)><rsup|\<top\>>>>|<cell|<math|I>>|<cell|>>|<row|<cell|<math|s<rsup|-1>>>|<cell|<math|<around*|(|1/s<rsub|1>,\<ldots\>,1/s<rsub|n>|)>>>|<cell|<math|S<rsup|-1>>>|<cell|>>|<row|<cell|<math|s<rsup|-1>c>>|<cell|<math|<around*|(|c<rsub|1>/s<rsub|1>,\<ldots\>,c<rsub|n>/s<rsub|n>|)>>>|<cell|<math|S<rsup|-1>C>>|<cell|>>|<row|<cell|<math|<around*|\<langle\>|c,x|\<rangle\>>>>|<cell|<math|c<rsup|\<top\>>x>>|<cell|<math|<around*|\<langle\>|C,X|\<rangle\>>>>|<cell|>>|<row|<cell|<math|<fa>x>>|<cell|<math|A
  x>>|<cell|<math|<around*|(|<around*|\<langle\>|A<rsub|1>,X|\<rangle\>>,\<ldots\>,<around*|\<langle\>|A<rsub|m>,X|\<rangle\>>|)>>>|<cell|>>|<row|<cell|<math|<fa><rsup|\<ast\>>y>>|<cell|<math|A<rsup|\<top\>>y>>|<cell|<math|<big|sum><rsub|i>y<rsub|i>A<rsub|i>>>|<cell|>>|<row|<cell|<math|<fa>s<rsup|-1>c
  s<rsup|-1>>>|<cell|<math|A s<rsup|-2>c>>|<cell|<math|<around*|(|<around*|\<langle\>|A<rsub|1>,S<rsup|-1>
  C S<rsup|-1>,\<ldots\>,A<rsub|m>,S<rsup|-1> C
  S<rsup|-1>|\<rangle\>>|)>>>|<cell|>>|<row|<cell|<math|<around*|\<langle\>|c,s<rsup|-1>c
  s<rsup|-1>|\<rangle\>>>>|<cell|<math|c<rsup|\<top\>>s<rsup|-2>c>>|<cell|<math|<around*|\<langle\>|C,S<rsup|-1>
  C S<rsup|-1>|\<rangle\>>>>|<cell|>>|<row|<cell|<math|M=<fa>s<rsup|-2><fa><rsup|\<ast\>>>>|<cell|<math|A
  s<rsup|-2>A<rsup|\<top\>>>>|<cell|<math|M<rsub|i
  j>=<around*|\<langle\>|A<rsub|i>,S<rsup|-1>A<rsub|j>S<rsup|-1>|\<rangle\>>>>|<cell|>>>>>|Conic
  operations>

  When more than one cone is present, by linearity of <math|\<cal-A\>>, we
  simply add the quantities, e.g.,

  <\eqnarray*>
    <tformat|<table|<row|<cell|M>|<cell|=>|<cell|M<rsub|<text|LP>>+M<rsub|<text|SDP>>>>|<row|<cell|<fa>s<rsup|-1>c
    s<rsup|-1>>|<cell|=>|<cell|A s<rsup|-2>c+<around*|(|<around*|\<langle\>|A<rsub|1>,S<rsup|-1>
    C S<rsup|-1>,\<ldots\>,A<rsub|m>,S<rsup|-1> C
    S<rsup|-1>|\<rangle\>>|)><rsup|\<top\>>>>|<row|<cell|<around*|\<langle\>|c,s<rsup|-1>c
    s<rsup|-1>|\<rangle\>>>|<cell|=>|<cell|c<rsup|\<top\>>s<rsup|-2>c+<around*|\<langle\>|C,S<rsup|-1>
    C S<rsup|-1>|\<rangle\>>.>>>>
  </eqnarray*>

  When other cones are added, we only need to define their conic operations
  and add them together.

  <subsection|Primal solution recovery and cone projection with
  infeasibility>

  Recall that in the dual formulation we have
  <math|<around*|(|y,s,\<tau\>|)>> and it is desirable to recover some
  <math|x> satisfying\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\>x-b\<tau\>>|<cell|=>|<cell|0>>|<row|<cell|\<cal-A\><rsup|\<ast\>>y+s-c\<tau\>>|<cell|=>|<cell|0>>|<row|<cell|b<rsup|\<top\>>y-c<rsup|\<top\>>x-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|x
    s>|<cell|=>|<cell|\<mu\>e>>|<row|<cell|\<kappa\>\<tau\>>|<cell|=>|<cell|\<mu\>>>|<row|<cell|<around*|(|x,s,\<kappa\>,\<tau\>|)>>|<cell|\<geq\>>|<cell|0>>>>
  </eqnarray*>

  Substituting <math|x=\<mu\>s<rsup|-1>> gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<mu\>\<cal-A\>s<rsup|-1>-b\<tau\>>|<cell|=>|<cell|0>>|<row|<cell|\<cal-A\><rsup|\<ast\>>y+s-c\<tau\>>|<cell|=>|<cell|0>>|<row|<cell|b<rsup|\<top\>>y-\<mu\>c<rsup|\<top\>>s<rsup|-1>-\<mu\>\<tau\><rsup|-1>>|<cell|=>|<cell|0>>|<row|<cell|<around*|(|s,\<tau\>|)>>|<cell|\<geq\>>|<cell|0>>>>
  </eqnarray*>

  Given <math|y,s,\<tau\>> such that <math|\<cal-A\><rsup|\<ast\>>y+s-c=r e>,
  where <math|r> is the residual of the dual feasibility system the proximity
  norm <math|\<delta\><rsub|\<mu\>><around*|(|y,s,\<tau\>|)>> and
  <math|\<delta\><rsub|\<mu\>><around*|(|y,s|)>> are given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<delta\><rsub|\<mu\>><around*|(|y,s,\<tau\>|)>>|<cell|\<assign\>>|<cell|min<rsub|x,\<kappa\>><frac|1|2><around*|\<\|\|\>|s
    x-\<mu\>e|\<\|\|\>><rsup|2>+<frac|1|2><around*|(|\<kappa\>\<tau\>-\<mu\>|)><rsup|2>,>>|<row|<cell|<text|subject
    to>>|<cell|>|<cell|\<cal-A\> x-b\<tau\>=0>>|<row|<cell|>|<cell|>|<cell|b<rsup|\<top\>>y-c<rsup|\<top\>>x-\<kappa\>=0>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<delta\><rsub|\<mu\>><around*|(|y,s|)>>|<cell|\<assign\>>|<cell|min<rsub|x,\<kappa\>><frac|1|2><around*|\<\|\|\>|s
    x-\<mu\>e|\<\|\|\>><rsup|2>>>|<row|<cell|<text|subject
    to>>|<cell|>|<cell|\<cal-A\> x-b=0,>>>>
  </eqnarray*>

  where we consider the subproblem

  <hrule>

  <\eqnarray*>
    <tformat|<table|<row|<cell|L<around*|(|x,y|)>>|<cell|=>|<cell|<frac|1|2><around*|\<\|\|\>|s
    x-\<mu\>e|\<\|\|\>><rsup|2>-\<mu\>\<Delta\>y<rsup|\<top\>><around*|(|\<cal-A\>
    x-b|)>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|s<rsup|2>x-\<mu\>s-\<mu\>\<cal-A\><rsup|*\<ast\>>\<Delta\>y>|<cell|=>|<cell|0>>|<row|<cell|\<cal-A\>x-b>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  <\equation*>
    x=\<mu\>s<rsup|-1>+\<mu\>s<rsup|-2>\<cal-A\><rsup|*\<ast\>>\<Delta\>y
  </equation*>

  <\equation*>
    \<cal-A\>s<rsup|-2>\<cal-A\><rsup|*\<ast\>>\<Delta\>y=\<mu\><rsup|-1>b-\<cal-A\>s<rsup|-1>\<Rightarrow\>\<Delta\>y=\<mu\><rsup|-1>\<Delta\>y<rsub|1>-\<Delta\>y<rsub|2>
  </equation*>

  \;

  <\eqnarray*>
    <tformat|<table|<row|<cell|x>|<cell|=>|<cell|\<mu\>s<rsup|-1>+\<mu\>s<rsup|-2>\<cal-A\><rsup|*\<ast\>>\<Delta\>y>>|<row|<cell|>|<cell|=>|<cell|\<mu\>s<rsup|-1><around*|(|s-\<cal-A\><rsup|*\<ast\>>\<Delta\>y|)>s<rsup|-1>>>|<row|<cell|>|<cell|=>|<cell|\<mu\>s<rsup|-1><around*|(|c+r
    e-\<cal-A\><rsup|\<ast\>><around*|(|y+\<mu\><rsup|-1>\<Delta\>y<rsub|1>-\<Delta\>y<rsub|2>|)>|)>s<rsup|-1>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|c+r
    e-\<cal-A\><rsup|\<ast\>><around*|(|y+<frac|1|\<mu\>>\<Delta\>y<rsub|1>-\<Delta\>y<rsub|2>|)>>>>>
  </eqnarray*>

  <hrule>

  <\eqnarray*>
    <tformat|<table|<row|<cell|L<around*|(|x,\<kappa\>,v,\<lambda\>|)>>|<cell|=>|<cell|<frac|1|2><around*|\<\|\|\>|S
    x-\<mu\>e|\<\|\|\>><rsup|2>+<frac|1|2><around*|(|\<kappa\>\<tau\>-\<mu\>|)><rsup|2>>>|<row|<cell|>|<cell|>|<cell|+v<rsup|\<top\>><around*|(|\<cal-A\>x-b\<tau\>|)>+\<lambda\><around*|(|b<rsup|\<top\>>y-c<rsup|\<top\>>x-\<kappa\>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|2><around*|\<\|\|\>|S
    x-\<mu\>e|\<\|\|\>><rsup|2>+<frac|1|2><around*|(|\<kappa\>\<tau\>-\<mu\>|)><rsup|2>>>|<row|<cell|>|<cell|>|<cell|+<around*|(|\<cal-A\><rsup|\<ast\>>v-\<lambda\>c|)><rsup|\<top\>>x-b<rsup|\<top\>>v\<tau\>-\<kappa\>\<lambda\>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\>x-b\<tau\>>|<cell|=>|<cell|0>>|<row|<cell|b<rsup|\<top\>>y-c<rsup|\<top\>>x-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|S<around*|(|S
    x-\<mu\>e|)>+\<cal-A\><rsup|\<ast\>>v-\<lambda\>c>|<cell|=>|<cell|0>>|<row|<cell|\<tau\><around*|(|\<kappa\>\<tau\>-\<mu\>|)>-\<lambda\>>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<kappa\>>|<cell|=>|<cell|\<lambda\>\<tau\><rsup|-2>+\<mu\>\<tau\><rsup|-1>>>|<row|<cell|S<rsup|2>x-\<mu\>s+\<cal-A\><rsup|\<ast\>>v-\<lambda\>c>|<cell|=>|<cell|0>>|<row|<cell|x>|<cell|=>|<cell|\<lambda\>S<rsup|-1>c-S<rsup|-1>\<cal-A\><rsup|\<ast\>>v+\<mu\>S<rsup|-1>e>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<cal-A\>x-b\<tau\>>>|<row|<cell|>|<cell|=>|<cell|\<lambda\>\<cal-A\>S<rsup|-1>c-\<cal-A\>S<rsup|-1>\<cal-A\><rsup|\<ast\>>v+\<mu\>\<cal-A\>S<rsup|-1>e-b\<tau\>>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|b<rsup|\<top\>>y-c<rsup|\<top\>>x-\<kappa\>>>|<row|<cell|>|<cell|=>|<cell|b<rsup|\<top\>>y-c<rsup|\<top\>><around*|(|\<lambda\>S<rsup|-1>c-S<rsup|-1>\<cal-A\><rsup|\<ast\>>v+\<mu\>S<rsup|-1>e|)>-\<lambda\>\<tau\><rsup|-2>-\<mu\>\<tau\><rsup|-1>>>|<row|<cell|>|<cell|=>|<cell|b<rsup|\<top\>>y-\<lambda\>c<rsup|\<top\>>S<rsup|-1>c+c<rsup|\<top\>>S<rsup|-1>\<cal-A\><rsup|\<ast\>>v-\<mu\>c<rsup|\<top\>>S<rsup|-1>e-\<lambda\>\<tau\><rsup|-2>-\<mu\>\<tau\><rsup|-1>>>|<row|<cell|>|<cell|=>|<cell|b<rsup|\<top\>>y+<around*|(|\<cal-A\>S<rsup|-1>c|)><rsup|\<top\>>v-\<mu\>c<rsup|\<top\>>S<rsup|-1>e-\<lambda\><around*|(|c<rsup|\<top\>>S<rsup|-1>c+\<tau\><rsup|-2>|)>-\<mu\>\<tau\><rsup|-1>>>|<row|<cell|>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  <hrule>
</body>

<\initial>
  <\collection>
    <associate|font-base-size|9>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|3.2|4>>
    <associate|auto-11|<tuple|3.2.1|5>>
    <associate|auto-12|<tuple|3.2.2|5>>
    <associate|auto-13|<tuple|3.2.3|5>>
    <associate|auto-14|<tuple|3.2.4|5>>
    <associate|auto-15|<tuple|3.3|5>>
    <associate|auto-16|<tuple|3.4|5>>
    <associate|auto-17|<tuple|4|5>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|1.2|1>>
    <associate|auto-4|<tuple|1.3|2>>
    <associate|auto-5|<tuple|1.4|3>>
    <associate|auto-6|<tuple|1|3>>
    <associate|auto-7|<tuple|1.5|3>>
    <associate|auto-8|<tuple|3|4>>
    <associate|auto-9|<tuple|3.1|4>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Conic
      operations>|<pageref|auto-6>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Algorithm>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Notations
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>Basic dual algorithm from a
      KKT view <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|1.3<space|2spc>Homogeneous dual method
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1tab>|1.4<space|2spc>Conic operations and mixed
      cones <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Algorithm
      and conic interface> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>SDP
      data structures> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Factorized data
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>SDP coefficient matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|2tab>|3.2.1<space|2spc>Sparse matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|2tab>|3.2.2<space|2spc>Dense matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|2tab>|3.2.3<space|2spc>Rank-one sparse matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|2tab>|3.2.4<space|2spc>Rank-one dense matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|1tab>|3.3<space|2spc>SDP variable and step
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|1tab>|3.4<space|2spc>Schur complement matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Contribution
      and formats> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>