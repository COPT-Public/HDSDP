<TeXmacs|2.1>

<style|generic>

<\body>
  <\hide-preamble>
    \;

    <assign|x|<macro|<math|<math-bf|x>>>>

    <assign|A|<macro|<math|<math-bf|A>>>>

    <assign|0|<macro|<math|<math-bf|0>>>>

    <assign|e|<macro|<math|<math-bf|e>>>>

    <assign|n|<macro|<math|\<nabla\>>>>

    <assign|X|<macro|<math|<math-bf|X>>>>

    <assign|i|<macro|<math|<math-bf|I>>>>

    <assign|I|<macro|<math|<math-bf|I>>>>

    <assign|H|<macro|<math|<math-bf|H>>>>

    <assign|h|<macro|<math|<math-bf|h>>>>

    <assign|P|<macro|<math|<math-bf|P>>>>

    <assign|m|<macro|<math|<math-bf|m>>>>

    <assign|g|<macro|<math|<math-bf|g>>>>

    <assign|M|<macro|<math-bf|M>>>

    <assign|d|<macro|<math-bf|d>>>

    <assign|G|<macro|<math-bf|G>>>

    <assign|a|<macro|<math-bf|a>>>

    <assign|b|<macro|<math-bf|b>>>

    <assign|c|<macro|<math-bf|c>>>

    <assign|y|<macro|<math-bf|y>>>

    <assign|s|<macro|<math-bf|s>>>

    <assign|u|<macro|<math|<math-bf|u>>>>

    <assign|HA|<macro|<wide|<math-bf|A>|^>>>

    <assign|bs|<macro|<math-bf|S>>>

    <assign|r|<macro|<math-bf|r>>>

    <assign|Y|<macro|<math-bf|Y>>>

    <assign|p|<macro|<math-bf|p>>>

    <assign|V|<macro|<math-bf|V>>>

    <assign|v|<macro|<math-bf|v>>>

    <assign|D|<macro|<math-bf|D>>>

    <assign|w|<macro|<math-bf|w>>>

    <assign|B|<macro|<math|<math-bf|B>>>>
  </hide-preamble>

  <doc-data|<doc-title|On the Practical Implementataion of a First-order
  Potential Reduction Algorithm for Linear
  Programming>|<doc-author|<author-data|<\author-affiliation>
    Huikang Liu<space|2em> Wenzhi Gao<space|2em>Yinyu Ye

    \;

    \;

    <date|>
  </author-affiliation>>>>

  <abstract-data|<abstract|In this report, we present the detailed
  implementation details for a first-order potential reduction algorithm for
  linear programming (LP) problems. The algorithm applies dimension-reduced
  method to reduce the potential function defined over the well-known
  homogeneous self-dual model for LP and leverages the negative curvature of
  the potential function to accelerate convergence. A complete recipe on
  algorithm design and implementation is depicted in this report and some
  preliminary experiment results are given.>>

  \;

  <section|Introduction>

  <subsection|First-order Potential Reduction Method for LP>

  In this report, we are interested in a first-order method for the standard
  LP problems.

  <strong|Standard Primal-dual LP>

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<x>\<in\>\<bbb-R\><rsup|n>>>|<cell|<c><rsup|\<top\>><x>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<A><x>=<b>>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|max<rsub|<y>\<in\>\<bbb-R\><rsup|m>>>|<cell|<b><rsup|\<top\>><y>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<A><rsup|\<top\>><y>+<s>=<c>>|<cell|>>|<row|<cell|>|<cell|<s>\<geq\><0>>|<cell|>>>>
  </eqnarray*>

  It is well-known that LPs admit a homogeneous self-dual (HSD) model

  <\eqnarray*>
    <tformat|<table|<row|<cell|<A><x>-<b>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|-<A><rsup|\<top\>><y>-<s>+<c>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|<b><rsup|\<top\>><y>-<c><rsup|\<top\>><x>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|<around*|(|<x>,<s>,\<kappa\>,\<tau\>|)>>|<cell|\<geq\>>|<cell|<0>,>>>>
  </eqnarray*>

  where the two homogenizing variables <math|\<kappa\>,\<tau\>> are
  introduced for infeasibility detection <cite|ye2011interior>. The
  first-order potential reduction method, initially proposed by
  <cite|ye2015first>, encodes the above HSD model into the following
  simplex-constrained QP

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<around*|(|<x>,<y>,<s>,\<kappa\>,\<tau\>|)>>>|<cell|<frac|1|2><around*|\<\|\|\>|<r><around*|(|<x>,<y>,<s>,\<kappa\>,\<tau\>|)>|\<\|\|\>><rsup|2>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsub|n><rsup|\<top\>><x>+<e><rsub|n><rsup|\<top\>><s>+\<kappa\>+\<tau\>=1,>|<cell|>>>>
  </eqnarray*>

  where\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<r><around*|(|<x>,<y>,<s>,\<kappa\>,\<tau\>|)>>|<cell|\<assign\>>|<cell|<matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|<A>>|<cell|<0><rsub|m\<times\>n>>|<cell|<0><rsub|m\<times\>1>>|<cell|-<b>>>|<row|<cell|-<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|<c>>>|<row|<cell|<b><rsup|\<top\>>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>|<cell|0>>>>><matrix|<tformat|<table|<row|<cell|<y>>>|<row|<cell|<x>>>|<row|<cell|<s>>>|<row|<cell|\<kappa\>>>|<row|<cell|\<tau\>>>>>>>>>>
  </eqnarray*>

  and for brevity, we simplify the notation by re-defining <math|<A>> and
  <math|<x>> and consider the formulation below

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<x>>>|<cell|<frac|1|2><around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>|<cell|\<backassign\>f<around*|(|<x>|)>>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsup|\<top\>><x>=1>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>.>|<cell|>>>>
  </eqnarray*>

  Given the re-formulation, the first-order potential reduction method adopts
  the potential function

  <\equation*>
    \<phi\><around*|(|<x>|)>\<assign\>\<rho\> log
    <around*|(|f<around*|(|<x>|)>|)>-<big|sum><rsub|i=1><rsup|n>log x<rsub|i>
  </equation*>

  and applies a conditional gradient method to drive <math|\<phi\>> to
  <math|-\<infty\>>. More detailedly, the gradient of <math|\<phi\>> is given
  by

  <\equation*>
    <n>\<phi\><around*|(|<x>|)>=<frac|\<rho\><n>f<around*|(|<x>|)>|f<around*|(|<x>|)>>-<X><rsup|-1><e>.
  </equation*>

  At each iteration, we evaluate the gradient
  <math|<n>\<phi\><around*|(|<x><rsup|k>|)>>, let
  <math|\<Delta\><rsup|k>\<assign\><x><rsup|k+1>-<x><rsup|k>> and solve
  following subproblem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|\<Delta\>>>|<cell|<around*|\<langle\>|<n>\<phi\><around*|(|<x><rsup|k>|)>,\<Delta\>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsup|\<top\>>\<Delta\><rsup|k>=0>|<cell|>>|<row|<cell|>|<cell|<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1>\<Delta\><rsup|k>|\<\|\|\>>\<leq\>\<beta\>>|<cell|>>>>
  </eqnarray*>

  to update <math|<x><rsup|k+1>\<leftarrow\><x><rsup|k>+\<Delta\><rsup|k>>.
  In the next section, we extend the basic potential reduction framwork by
  incorporating momentum term from the Dimension-reduced method proposed in
  <cite|zhang2022drsom>.

  <subsection|Dimension-reduced Potential Reduction>

  In this section, we consider two direction extension of the potential
  reduction framework. In a word, by keeping track of one recent history
  iterate, we update

  <\eqnarray*>
    <tformat|<table|<row|<cell|<d><rsup|k>>|<cell|\<leftarrow\>>|<cell|\<alpha\><rsup|g><P><rsub|\<Delta\>><around*|[|\<nabla\>\<phi\><around*|(|<x><rsup|k>|)>|]>+\<alpha\><rsup|m><around*|(|<x><rsup|k>-<x><rsup|k-1>|)>>>|<row|<cell|<x><rsup|k>>|<cell|\<leftarrow\>>|<cell|<x><rsup|k>+<d><rsup|k>>>>>
  </eqnarray*>

  where <math|<P><rsub|\<Delta\>><around*|[|\<cdummy\>|]>> is the orthogonal
  projection onto null space of the simplex constraint
  <math|<e><rsup|\<top\>><x>=0>. Since we leverage the dimension-reduced
  method, \ <math|\<alpha\><rsup|g>,\<alpha\><rsup|d>> are evaluated through
  the following model

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<d>,\<alpha\><rsup|g>,\<alpha\><rsup|m>>>|<cell|<frac|1|2><d><rsup|\<top\>>\<nabla\><rsup|2>\<phi\><around*|(|<x>|)><d>+\<nabla\>\<phi\><around*|(|<x>|)><rsup|\<top\>><d>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<around*|\<\|\|\>|<X><rsup|-1><d>|\<\|\|\>>\<leq\>\<Delta\>>|<cell|>>|<row|<cell|>|<cell|<d>=\<alpha\><rsup|g><g><rsup|k>+\<alpha\><rsup|m><m><rsup|k>>|<cell|>>>>
  </eqnarray*>

  where <math|<g><rsup|k>\<assign\><P><rsub|\<Delta\>><around*|[|\<nabla\>\<phi\><around*|(|<x><rsup|k>|)>|]>>,
  <math|<m><rsup|k>\<assign\><x><rsup|k>-<x><rsup|k-1>>. If we define\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|<H>|~>>|<cell|\<assign\>>|<cell|<matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<g><rsup|k>,\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<langle\>|<g><rsup|k>,\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>|<row|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>>>>>>|<row|<cell|<wide|<h>|~>>|<cell|\<assign\>>|<cell|<matrix|<tformat|<table|<row|<cell|<around*|\<\|\|\>|<g><rsup|k>|\<\|\|\>><rsup|2>>>|<row|<cell|<around*|\<langle\>|<g><rsup|k>,<m><rsup|k>|\<rangle\>>>>>>>>>|<row|<cell|<M>>|<cell|\<assign\>>|<cell|<matrix|<tformat|<table|<row|<cell|<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><g><rsup|k>|\<\|\|\>><rsup|2>>|<cell|<around*|\<langle\>|<g><rsup|k>,<around*|(|<X><rsup|k>|)><rsup|-2><m><rsup|k>|\<rangle\>>>>|<row|<cell|<around*|\<langle\>|<m><rsup|k>,<around*|(|<X><rsup|k>|)><rsup|-2><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><m><rsup|k>|\<\|\|\>><rsup|2>>>>>>,>>>>
  </eqnarray*>

  the above model simplies into a two-dimensional QCQP.

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|\<alpha\>>>|<cell|<frac|1|2>\<alpha\><rsup|\<top\>><wide|<H>|~>\<alpha\>+<wide|<h>|~>\<alpha\>>|<cell|\<backassign\>m<around*|(|\<alpha\>|)>>>|<row|<cell|<text|subject
    to>>|<cell|<around*|\<\|\|\>|<M>\<alpha\>|\<\|\|\>>\<leq\>\<Delta\>>|<cell|>>>>
  </eqnarray*>

  \;

  Note that <math|\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)>=-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>>|f<around*|(|<x><rsup|k>|)><rsup|2>>+\<rho\><frac|<A><rsup|\<top\>><A>|f<around*|(|<x><rsup|k>|)>>+<around*|(|<X><rsup|k>|)><rsup|-2>>
  and we evaluate the above relations via

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|\<langle\>|<a>,\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><a>|\<rangle\>>>|<cell|=>|<cell|<around*|\<langle\>|<a>,-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><a>|f<around*|(|<x><rsup|k>|)><rsup|2>>|\<rangle\>>+<frac|<around*|\<\|\|\>|<A><a>|\<\|\|\>><rsup|2>|f<around*|(|<x><rsup|k>|)>>+<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><a>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|-\<rho\><around*|(|<frac|<n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><a>|f<around*|(|<x><rsup|k>|)>>|)><rsup|2>+<frac|<around*|\<\|\|\>|<A><a>|\<\|\|\>><rsup|2>|f<around*|(|<x><rsup|k>|)>>+<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><a>|\<\|\|\>><rsup|2>>>|<row|<cell|<around*|\<langle\>|<a>,\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><b>|\<rangle\>>>|<cell|=>|<cell|<around*|\<langle\>|<a>,-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><b>|f<around*|(|<x><rsup|k>|)><rsup|2>>|\<rangle\>>+<frac|<around*|\<langle\>|<A><a>,<A><b>|\<rangle\>>|f<around*|(|<x><rsup|k>|)>>+<around*|\<langle\>|<a>,<around*|(|<X><rsup|k>|)><rsup|-2><b>|\<rangle\>>>>|<row|<cell|>|<cell|=>|<cell|-\<rho\><around*|(|<frac|<n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><a>|f<around*|(|<x><rsup|k>|)>>|)><around*|(|<frac|<n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><b>|f<around*|(|<x><rsup|k>|)>>|)>+<frac|<around*|\<langle\>|<A><a>,<A><b>|\<rangle\>>|f<around*|(|<x><rsup|k>|)>>+<around*|\<langle\>|<a>,<around*|(|<X><rsup|k>|)><rsup|-2><b>|\<rangle\>>.>>>>
  </eqnarray*>

  To ensure feasibility, we always choose <math|\<Delta\>\<leq\>1> and adjust
  it based on the trust-region rule.

  <section|Accelerating the Dimension-reduced Potential Reduction>

  In this section, we summarize several techniques applied to improve the
  potential reduction method.

  <subsection|Scaling>

  As is often observed in the first-order type methods, proper scaling
  accelerates the performance of the algorithm. In practice, we scale\ 

  <\equation*>
    <b>\<leftarrow\><frac|<b>|<around*|\<\|\|\>|<b>|\<\|\|\>><rsub|1>+1><space|3em><c>\<leftarrow\><frac|<c>|<around*|\<\|\|\>|<c>|\<\|\|\>><rsub|1>+1>
  </equation*>

  and then apply Ruiz scaling <cite|ruiz2001scaling> to
  <math|<matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|<A>>|<cell|<0><rsub|m\<times\>n>>|<cell|<0><rsub|m\<times\>1>>|<cell|-<b>>>|<row|<cell|-<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|<c>>>|<row|<cell|<b><rsup|\<top\>>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>|<cell|0>>>>>>
  to improve conditioning of the matrix.

  <subsection|Line-search>

  When a direction is assembled from the trust-region subproblem, instead of
  directly updating

  <\equation*>
    <x><rsup|k+1>\<leftarrow\><x><rsup|k>+<d><rsup|k>,
  </equation*>

  we allow a more aggressive exploitation of the direction by performing a
  line-search over\ 

  <\equation*>
    \<phi\><around*|(|<x>+\<alpha\><d>|)>,\<alpha\>\<in\><around*|[|0,0.9\<alpha\><rsub|<text|max>>|)>,
  </equation*>

  where <math|\<alpha\><rsub|<text|max>>=max<around*|{|\<alpha\>\<geq\>0,<x>+\<alpha\><d>\<geq\><0>|}>>.
  The line-search sometimes help accelerate convergence when the algorithm
  approaches optimality.

  <subsection|Escaping the Local Optimum>

  One most important accleration trick is to introduce the negative curvature
  as a search direction. Since potential function is nonconvex in nature,
  it's quite common that the algorithm stagates at a local solution. To help
  escape such local optimum, we make use of the negative curvature of
  <math|<n><rsup|2>\<phi\><around*|(|<x>|)>>. In our case this can be done by
  finding the (minimal) negative eigenvalue and eigenvector

  <\equation*>
    \<lambda\><rsub|min><around*|{|\<nabla\><rsup|2>\<phi\><around*|(|<x>|)>=<frac|2\<rho\><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|4>>+<small|<X><rsup|-2>>|}>
  </equation*>

  and we wish to solve the eigen-problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<around*|\<\|\|\>|<v>|\<\|\|\>>=1>>|<cell|<v><rsup|\<top\>><around*|{|<frac|2\<rho\><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|4>>+<small|<X><rsup|-2>>|}><v>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsup|\<top\>><v>=0.>|<cell|>>>>
  </eqnarray*>

  In general there are two ways to compute a valid direction. The first
  method approaches the problem directly and uses Lanczos iteration to find
  the negative eigen-value of <math|<n><rsup|2>\<phi\>>. As for the second
  approach, we apply the scaling matrix <math|<X>> and solve

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<around*|\<\|\|\>|<X><v>|\<\|\|\>>=1>>|<cell|<v><rsup|\<top\>><around*|{|<frac|2\<rho\><X><A><rsup|\<top\>><A><X>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><X><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A><X>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4>>+<small|<I>>|}><small|><v>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<x><rsup|\<top\>><v>=0.>|<cell|>>>>
  </eqnarray*>

  Since we are to find any negative curvature, it is safe to replace
  <math|<around*|\<\|\|\>|<X><v>|\<\|\|\>>=1> by
  <math|<around*|\<\|\|\>|<v>|\<\|\|\>>=1> and arrive at

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<around*|\<\|\|\>|<v>|\<\|\|\>>=1>>|<cell|<v><rsup|\<top\>><around*|{|<frac|2\<rho\><X><A><rsup|\<top\>><A><X>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><X><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A><X>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4>>+<small|<I>>|}><small|><v>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<x><rsup|\<top\>><v>=0.>|<cell|>>>>
  </eqnarray*>

  Another useful technique when evaluating the curvature is to reduce the
  support of the curvature. Since it's likely that <math|v<rsub|j>>,
  <math|j\<in\><around*|{|i:x<rsub|i>\<rightarrow\>0|}>> will contribute a
  lot in the negative curvature, we can restrict the support of <math|<v>> to
  <math|<around*|{|i:x<rsub|i>\<geq\>\<varepsilon\>|}>> for some
  <math|\<varepsilon\>\<gtr\>0>.

  <section|Algorithm Design>

  In this section, we discuss the design of the potential-reduction based
  solver.\ 

  <subsection|Abstract Function Class>

  To allow further extension, we design the solver to solve general problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<x>>>|<cell|f<around*|(|<x>|)>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<A><x>=<b>>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>>|<cell|>>>>
  </eqnarray*>

  using potential reduction

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<phi\><around*|(|<x>|)>>|<cell|\<assign\>>|<cell|\<rho\>log<around*|(|f<around*|(|<x>|)>-z|)>+log
    <big|sum><rsub|i=1><rsup|n>x<rsub|i>>>>>
  </eqnarray*>

  over Null space of <math|<A>>. To drive the method to work, the following
  methods should be provided by <math|f>.

  <\itemize>
    <item>Gradient evaluation <math|\<nabla\>f<around*|(|<x>|)>>

    <item>Hessian vector product <math|<n><rsup|2>f<around*|(|<x>|)><u>>

    <item>Progress monitor (Optional)
  </itemize>

  The potential reduction framework requries <math|<A>> and maintains
  <math|<x>,z,\<rho\>> to run\ 

  <\itemize>
    <item>Potential gradient evaluation\ 

    <\equation*>
      <n>\<phi\><around*|(|<x>|)>=<frac|\<rho\>|f<around*|(|<x>|)>-z><n>f<around*|(|<x>|)>-<X><rsup|-1><e>
    </equation*>

    <item>Potential (scaled) Hessian-vector product

    <\eqnarray*>
      <tformat|<table|<row|<cell|<n><rsup|2>\<phi\><around*|(|<x>|)>>|<cell|=>|<cell|-<frac|\<rho\><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>>|<around*|(|f<around*|(|<x>|)>-z|)><rsup|2>>+<frac|<n><rsup|2>f<around*|(|<x>|)>|f<around*|(|<x>|)>-z>+<X><rsup|-2>>>|<row|<cell|<X><n><rsup|2>\<phi\><around*|(|<x>|)><X>>|<cell|=>|<cell|-<frac|\<rho\><X><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>><X><u>|<around*|(|f<around*|(|<x>|)>-z|)><rsup|2>>+<frac|<X><n><rsup|2>f<around*|(|<x>|)><X><u>|f<around*|(|<x>|)>-z>+<u>>>>>
    </eqnarray*>

    <item>(Scaled) projection onto Null space

    <\equation*>
      <around*|(|<I>-<A><rsup|\<top\>><around*|(|<A><A><rsup|\<top\>>|)><rsup|-1><A>|)><u>
    </equation*>

    <\equation*>
      <around*|(|<I>-<X><A><rsup|\<top\>><around*|(|<A><X><rsup|2><A><rsup|\<top\>>|)><rsup|-1><A><X>|)><u>
    </equation*>

    <item>Scaled projected gradient and negative curvature\ 

    <item>Trust-region subproblem

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|\<alpha\>>>|<cell|<frac|1|2>\<alpha\><rsup|\<top\>><H>\<alpha\>+<g><rsup|\<top\>>\<alpha\>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|\<alpha\><rsup|\<top\>><G>\<alpha\>\<leq\>\<beta\>>|<cell|>>>>
    </eqnarray*>

    <item>Heuristic routines

    Line search, Curvature frequency, lower bound update
  </itemize>

  <subsection|Numerical Operations>

  In this section, we introduce how to implement the numerical operations
  from the potential reduction method. Here we define
  <math|<wide|<A>|~>\<assign\><matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|<A>>|<cell|<0><rsub|m\<times\>n>>|<cell|<0><rsub|m\<times\>1>>|<cell|-<b>>>|<row|<cell|-<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|<c>>>|<row|<cell|<b><rsup|\<top\>>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>|<cell|0>>>>>>.

  <strong|Residual setup>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<r><rsub|1>>|<cell|=>|<cell|<A><x>-<b>\<tau\>>>|<row|<cell|<r><rsub|2>>|<cell|=>|<cell|-<A><rsup|\<top\>><y>-<s>+<c>\<tau\>>>|<row|<cell|r<rsub|3>>|<cell|=>|<cell|<b><rsup|\<top\>><y>-<c><rsup|\<top\>><x>-\<kappa\>.>>>>
  </eqnarray*>

  <strong|Objective value>

  <\equation*>
    f=<frac|1|2><around*|[|<around*|\<\|\|\>|<r><rsub|1>|\<\|\|\>><rsup|2>+<around*|\<\|\|\>|<r><rsub|2>|\<\|\|\>><rsup|2>+r<rsub|3><rsup|2>|]>
  </equation*>

  <strong|Gradient setup>

  <\equation*>
    \<nabla\>f=<matrix|<tformat|<table|<row|<cell|-<A><r><rsub|2>+<b>r<rsub|3>>>|<row|<cell|<A><rsup|\<top\>><r><rsub|1>-<c>r<rsub|3>>>|<row|<cell|-<r><rsub|2>>>|<row|<cell|-r<rsub|3>>>|<row|<cell|-<b><rsup|\<top\>><r><rsub|1>+<c><rsup|\<top\>><r><rsub|2>>>>>>
  </equation*>

  <\equation*>
    \<nabla\>\<varphi\>=<frac|\<rho\><n>f|f>-<matrix|<tformat|<table|<row|<cell|<X><rsup|-1><e>>>|<row|<cell|<0><rsub|m>>>|<row|<cell|<bs><rsup|-1><e>>>|<row|<cell|\<kappa\><rsup|-1>>>|<row|<cell|\<tau\><rsup|-1>>>>>>
  </equation*>

  <strong|Hessian-vector (with projection)>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<u>>|<cell|=>|<cell|<x>-<frac|<e><rsup|\<top\>><x>|n>\<cdummy\><e>>>|<row|<cell|<n><rsup|2>\<phi\><u>>|<cell|=>|<cell|-<frac|\<rho\><around*|(|<n>f<rsup|\<top\>><u>|)>|f<rsup|2>><n>f+<frac|\<rho\>|f><wide|<A>|~><rsup|\<top\>><around*|(|<wide|<A>|~><u>|)>+<small|<matrix|<tformat|<table|<row|<cell|<X><rsup|-2>>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|<0><rsub|m\<times\>m>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|<bs><rsup|-2>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|\<kappa\><rsup|-2>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|\<tau\><rsup|-2>>>>>>><u>.>>>>
  </eqnarray*>

  \;

  <strong|Lanczos Hessian-vector (with projection)>

  <\equation*>
    <M>\<assign\><small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><x><rsup|\<top\>>/<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>>>>>><around*|[|<frac|2\<rho\><bs><wide|<A>|~><rsup|\<top\>><wide|<A>|~><bs>|<around*|\<\|\|\>|<wide|<A>|~><u>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><bs><wide|<A>|~><rsup|\<top\>><wide|<A>|~><u><u><rsup|\<top\>><wide|<A>|~><rsup|\<top\>><wide|<A>|~><bs>|<around*|\<\|\|\>|<wide|<A>|~><u>|\<\|\|\>><rsup|4>>+<small|<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>>>>>>>|]><small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><x><rsup|\<top\>>/<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>>>>>>
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<x><rprime|'>>|<cell|\<leftarrow\>>|<cell|<frac|<x>|<around*|\<\|\|\>|<x>|\<\|\|\>>>>>|<row|<cell|<v>>|<cell|\<leftarrow\>>|<cell|<small|<matrix|<tformat|<table|<row|<cell|<v><rsub|<y>>>>|<row|<cell|<v><rsub|<x>>-<around*|(|<x><rprime|'><rsup|\<top\>><v><rsub|<x>>|)><x><rprime|'>>>>>>>>>|<row|<cell|<u><rsub|1>>|<cell|\<leftarrow\>>|<cell|<small|<matrix|<tformat|<table|<row|<cell|<0>>>|<row|<cell|<v><rsub|<x>>-<around*|(|<x><rprime|'><rsup|\<top\>><v><rsub|<x>>|)><x><rprime|'>>>>>>>>>|<row|<cell|<u><rsub|2>>|<cell|\<leftarrow\>>|<cell|<small|<small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><rprime|'><x><rprime|'><rsup|\<top\>>>>>>>><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><wide|<A>|~><rsup|\<top\>><wide|<A>|~><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><v>>>>|<row|<cell|<u><rsub|3>>|<cell|\<leftarrow\>>|<cell|<small|<g><rsup|\<top\>><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><v><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><rprime|'><x><rprime|'><rsup|\<top\>>>>>>><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><g>>>>|<row|<cell|<M><v>>|<cell|\<leftarrow\>>|<cell|<frac|f<rsup|2><u><rsub|1>+\<rho\>f<u><rsub|2>-\<rho\><u><rsub|3>|<around*|\<\|\|\>|<wide|<A>|~><u>|\<\|\|\>><rsup|4>>>>>>
  </eqnarray*>

  \;

  <section|Numerical Experiments>

  We provide some preliminary computational results on the NETLIB LP
  problems. The results are obtained using <verbatim|MATLAB> after 1000
  iterations.

  <big-table|<small|<block|<tformat|<cwith|2|13|1|4|cell-valign|c>|<cwith|1|-1|1|-1|cell-halign|c>|<cwith|1|13|5|8|cell-halign|c>|<cwith|2|7|6|8|cell-valign|c>|<cwith|2|7|5|5|cell-valign|c>|<cwith|14|18|1|4|cell-halign|c>|<cwith|14|18|1|4|cell-halign|c>|<cwith|14|18|2|4|cell-valign|c>|<cwith|14|18|1|1|cell-valign|c>|<cwith|2|3|5|8|cell-halign|c>|<cwith|2|3|5|8|cell-halign|c>|<cwith|2|2|6|8|cell-valign|c>|<cwith|2|2|5|5|cell-valign|c>|<cwith|4|13|5|5|cell-halign|c>|<cwith|4|8|5|5|cell-halign|c>|<cwith|4|13|5|5|cell-valign|c>|<cwith|4|13|6|8|cell-valign|c>|<cwith|14|21|5|5|cell-rborder|0ln>|<cwith|14|21|5|5|cell-valign|c>|<cwith|14|18|6|8|cell-valign|c>|<cwith|19|20|1|4|cell-valign|c>|<cwith|19|19|5|8|cell-halign|c>|<cwith|19|19|5|8|cell-valign|c>|<cwith|20|20|1|4|cell-valign|c>|<cwith|20|20|5|8|cell-valign|c>|<cwith|1|-1|1|-1|cell-tborder|1ln>|<cwith|1|-1|1|-1|cell-bborder|1ln>|<cwith|1|-1|1|-1|cell-lborder|0ln>|<cwith|1|-1|1|-1|cell-rborder|0ln>|<cwith|1|-1|5|5|cell-tborder|1ln>|<cwith|1|-1|5|5|cell-bborder|1ln>|<cwith|1|-1|5|5|cell-lborder|1ln>|<cwith|1|-1|5|5|cell-rborder|1ln>|<cwith|1|-1|4|4|cell-rborder|1ln>|<cwith|1|-1|6|6|cell-lborder|1ln>|<table|<row|<cell|Problem>|<cell|PInfeas>|<cell|DInfeas.>|<cell|Compl.>|<cell|Problem>|<cell|PInfeas>|<cell|DInfeas.>|<cell|Compl.>>|<row|<cell|DLITTLE>|<cell|1.347e-10>|<cell|2.308e-10>|<cell|2.960e-09>|<cell|KB2>|<cell|5.455e-11>|<cell|6.417e-10>|<cell|7.562e-11>>|<row|<cell|AFIRO>|<cell|7.641e-11>|<cell|7.375e-11>|<cell|3.130e-10>|<cell|LOTFI>|<cell|2.164e-09>|<cell|4.155e-09>|<cell|8.663e-08>>|<row|<cell|AGG2>|<cell|3.374e-08>|<cell|4.859e-08>|<cell|6.286e-07>|<cell|MODSZK1>|<cell|1.527e-06>|<cell|5.415e-05>|<cell|2.597e-04>>|<row|<cell|AGG3>|<cell|2.248e-05>|<cell|1.151e-06>|<cell|1.518e-05>|<cell|RECIPELP>|<cell|5.868e-08>|<cell|6.300e-08>|<cell|1.285e-07>>|<row|<cell|BANDM>|<cell|2.444e-09>|<cell|4.886e-09>|<cell|3.769e-08>|<cell|SC105>|<cell|7.315e-11>|<cell|5.970e-11>|<cell|2.435e-10>>|<row|<cell|BEACONFD>|<cell|5.765e-12>|<cell|9.853e-12>|<cell|1.022e-10>|<cell|SC205>|<cell|6.392e-11>|<cell|5.710e-11>|<cell|2.650e-10>>|<row|<cell|BLEND>|<cell|2.018e-10>|<cell|3.729e-10>|<cell|1.179e-09>|<cell|SC50A>|<cell|1.078e-05>|<cell|6.098e-06>|<cell|4.279e-05>>|<row|<cell|BOEING2>|<cell|1.144e-07>|<cell|1.110e-08>|<cell|2.307e-07>|<cell|SC50B>|<cell|4.647e-11>|<cell|3.269e-11>|<cell|1.747e-10>>|<row|<cell|BORE3D>|<cell|2.389e-08>|<cell|5.013e-08>|<cell|1.165e-07>|<cell|SCAGR25>|<cell|1.048e-07>|<cell|5.298e-08>|<cell|1.289e-06>>|<row|<cell|BRANDY>|<cell|2.702e-05>|<cell|7.818e-06>|<cell|1.849e-05>|<cell|SCAGR7>|<cell|1.087e-07>|<cell|1.173e-08>|<cell|2.601e-07>>|<row|<cell|CAPRI>|<cell|7.575e-05>|<cell|4.488e-05>|<cell|4.880e-05>|<cell|SCFXM1>|<cell|4.323e-06>|<cell|5.244e-06>|<cell|8.681e-06>>|<row|<cell|E226>|<cell|2.656e-06>|<cell|4.742e-06>|<cell|2.512e-05>|<cell|SCORPION>|<cell|1.674e-09>|<cell|1.892e-09>|<cell|1.737e-08>>|<row|<cell|FINNIS>|<cell|8.577e-07>|<cell|8.367e-07>|<cell|1.001e-05>|<cell|SCTAP1>|<cell|5.567e-07>|<cell|8.430e-07>|<cell|5.081e-06>>|<row|<cell|FORPLAN>|<cell|5.874e-07>|<cell|2.084e-07>|<cell|4.979e-06>|<cell|SEBA>|<cell|2.919e-11>|<cell|5.729e-11>|<cell|1.448e-10>>|<row|<cell|GFRD-PNC>|<cell|4.558e-05>|<cell|1.052e-05>|<cell|4.363e-05>|<cell|SHARE1B>|<cell|3.367e-07>|<cell|1.339e-06>|<cell|3.578e-06>>|<row|<cell|GROW7>|<cell|1.276e-04>|<cell|4.906e-06>|<cell|1.024e-04>|<cell|SHARE2B>|<cell|2.142e-04>|<cell|2.014e-05>|<cell|6.146e-05>>|<row|<cell|ISRAEL>|<cell|1.422e-06>|<cell|1.336e-06>|<cell|1.404e-05>|<cell|STAIR>|<cell|5.549e-04>|<cell|8.566e-06>|<cell|2.861e-05>>|<row|<cell|STANDATA>|<cell|5.645e-08>|<cell|2.735e-07>|<cell|5.130e-06>|<cell|STANDGUB>|<cell|2.934e-08>|<cell|1.467e-07>|<cell|2.753e-06>>|<row|<cell|STOCFOR1>|<cell|6.633e-09>|<cell|9.701e-09>|<cell|4.811e-08>|<cell|VTP-BASE>|<cell|1.349e-10>|<cell|5.098e-11>|<cell|2.342e-10>>>>>>|Solving
  NETLIB LPs in 1000 iterations>

  <\bibliography|bib|tm-plain|ref>
    <\bib-list|4>
      <bibitem*|1><label|bib-ruiz2001scaling>Daniel Ruiz. <newblock>A scaling
      algorithm to equilibrate both rows and columns norms in matrices.
      <newblock><localize|Technical Report>, CM-P00040415, 2001.<newblock>

      <bibitem*|2><label|bib-ye2011interior>Yinyu Ye.
      <newblock><with|font-shape|italic|Interior point algorithms: theory and
      analysis>. <newblock>John Wiley & Sons, 2011.<newblock>

      <bibitem*|3><label|bib-ye2015first>Yinyu Ye. <newblock>On a first-order
      potential reduction algorithm for linear programming.
      <newblock>2015.<newblock>

      <bibitem*|4><label|bib-zhang2022drsom>Chuwen Zhang, Dongdong Ge, Bo
      Jiang<localize|, and >Yinyu Ye. <newblock>Drsom: a dimension reduced
      second-order method and preliminary analyses.
      <newblock><with|font-shape|italic|ArXiv preprint arXiv:2208.00208>,
      2022.<newblock>
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|3.2|4>>
    <associate|auto-11|<tuple|4|4>>
    <associate|auto-12|<tuple|1|5>>
    <associate|auto-13|<tuple|1|8>>
    <associate|auto-14|<tuple|1|8>>
    <associate|auto-15|<tuple|1|8>>
    <associate|auto-16|<tuple|1|?>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|1.2|2>>
    <associate|auto-4|<tuple|2|3>>
    <associate|auto-5|<tuple|2.1|3>>
    <associate|auto-6|<tuple|2.2|3>>
    <associate|auto-7|<tuple|2.3|3>>
    <associate|auto-8|<tuple|3|3>>
    <associate|auto-9|<tuple|3.1|3>>
    <associate|bib-ruiz2001scaling|<tuple|1|8>>
    <associate|bib-ye2011interior|<tuple|2|8>>
    <associate|bib-ye2015first|<tuple|3|8>>
    <associate|bib-zhang2022drsom|<tuple|4|8>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      ye2011interior

      ye2015first

      zhang2022drsom

      ruiz2001scaling
    </associate>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Solving NETLIB LPs
      in <with|color|<quote|red>|<with|font-series|<quote|bold>|1000>>
      iterations>|<pageref|auto-14>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>First-order Potential
      Reduction Method for LP <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>Dimension-reduced Potential
      Reduction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Accelerating
      the Dimension-reduced Potential Reduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Scaling
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>Line-search
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|2.3<space|2spc>Escaping the Local Optimum
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Algorithm
      Design> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Acceleration by negative
      curvature <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Direct computation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|1tab>|3.3<space|2spc>Scaled Hessian
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Algorithm
      Design> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Numerical
      Experiments> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>