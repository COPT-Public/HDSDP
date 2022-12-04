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

  <doc-data|<doc-title|First-order Potential Reduction
  Method>|<doc-author|<author-data|<author-name|Gwz>|<\author-affiliation>
    \;
  </author-affiliation>|<\author-affiliation>
    <date|>
  </author-affiliation>>>>

  <section|Dimension-reduced Method for Potential Reduction>

  In this section, we discuss the application of dimension-reduced method to
  potential reduction. For brevity, we for now only consider the primal
  potential reduction and focus on the simplex-constrained QP.

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<x>>>|<cell|<frac|1|2><around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>|<cell|\<backassign\>f<around*|(|<x>|)>>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsup|\<top\>><x>=1>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>>|<cell|>>>>
  </eqnarray*>

  and we adopt the potential function

  <\equation*>
    \<varphi\><around*|(|<x>|)>\<assign\>\<rho\> log
    <around*|(|f<around*|(|<x>|)>|)>-<big|sum><rsub|i=1><rsup|n>log
    x<rsub|i>,
  </equation*>

  whose gradient is given by

  <\equation*>
    <n>\<varphi\><around*|(|<x>|)>=<frac|\<rho\><n>f<around*|(|<x>|)>|f<around*|(|<x>|)>>-<X><rsup|-1><e>.
  </equation*>

  At each iteration, we evaluate the gradient
  <math|<n>\<varphi\><around*|(|<x><rsup|k>|)>>, let
  <math|\<Delta\>\<assign\><x><rsup|k+1>-<x><rsup|k>> and solve following
  subproblem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|\<Delta\>>>|<cell|<around*|\<langle\>|<n>\<varphi\><around*|(|<x><rsup|k>|)>,\<Delta\>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsup|\<top\>>\<Delta\>=0>|<cell|>>|<row|<cell|>|<cell|<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1>\<Delta\>|\<\|\|\>>\<leq\>\<beta\>.>|<cell|>>>>
  </eqnarray*>

  Starting from the basic potential reduction, we extend it by incorporating
  momentum term for faster convergence.

  <subsection|Two directions>

  In this section, we consider two direction extension of the potential
  reduction framework. In a word, by keeping track of one recent history
  iterate, we update

  <\eqnarray*>
    <tformat|<table|<row|<cell|<d><rsup|k>>|<cell|\<leftarrow\>>|<cell|\<alpha\><rsup|g><P><rsub|\<Delta\>><around*|[|\<nabla\>\<varphi\><around*|(|<x><rsup|k>|)>|]>+\<alpha\><rsup|m><around*|(|<x><rsup|k>-<x><rsup|k-1>|)>>>|<row|<cell|<x><rsup|k>>|<cell|\<leftarrow\>>|<cell|<x><rsup|k>+<d><rsup|k>>>>>
  </eqnarray*>

  where <math|<P><rsub|\<Delta\>><around*|[|\<cdummy\>|]>> is the orthogonal
  projection onto <math|<e><rsup|\<top\>><x>=0>. Note that we compute
  <math|\<alpha\><rsup|g>,\<alpha\><rsup|d>> through the following model

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<d>,\<alpha\><rsup|g>,\<alpha\><rsup|m>>>|<cell|<frac|1|2><d><rsup|\<top\>><H><d>+<h><rsup|\<top\>><d>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<around*|\<\|\|\>|<X><rsup|-1><d>|\<\|\|\>>\<leq\>\<Delta\>>|<cell|>>|<row|<cell|>|<cell|<d>=\<alpha\><rsup|g><g><rsup|k>+\<alpha\><rsup|m><m><rsup|k>>|<cell|>>>>
  </eqnarray*>

  where <math|<g><rsup|k>\<assign\><P><rsub|\<Delta\>><around*|[|\<nabla\>\<varphi\><around*|(|<x><rsup|k>|)>|]>>,
  <math|<m><rsup|k>\<assign\><x><rsup|k>-<x><rsup|k-1>>. Alternatively, we
  define <math|<G>\<assign\><matrix|<tformat|<table|<row|<cell|\|>|<cell|\|>>|<row|<cell|<g><rsup|k>>|<cell|<m><rsup|k>>>|<row|<cell|\|>|<cell|\|>>>>>,\<alpha\>=<matrix|<tformat|<table|<row|<cell|\<alpha\><rsup|g>>>|<row|<cell|\<alpha\><rsup|m>>>>>>>
  and <math|<d>=<G>\<alpha\>>, giving

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|\<alpha\>>>|<cell|<frac|1|2>\<alpha\><rsup|\<top\>><G><rsup|\<top\>><H><G><rsup|\<top\>>\<alpha\>+<h><rsup|\<top\>><G>\<alpha\>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<around*|\<\|\|\>|<X><rsup|-1><G>\<alpha\>|\<\|\|\>>\<leq\>\<Delta\>,>|<cell|>>>>
  </eqnarray*>

  or

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|\<alpha\>>>|<cell|<frac|1|2>\<alpha\><rsup|\<top\>><wide|<H>|~>\<alpha\>+<wide|<h>|~>\<alpha\>>|<cell|\<backassign\>m<around*|(|\<alpha\>|)>>>|<row|<cell|<text|subject
    to>>|<cell|<around*|\<\|\|\>|<M>\<alpha\>|\<\|\|\>>\<leq\>\<Delta\>>|<cell|>>>>
  </eqnarray*>

  for

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|<H>|~>>|<cell|\<assign\>>|<cell|<matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<g><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<langle\>|<g><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>|<row|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>>>>>>|<row|<cell|<wide|<h>|~>>|<cell|\<assign\>>|<cell|<matrix|<tformat|<table|<row|<cell|<around*|\<\|\|\>|<g><rsup|k>|\<\|\|\>><rsup|2>>>|<row|<cell|<around*|\<langle\>|<g><rsup|k>,<m><rsup|k>|\<rangle\>>>>>>>>>|<row|<cell|<M>>|<cell|\<assign\>>|<cell|<matrix|<tformat|<table|<row|<cell|<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><g><rsup|k>|\<\|\|\>><rsup|2>>|<cell|<around*|\<langle\>|<g><rsup|k>,<around*|(|<X><rsup|k>|)><rsup|-2><m><rsup|k>|\<rangle\>>>>|<row|<cell|<around*|\<langle\>|<m><rsup|k>,<around*|(|<X><rsup|k>|)><rsup|-2><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><m><rsup|k>|\<\|\|\>><rsup|2>>>>>>.>>>>
  </eqnarray*>

  \;

  Note that <math|\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)>=-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>>|f<around*|(|<x><rsup|k>|)><rsup|2>>+\<rho\><frac|<A><rsup|\<top\>><A>|f<around*|(|<x><rsup|k>|)>>+<around*|(|<X><rsup|k>|)><rsup|-2>>
  and we evaluate the above relations via

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|\<langle\>|<a>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><a>|\<rangle\>>>|<cell|=>|<cell|<around*|\<langle\>|<a>,-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><a>|f<around*|(|<x><rsup|k>|)><rsup|2>>|\<rangle\>>+<frac|<around*|\<\|\|\>|<A><a>|\<\|\|\>><rsup|2>|f<around*|(|<x><rsup|k>|)>>+<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><a>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|-\<rho\><around*|(|<frac|<n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><a>|f<around*|(|<x><rsup|k>|)>>|)><rsup|2>+<frac|<around*|\<\|\|\>|<A><a>|\<\|\|\>><rsup|2>|f<around*|(|<x><rsup|k>|)>>+<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><a>|\<\|\|\>><rsup|2>>>|<row|<cell|<around*|\<langle\>|<a>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><b>|\<rangle\>>>|<cell|=>|<cell|<around*|\<langle\>|<a>,-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><b>|f<around*|(|<x><rsup|k>|)><rsup|2>>|\<rangle\>>+<frac|<around*|\<langle\>|<A><a>,<A><b>|\<rangle\>>|f<around*|(|<x><rsup|k>|)>>+<around*|\<langle\>|<a>,<around*|(|<X><rsup|k>|)><rsup|-2><b>|\<rangle\>>>>|<row|<cell|>|<cell|=>|<cell|-\<rho\><around*|(|<frac|<n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><a>|f<around*|(|<x><rsup|k>|)>>|)><around*|(|<frac|<n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><b>|f<around*|(|<x><rsup|k>|)>>|)>+<frac|<around*|\<langle\>|<A><a>,<A><b>|\<rangle\>>|f<around*|(|<x><rsup|k>|)>>+<around*|\<langle\>|<a>,<around*|(|<X><rsup|k>|)><rsup|-2><b>|\<rangle\>>.>>>>
  </eqnarray*>

  To ensure feasibility, we always choose <math|\<Delta\>\<leq\>1> and adjust
  it based on the trust-region rule.

  <section|Potential Reduction for LP>

  In this section, we discuss the potential reduction method on LP HSD model.

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<x>\<in\>\<bbb-R\><rsup|n>>>|<cell|<c><rsup|\<top\>><x>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<A><x>=<b>>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|max<rsub|<y>\<in\>\<bbb-R\><rsup|m>>>|<cell|<b><rsup|\<top\>><y>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<A><rsup|\<top\>><y>+<s>=<c>>|<cell|>>|<row|<cell|>|<cell|<s>\<geq\><0>>|<cell|>>>>
  </eqnarray*>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|<A><x>-<b>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|-<A><rsup|\<top\>><y>-<s>+<c>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|<b><rsup|\<top\>><y>-<c><rsup|\<top\>><x>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|<e><rsub|n><rsup|\<top\>><x>+<e><rsub|n><rsup|\<top\>><s>+\<kappa\>+\<tau\>>|<cell|=>|<cell|1>>>>
  </eqnarray*>

  <subsection|Potential Reduction for HSD>

  In this section we consider the original HSD formulation

  <\eqnarray*>
    <tformat|<table|<row|<cell|<A><x>-<b>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|-<A><rsup|\<top\>><y>-<s>+<c>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|<b><rsup|\<top\>><y>-<c><rsup|\<top\>><x>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|<e><rsub|n><rsup|\<top\>><x>+<e><rsub|n><rsup|\<top\>><s>+\<kappa\>+\<tau\>>|<cell|=>|<cell|1>>>>
  </eqnarray*>

  and we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|<D><matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|<A>>|<cell|<0><rsub|m\<times\>n>>|<cell|<0><rsub|m\<times\>1>>|<cell|-<b>>>|<row|<cell|-<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|<c>>>|<row|<cell|<b><rsup|\<top\>>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>|<cell|0>>>>><math-bf|E><matrix|<tformat|<table|<row|<cell|<y>>>|<row|<cell|<x>>>|<row|<cell|<s>>>|<row|<cell|\<kappa\>>>|<row|<cell|\<tau\>>>>>>>|<cell|=>|<cell|<0>>>|<row|<cell|<e><rsub|n><rsup|\<top\>><x>+<e><rsub|n><rsup|\<top\>><s>+\<kappa\>+\<tau\>>|<cell|=>|<cell|1.>>>>
  </eqnarray*>

  In this method, the dual variable <math|<y>> is free and needs special
  treatment. First we consider the potential function\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|f<around*|(|<x>,<y>,<s>,\<kappa\>,\<tau\>|)>>|<cell|\<assign\>>|<cell|<frac|1|2><around*|\<\|\|\>|<wide|<A>|~><u>|\<\|\|\>><rsup|2>>>|<row|<cell|\<varphi\><around*|(|<x>,<y>,<s>,\<kappa\>,\<tau\>|)>>|<cell|\<assign\>>|<cell|\<rho\>
    log <around*|(|f<around*|(|<u>|)>|)>-B<around*|(|<x>|)>-B<around*|(|<s>|)>-log
    \<kappa\>-log \<tau\>>>>>
  </eqnarray*>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<nabla\>f<around*|(|<x>,<y>,<s>,\<kappa\>,\<tau\>|)>>>|<row|<cell|>|<cell|=>|<cell|<wide|<A>|~><rsup|\<top\>><wide|<A>|~><u>>>|<row|<cell|>|<cell|=>|<cell|<matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|<A>>|<cell|<0><rsub|m\<times\>n>>|<cell|<0><rsub|m\<times\>1>>|<cell|-<b>>>|<row|<cell|-<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|<c>>>|<row|<cell|<b><rsup|\<top\>>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>|<cell|0>>>>><rsup|\<top\>><matrix|<tformat|<table|<row|<cell|<A><x>-<b>\<tau\>>|<cell|\<backassign\><r><rsub|1>>>|<row|<cell|-<A><rsup|\<top\>><y>-<s>+<c>\<tau\>>|<cell|\<backassign\><r><rsub|2>>>|<row|<cell|<b><rsup|\<top\>><y>-<c><rsup|\<top\>><x>-\<kappa\>>|<cell|\<backassign\>r<rsub|3>>>>>>>>|<row|<cell|>|<cell|=>|<cell|<matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|-<A>>|<cell|<b>>>|<row|<cell|<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<c>>>|<row|<cell|<0><rsub|n\<times\>m>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>>|<row|<cell|<0><rsub|1\<times\>m>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>>|<row|<cell|-<b><rsup|\<top\>>>|<cell|<c><rsup|\<top\>>>|<cell|0>>>>><matrix|<tformat|<table|<row|<cell|<r><rsub|1>>>|<row|<cell|<r><rsub|2>>>|<row|<cell|r<rsub|3>>>>>>=<matrix|<tformat|<table|<row|<cell|-<A><r><rsub|2>+<b>r<rsub|3>>>|<row|<cell|<A><rsup|\<top\>><r><rsub|1>-<c>r<rsub|3>>>|<row|<cell|-<r><rsub|2>>>|<row|<cell|-r<rsub|3>>>|<row|<cell|-<b><rsup|\<top\>><r><rsub|1>+<c><rsup|\<top\>><r><rsub|2>>>>>>.>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\>\<varphi\><around*|(|<u>|)>>|<cell|=>|<cell|<frac|\<rho\><n>f<around*|(|<u>|)>|f<around*|(|<u>|)>>-<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>>|<row|<cell|<X><rsup|-1><e>>>|<row|<cell|<bs><rsup|-1><e>>>|<row|<cell|\<kappa\><rsup|-1>>>|<row|<cell|\<tau\><rsup|-1>>>>>>>>|<row|<cell|\<nabla\><rsup|2><rsub|<u>,<u>>\<varphi\><around*|(|<u>|)>>|<cell|=>|<cell|-<frac|\<rho\><n>f<around*|(|<u>|)><n>f<around*|(|<u>|)><rsup|\<top\>>|f<around*|(|<u>|)><rsup|2>>+\<rho\><frac|<wide|<A>|~><rsup|\<top\>><wide|<A>|~>|f<around*|(|<u>|)>>+diag<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>>|<row|<cell|<X><rsup|-2><e>>>|<row|<cell|<bs><rsup|-2><e>>>|<row|<cell|\<kappa\><rsup|-2>>>|<row|<cell|\<tau\><rsup|-2>>>>>>.>>>>
  </eqnarray*>

  <subsection|Acceleration by negative curvature>

  In this section, we discuss how to find the negative curvature of the
  Hessian to help accelerate algorithm convergence. More specifically, we
  consider the following problem

  <\equation*>
    \<lambda\><rsub|min><around*|{|\<nabla\><rsup|2>\<varphi\><around*|(|<u>=<around*|(|<x>,<y>|)>|)>=<frac|2\<rho\><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><A><rsup|\<top\>><A><u><u><rsup|\<top\>><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4>>+<small|<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X><rsup|-2>>>>>>>|}>.
  </equation*>

  And we wish to solve the eigen-problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<around*|\<\|\|\>|<v>|\<\|\|\>>=1>>|<cell|<v><rsup|\<top\>><around*|{|<frac|2\<rho\><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><A><rsup|\<top\>><A><u><u><rsup|\<top\>><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4>>+<small|<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X><rsup|-2>>>>>>>|}><v>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsup|\<top\>><v><rsub|<x>>=0.>|<cell|>>>>
  </eqnarray*>

  In general there are two ways to compute a valid direction. The first
  method approaches the problem directly and uses Lanczos iteration to find
  the negative eigen-value of <math|<n><rsup|2>\<varphi\>>. As for the second
  approach, we apply the scaling matrix <math|<bs>\<assign\><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>>>
  and solve

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<around*|\<\|\|\>|<bs><v>|\<\|\|\>>=1>>|<cell|<v><rsup|\<top\>><small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>>><around*|{|<frac|2\<rho\><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><A><rsup|\<top\>><A><u><u><rsup|\<top\>><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4>>+<small|<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X><rsup|-2>>>>>>>|}><small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>>><v>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<x><rsup|\<top\>><v><rsub|<x>>=0.>|<cell|>>>>
  </eqnarray*>

  To improve the conditioning of the Hessian, we replace
  <math|<around*|\<\|\|\>|<bs><v>|\<\|\|\>>=1> by
  <math|<around*|\<\|\|\>|<v>|\<\|\|\>>=1> and arrive at

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<around*|\<\|\|\>|<v>|\<\|\|\>>=1>>|<cell|<v><rsup|\<top\>><small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>>><around*|{|<frac|2\<rho\><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><A><rsup|\<top\>><A><u><u><rsup|\<top\>><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4>>+<small|<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X><rsup|-2>>>>>>>|}><small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>>><v>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<x><rsup|\<top\>><v><rsub|<x>>=0.>|<cell|>>>>
  </eqnarray*>

  Another trick we apply is to ignore the variables which are predicted to be
  nonbasic in the optimal solution so that the Hessian computation can be
  greatly simplified.

  <subsection|Direct computation>

  When evaluating the Hessian, it is possible that the matrix is
  ill-conditioned. Hence we need to consider the following relation

  <\equation*>
    <small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<e><rsub|n><e><rsub|n><rsup|\<top\>>/n>>>>>><around*|[|<frac|\<rho\><A><rsup|\<top\>><A>|f>-<frac|\<rho\><A><rsup|\<top\>><A><u><u><rsup|\<top\>><A><rsup|\<top\>><A>|f<rsup|2>>+<small|<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X><rsup|-2>>>>>>>|]><small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<e><rsub|n><e><rsub|n><rsup|\<top\>>/n>>>>>>
  </equation*>

  <math|<around*|(|<I><rsub|n>-<e><rsub|n><e><rsub|n><rsup|\<top\>>/n|)><X><rsup|-2><around*|(|<I><rsub|n>-<e><rsub|n><e><rsub|n><rsup|\<top\>>/n|)>=<X><rsup|-2><v>-<frac|<X><rsup|-1><e><rsub|n><rsup|\<top\>>|n><v>-<frac|<e><rsub|n><e><rsub|n><rsup|\<top\>>|n><X><rsup|-1><v>+<frac|<e><rsub|n><e><rsub|n><rsup|\<top\>>|n<rsup|2>><e><rsub|n><rsup|\<top\>><X><rsup|-2><e><rsub|n>>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<v>>|<cell|\<leftarrow\>>|<cell|<matrix|<tformat|<table|<row|<cell|<v><rsub|<y>>>>|<row|<cell|<v><rsub|<x>>-<around*|(|<e><rsub|n><rsup|\<top\>><v>|)><v>/n>>>>>>>|<row|<cell|<u><rsub|1>>|<cell|\<leftarrow\>>|<cell|<matrix|<tformat|<table|<row|<cell|<0>>>|<row|<cell|<v><rsub|<x>>-<e><rsub|n><e><rsub|n><rsup|\<top\>>/n>>>>>>>|<row|<cell|<u><rsub|2>>|<cell|\<leftarrow\>>|<cell|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<e><rsub|n><e><rsub|n><rsup|\<top\>>/n>>>>><A><rsup|\<top\>><A><v>>>|<row|<cell|<u><rsub|3>>|<cell|\<leftarrow\>>|<cell|<around*|(|<g><rsup|\<top\>><v>|)><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<e><rsub|n><e><rsub|n><rsup|\<top\>>/n>>>>><g>>>|<row|<cell|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4><M><v>>|<cell|\<leftarrow\>>|<cell|f<rsup|2><u><rsub|1>+2\<rho\>f<u><rsub|2>-4\<rho\><u><rsub|3>>>>>
  </eqnarray*>

  <subsection|Scaled Hessian>

  <\equation*>
    <M>\<assign\><small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><x><rsup|\<top\>>/<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>>>>>><around*|[|<frac|2\<rho\><bs><A><rsup|\<top\>><A><bs>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><bs><A><rsup|\<top\>><A><u><u><rsup|\<top\>><A><rsup|\<top\>><A><bs>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4>>+<small|<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>>>>>>>|]><small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><x><rsup|\<top\>>/<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>>>>>>
  </equation*>

  In the computation of scaled Hessian, we implement matrix-vector product
  <math|<M><v>> as follows

  <\eqnarray*>
    <tformat|<table|<row|<cell|<x><rprime|'>>|<cell|\<leftarrow\>>|<cell|<frac|<x>|<around*|\<\|\|\>|<x>|\<\|\|\>>>>>|<row|<cell|<v>>|<cell|\<leftarrow\>>|<cell|<matrix|<tformat|<table|<row|<cell|<v><rsub|<y>>>>|<row|<cell|<v><rsub|<x>>-<around*|(|<x><rprime|'><rsup|\<top\>><v><rsub|<x>>|)><x><rprime|'>>>>>>>>|<row|<cell|<u><rsub|1>>|<cell|\<leftarrow\>>|<cell|<matrix|<tformat|<table|<row|<cell|<0>>>|<row|<cell|<v><rsub|<x>>-<around*|(|<x><rprime|'><rsup|\<top\>><v><rsub|<x>>|)><x><rprime|'>>>>>>>>|<row|<cell|<u><rsub|2>>|<cell|\<leftarrow\>>|<cell|<small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><rprime|'><x><rprime|'><rsup|\<top\>>>>>>>><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><A><rsup|\<top\>><A><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><v>>>|<row|<cell|<u><rsub|3>>|<cell|\<leftarrow\>>|<cell|<g><rsup|\<top\>><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><v><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><rprime|'><x><rprime|'><rsup|\<top\>>>>>>><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><g>>>|<row|<cell|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4><M><v>>|<cell|\<leftarrow\>>|<cell|f<rsup|2><u><rsub|1>+\<rho\>f<u><rsub|2>-\<rho\><u><rsub|3>>>>>
  </eqnarray*>

  <section|Numerical Experiments>

  <big-table|<small|<block|<tformat|<cwith|2|13|1|4|cell-valign|c>|<cwith|1|-1|1|-1|cell-halign|c>|<cwith|1|13|5|8|cell-halign|c>|<cwith|2|7|6|8|cell-valign|c>|<cwith|2|7|5|5|cell-valign|c>|<cwith|14|18|1|4|cell-halign|c>|<cwith|14|18|1|4|cell-halign|c>|<cwith|14|18|2|4|cell-valign|c>|<cwith|14|18|1|1|cell-valign|c>|<cwith|2|3|5|8|cell-halign|c>|<cwith|2|3|5|8|cell-halign|c>|<cwith|2|2|6|8|cell-valign|c>|<cwith|2|2|5|5|cell-valign|c>|<cwith|4|13|5|5|cell-halign|c>|<cwith|4|8|5|5|cell-halign|c>|<cwith|4|13|5|5|cell-valign|c>|<cwith|4|13|6|8|cell-valign|c>|<cwith|14|21|5|5|cell-rborder|0ln>|<cwith|14|21|5|5|cell-valign|c>|<cwith|14|18|6|8|cell-valign|c>|<cwith|19|20|1|4|cell-valign|c>|<cwith|19|19|5|8|cell-halign|c>|<cwith|19|19|5|8|cell-valign|c>|<cwith|20|20|1|4|cell-valign|c>|<cwith|20|20|5|8|cell-valign|c>|<cwith|1|-1|1|-1|cell-tborder|1ln>|<cwith|1|-1|1|-1|cell-bborder|1ln>|<cwith|1|-1|1|-1|cell-lborder|0ln>|<cwith|1|-1|1|-1|cell-rborder|0ln>|<cwith|1|-1|5|5|cell-tborder|1ln>|<cwith|1|-1|5|5|cell-bborder|1ln>|<cwith|1|-1|5|5|cell-lborder|1ln>|<cwith|1|-1|5|5|cell-rborder|1ln>|<cwith|1|-1|4|4|cell-rborder|1ln>|<cwith|1|-1|6|6|cell-lborder|1ln>|<table|<row|<cell|Problem>|<cell|PInfeas>|<cell|DInfeas.>|<cell|Compl.>|<cell|Problem>|<cell|PInfeas>|<cell|DInfeas.>|<cell|Compl.>>|<row|<cell|DLITTLE>|<cell|1.347e-10>|<cell|2.308e-10>|<cell|2.960e-09>|<cell|KB2>|<cell|5.455e-11>|<cell|6.417e-10>|<cell|7.562e-11>>|<row|<cell|AFIRO>|<cell|7.641e-11>|<cell|7.375e-11>|<cell|3.130e-10>|<cell|LOTFI>|<cell|2.164e-09>|<cell|4.155e-09>|<cell|8.663e-08>>|<row|<cell|AGG2>|<cell|3.374e-08>|<cell|4.859e-08>|<cell|6.286e-07>|<cell|MODSZK1>|<cell|1.527e-06>|<cell|5.415e-05>|<cell|2.597e-04>>|<row|<cell|AGG3>|<cell|2.248e-05>|<cell|1.151e-06>|<cell|1.518e-05>|<cell|RECIPELP>|<cell|5.868e-08>|<cell|6.300e-08>|<cell|1.285e-07>>|<row|<cell|BANDM>|<cell|2.444e-09>|<cell|4.886e-09>|<cell|3.769e-08>|<cell|SC105>|<cell|7.315e-11>|<cell|5.970e-11>|<cell|2.435e-10>>|<row|<cell|BEACONFD>|<cell|5.765e-12>|<cell|9.853e-12>|<cell|1.022e-10>|<cell|SC205>|<cell|6.392e-11>|<cell|5.710e-11>|<cell|2.650e-10>>|<row|<cell|BLEND>|<cell|2.018e-10>|<cell|3.729e-10>|<cell|1.179e-09>|<cell|SC50A>|<cell|1.078e-05>|<cell|6.098e-06>|<cell|4.279e-05>>|<row|<cell|BOEING2>|<cell|1.144e-07>|<cell|1.110e-08>|<cell|2.307e-07>|<cell|SC50B>|<cell|4.647e-11>|<cell|3.269e-11>|<cell|1.747e-10>>|<row|<cell|BORE3D>|<cell|2.389e-08>|<cell|5.013e-08>|<cell|1.165e-07>|<cell|SCAGR25>|<cell|1.048e-07>|<cell|5.298e-08>|<cell|1.289e-06>>|<row|<cell|BRANDY>|<cell|2.702e-05>|<cell|7.818e-06>|<cell|1.849e-05>|<cell|SCAGR7>|<cell|1.087e-07>|<cell|1.173e-08>|<cell|2.601e-07>>|<row|<cell|CAPRI>|<cell|7.575e-05>|<cell|4.488e-05>|<cell|4.880e-05>|<cell|SCFXM1>|<cell|4.323e-06>|<cell|5.244e-06>|<cell|8.681e-06>>|<row|<cell|E226>|<cell|2.656e-06>|<cell|4.742e-06>|<cell|2.512e-05>|<cell|SCORPION>|<cell|1.674e-09>|<cell|1.892e-09>|<cell|1.737e-08>>|<row|<cell|FINNIS>|<cell|8.577e-07>|<cell|8.367e-07>|<cell|1.001e-05>|<cell|SCTAP1>|<cell|5.567e-07>|<cell|8.430e-07>|<cell|5.081e-06>>|<row|<cell|FORPLAN>|<cell|5.874e-07>|<cell|2.084e-07>|<cell|4.979e-06>|<cell|SEBA>|<cell|2.919e-11>|<cell|5.729e-11>|<cell|1.448e-10>>|<row|<cell|GFRD-PNC>|<cell|4.558e-05>|<cell|1.052e-05>|<cell|4.363e-05>|<cell|SHARE1B>|<cell|3.367e-07>|<cell|1.339e-06>|<cell|3.578e-06>>|<row|<cell|GROW7>|<cell|1.276e-04>|<cell|4.906e-06>|<cell|1.024e-04>|<cell|SHARE2B>|<cell|2.142e-04>|<cell|2.014e-05>|<cell|6.146e-05>>|<row|<cell|ISRAEL>|<cell|1.422e-06>|<cell|1.336e-06>|<cell|1.404e-05>|<cell|STAIR>|<cell|5.549e-04>|<cell|8.566e-06>|<cell|2.861e-05>>|<row|<cell|STANDATA>|<cell|5.645e-08>|<cell|2.735e-07>|<cell|5.130e-06>|<cell|STANDGUB>|<cell|2.934e-08>|<cell|1.467e-07>|<cell|2.753e-06>>|<row|<cell|STOCFOR1>|<cell|6.633e-09>|<cell|9.701e-09>|<cell|4.811e-08>|<cell|VTP-BASE>|<cell|1.349e-10>|<cell|5.098e-11>|<cell|2.342e-10>>>>>>|Solving
  NETLIB LPs in <with|color|red|<with|font-series|bold|1000>> iterations>

  \;

  <section|Algorithm Design>

  In this section, we discuss the design of the potential-reduction based
  solver. To allow further extension, we design the solver to solve general
  problem

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

    Used for problem solving.
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

  \;

  <\eqnarray*>
    <tformat|<table|<row|<cell|<A><x>-<b>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|-<A><rsup|\<top\>><y>-<s>+<c>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|<b><rsup|\<top\>><y>-<c><rsup|\<top\>><x>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|<e><rsub|n><rsup|\<top\>><x>+<e><rsub|n><rsup|\<top\>><s>+\<kappa\>+\<tau\>>|<cell|=>|<cell|2n+2>>>>
  </eqnarray*>

  and we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|<matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|<A>>|<cell|<0><rsub|m\<times\>n>>|<cell|<0><rsub|m\<times\>1>>|<cell|-<b>>>|<row|<cell|-<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|<c>>>|<row|<cell|<b><rsup|\<top\>>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>|<cell|0>>>>><matrix|<tformat|<table|<row|<cell|<y>>>|<row|<cell|<x>>>|<row|<cell|<s>>>|<row|<cell|\<kappa\>>>|<row|<cell|\<tau\>>>>>>>|<cell|=>|<cell|<0>>>|<row|<cell|<e><rsub|n><rsup|\<top\>><x>+<e><rsub|n><rsup|\<top\>><s>+\<kappa\>+\<tau\>>|<cell|=>|<cell|1.>>>>
  </eqnarray*>

  There are some basic operations to implement

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

  <strong|Hessian (no actual setup)>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|<A>|~><rsup|\<top\>><wide|<A>|~>>|<cell|=>|<cell|<matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|-<A>>|<cell|<b>>>|<row|<cell|<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<c>>>|<row|<cell|<0><rsub|n\<times\>m>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>>|<row|<cell|<0><rsub|1\<times\>m>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>>|<row|<cell|-<b><rsup|\<top\>>>|<cell|<c><rsup|\<top\>>>|<cell|0>>>>><matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|<A>>|<cell|<0><rsub|m\<times\>n>>|<cell|<0><rsub|m\<times\>1>>|<cell|-<b>>>|<row|<cell|-<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|<c>>>|<row|<cell|<b><rsup|\<top\>>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>|<cell|0>>>>>>>|<row|<cell|>|<cell|=>|<cell|<matrix|<tformat|<table|<row|<cell|<A><A><rsup|\<top\>>+<b><b><rsup|\<top\>>>|<cell|-<b><c><rsup|\<top\>>>|<cell|<A>>|<cell|-<b>>|<cell|-<A><c>>>|<row|<cell|-<c><b><rsup|\<top\>>>|<cell|<A><rsup|\<top\>><A>>|<cell|<0><rsub|n\<times\>n>>|<cell|<c>>|<cell|-<A><rsup|\<top\>><b>>>|<row|<cell|<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|-<c>>>|<row|<cell|-<b><rsup|\<top\>>>|<cell|<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|<0><rsub|1\<times\>1>>|<cell|<0><rsub|1\<times\>1>>>|<row|<cell|-<c><rsup|\<top\>><A><rsup|\<top\>>>|<cell|-<b><rsup|\<top\>><A>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>1>>|<cell|<around*|\<\|\|\>|<b>|\<\|\|\>><rsup|2>+<around*|\<\|\|\>|<c>|\<\|\|\>><rsup|2>>>>>>>>>>
  </eqnarray*>

  <\equation*>
    <n><rsup|2>\<varphi\>=-<frac|\<rho\><n>f<n>f<rsup|\<top\>>|f<rsup|2>>+<frac|\<rho\><wide|<A>|~><rsup|\<top\>><wide|<A>|~>|f>+<small|<matrix|<tformat|<table|<row|<cell|<X><rsup|-2>>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|<0><rsub|m>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|<bs><rsup|-2>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|\<kappa\><rsup|-2>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|\<tau\><rsup|-2>>>>>>>.
  </equation*>

  <strong|Hessian-vector (with projection)>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<u>>|<cell|=>|<cell|<x>-<frac|<e><rsup|\<top\>><x>|n>\<cdummy\><e>>>|<row|<cell|<n><rsup|2>\<varphi\><u>>|<cell|=>|<cell|-<frac|\<rho\><around*|(|<n>f<rsup|\<top\>><u>|)>|f<rsup|2>><n>f+<frac|\<rho\>|f><wide|<A>|~><rsup|\<top\>><around*|(|<wide|<A>|~><u>|)>+<small|<matrix|<tformat|<table|<row|<cell|<X><rsup|-2>>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|<0><rsub|m\<times\>m>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|<bs><rsup|-2>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|\<kappa\><rsup|-2>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|\<tau\><rsup|-2>>>>>>><u>.>>>>
  </eqnarray*>

  <strong|Minimal eigenvalue>

  To evaluate the minimum eigen-value of <math|<P><rsub|\<Delta\>><n><rsup|2>\<varphi\><P><rsub|\<Delta\>>>
  and the corresponding eigen-vector\ 

  \;

  <\equation*>
    <X><n><rsup|2>\<varphi\><X>=-<frac|4\<rho\><X><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A><X>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|4>>+<frac|2\<rho\><X><A><rsup|\<top\>><A><X>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>+<I>
  </equation*>

  Note that\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<X><n><rsup|2>\<varphi\><around*|(|<u>|)><X>>>|<row|<cell|>|<cell|=>|<cell|-<frac|\<rho\><X><n>f<around*|(|<u>|)><n>f<around*|(|<u>|)><rsup|\<top\>><X>|f<around*|(|<u>|)><rsup|2>>+<frac|\<rho\><X><A><rsup|\<top\>><A><X>|f<around*|(|<u>|)>>+<D>.>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<v>>>|<cell|<around*|\<langle\>|<X><v>,<n><rsup|2>\<varphi\><around*|(|<u>|)><X><v>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsup|\<top\>><X><v>=0>|<cell|>>|<row|<cell|>|<cell|<around*|\<\|\|\>|<X><v>|\<\|\|\>>=1>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<v>>>|<cell|<around*|\<langle\>|<v>,<around*|(|<I>-<frac|<x><x><rsup|\<top\>>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>|)><around*|(|<X><n><rsup|2>\<varphi\><around*|(|<u>|)><X>|)><around*|(|<I>-<frac|<x><x><rsup|\<top\>>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>|)><v>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<around*|\<\|\|\>|<X><v>|\<\|\|\>>=1>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<around*|(|<I>-<frac|<x><x><rsup|\<top\>>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>|)><H><around*|(|<I>-<frac|<x><x><rsup|\<top\>>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>|)>>>|<row|<cell|>|<cell|=>|<cell|<H>-<frac|<x><x><rsup|\<top\>>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>><H>-<H><frac|<x><x><rsup|\<top\>>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>+<x><rsup|\<top\>><H><x><frac|<x><x><rsup|\<top\>>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|4>>>>>>
  </eqnarray*>

  <section|Analysis>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\><rsup|2>\<varphi\><around*|(|<x>|)>>|<cell|=>|<cell|-<frac|\<rho\><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>>|f<around*|(|<x>|)><rsup|2>>+\<rho\><frac|<A><rsup|\<top\>><A>|f<around*|(|<x>|)>>+<X><rsup|-2>>>|<row|<cell|<n>\<varphi\><around*|(|<x>|)>>|<cell|=>|<cell|<frac|\<rho\><n>f<around*|(|<x>|)>|f<around*|(|<x>|)>>-<X><rsup|-1><e>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|<A>>|<cell|<0><rsub|m\<times\>n>>|<cell|<0><rsub|m\<times\>1>>|<cell|-<b>>>|<row|<cell|-<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|<c>>>|<row|<cell|<b><rsup|\<top\>>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>|<cell|0>>>>><matrix|<tformat|<table|<row|<cell|<y>>>|<row|<cell|<x>>>|<row|<cell|<s>>>|<row|<cell|\<kappa\>>>|<row|<cell|\<tau\>>>>>>>|<cell|=>|<cell|<0>>>|<row|<cell|<e><rsub|n><rsup|\<top\>><x>+<e><rsub|n><rsup|\<top\>><s>+\<kappa\>+\<tau\>>|<cell|=>|<cell|1.>>>>
  </eqnarray*>

  <\equation*>
    \<nabla\><rsup|2>\<varphi\><around*|(|<x>|)>=-<frac|4\<rho\><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|4>>+<frac|2\<rho\><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>+<X><rsup|-2>
  </equation*>

  <\equation*>
    <X><n><rsup|2>\<varphi\><around*|(|<x>|)><X>=-<frac|4\<rho\><X><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A><X>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|4>>+<frac|2\<rho\><X><A><rsup|\<top\>><A><X>|<around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>+<I>
  </equation*>

  <new-page*>

  <section|General Potential Method>

  <subsection|Second-order Potential Reduction>

  In this section, we consider the second order potential reduction method,
  where we update the iterates by\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<d><rsup|k>>|<cell|=>|<cell|<below|arg
    min|<around*|\<\|\|\>|<d>|\<\|\|\>>\<leq\>\<beta\>,<x><rsup|\<top\>><d>=0>
    <around*|{|<around*|\<langle\>|<X><rsup|k><n>\<phi\><around*|(|<x><rsup|k>|)>,<d>|\<rangle\>>+<frac|1|2><around*|\<langle\>|<d>,<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k><d>|\<rangle\>>|}>>>|<row|<cell|<x><rsup|k+1>>|<cell|=>|<cell|<x><rsup|k>+<X><rsup|k><d><rsup|k>>>>>
  </eqnarray*>

  First, by the optimality condition of the trust region subproblem, we have,
  for some <math|\<lambda\><rsup|k>\<geq\>0,\<mu\><rsup|k>> that

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|(|<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k>+\<lambda\><rsup|k><I>|)><d><rsup|k>-\<mu\><x><rsup|k>>|<cell|=>|<cell|-<X><rsup|k><n>\<phi\><around*|(|<x><rsup|k>|)>>>|<row|<cell|\<mu\><rsup|k><around*|(|<around*|\<\|\|\>|<d><rsup|k>|\<\|\|\>>-\<beta\>|)>>|<cell|=>|<cell|0>>|<row|<cell|<x><rsup|\<top\>><d><rsup|k>>|<cell|=>|<cell|0>>|<row|<cell|<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k>+\<lambda\><rsup|k><I>>|<cell|\<succeq\><rsub|<x>>>|<cell|<0>>>>>
  </eqnarray*>

  Assume that <math|<around*|\<\|\|\>|<d><rsup|k>|\<\|\|\>>=\<beta\>> and
  define\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|<x>,\<mu\>|)>>|<cell|\<assign\>>|<cell|<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k><d><rsup|k>+<X><rsup|k><n>\<phi\><around*|(|<x><rsup|k>|)>-\<mu\><rsup|k><x><rsup|k>.>>>>
  </eqnarray*>

  Then it follows that

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<lambda\><rsup|k><d><rsup|k>>|<cell|=>|<cell|-p<around*|(|<x><rsup|k>,\<mu\><rsup|k>|)>>>>>
  </eqnarray*>

  and we successively deduce that

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<around*|\<langle\>|<X><rsup|k><n>\<phi\><around*|(|<x><rsup|k>|)>,<d><rsup|k>|\<rangle\>>+<frac|1|2><around*|\<langle\>|<d><rsup|k>,<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k><d><rsup|k>|\<rangle\>>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<langle\>|-\<lambda\><rsup|k><d><rsup|k>-<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k><d><rsup|k>+\<mu\><rsup|k><x><rsup|k>,<d><rsup|k>|\<rangle\>>+<frac|1|2><around*|\<langle\>|<d><rsup|k>,<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k><d><rsup|k>|\<rangle\>>>>|<row|<cell|>|<cell|=>|<cell|-\<lambda\><rsup|k><around*|\<\|\|\>|<d><rsup|k>|\<\|\|\>><rsup|2>-<frac|1|2><around*|\<langle\>|<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k><d><rsup|k>,<d><rsup|k>|\<rangle\>>.>>>>
  </eqnarray*>

  Since <math|<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k>\<succeq\><rsub|<x>>-\<lambda\><rsup|k><I>>,
  we have

  <\equation*>
    <around*|\<langle\>|<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k><d><rsup|k>,<d><rsup|k>|\<rangle\>>\<geq\>-<around*|\<\|\|\>|<d><rsup|k>|\<\|\|\>><rsup|2>
  </equation*>

  and that

  <\equation*>
    <around*|\<langle\>|<X><rsup|k><n>\<phi\><around*|(|<x><rsup|k>|)>,<d><rsup|k>|\<rangle\>>+<frac|1|2><around*|\<langle\>|<d><rsup|k>,<X><rsup|k>\<nabla\><rsup|2>\<phi\><around*|(|<x><rsup|k>|)><X><rsup|k><d><rsup|k>|\<rangle\>>\<leq\>-<frac|\<lambda\><rsup|k>|2><around*|\<\|\|\>|<d><rsup|k>|\<\|\|\>><rsup|2>=-<frac|\<lambda\><rsup|k>\<beta\><rsup|2>|2>.
  </equation*>

  Next we derive the reduction of the potential function. It follows
  naturally that

  <\equation*>
    <big|sum><rsub|i=1><rsup|n>log x<rsub|i>-<big|sum><rsub|i=1><rsup|n>log
    <around*|(|x<rsub|i>+x<rsub|i>d<rsub|i>|)>\<leq\>-<around*|\<langle\>|<e>,<d>|\<rangle\>>+<frac|\<beta\><rsup|2>|2<around*|(|1-\<beta\>|)>>.
  </equation*>

  \;

  Then we bound the reduction in <math|log<around*|(|f<around*|(|<x>|)>|)>>
  by

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|log<around*|(|<frac|f<around*|(|<x>+<X><d>|)>|f<around*|(|<x>|)>>|)>>>|<row|<cell|>|<cell|\<leq\>>|<cell|log
    <around*|(|1+<frac|<around*|\<langle\>|<n>f<around*|(|<x>|)>,<X><d>|\<rangle\>>+<frac|1|2><around*|\<langle\>|<d><rsup|\<top\>><X><n><rsup|2>f<around*|(|<x>|)>,<X><d>|\<rangle\>>|f<around*|(|<x>|)>>|)>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|<around*|\<langle\>|<n>f<around*|(|<x>|)>,<X><d>|\<rangle\>>+<frac|1|2><around*|\<langle\>|<d><rsup|\<top\>><X><n><rsup|2>f<around*|(|<x>|)>,<X><d>|\<rangle\>>|f<around*|(|<x>|)>>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<langle\>|<X><n>\<phi\><around*|(|<x>|)>,<d>|\<rangle\>>+<frac|1|2><around*|\<langle\>|<d><rsup|\<top\>><X>\<nabla\><rsup|2>\<phi\><around*|(|<x>|)><X>,<d>|\<rangle\>>+<around*|\<langle\>|<e>,<d>|\<rangle\>>-<around*|\<\|\|\>|<d>|\<\|\|\>><rsup|2>+<frac|\<rho\><d><rsup|\<top\>><X><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>><X><d>|2f<around*|(|<x>|)><rsup|2>>.>>>>
  </eqnarray*>

  Combining the above relations, we deduce that

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<phi\><around*|(|<x>+<X><d>|)>-\<phi\><around*|(|<x>|)>>>|<row|<cell|>|<cell|=>|<cell|\<rho\>log<around*|(|<frac|f<around*|(|<x>+<X><d>|)>|f<around*|(|<x>|)>>|)>+<big|sum><rsub|i=1><rsup|n>log
    x<rsub|i>-<big|sum><rsub|i=1><rsup|n>log
    <around*|(|x<rsub|i>+x<rsub|i>d<rsub|i>|)>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<langle\>|<X><n>\<phi\><around*|(|<x>|)>,<d>|\<rangle\>>+<frac|1|2><around*|\<langle\>|<d><rsup|\<top\>><X>\<nabla\><rsup|2>\<phi\><around*|(|<x>|)><X>,<d>|\<rangle\>>>>|<row|<cell|>|<cell|>|<cell|+<around*|\<langle\>|<e>,<d>|\<rangle\>>-<around*|\<\|\|\>|<d>|\<\|\|\>><rsup|2>-<around*|\<langle\>|<e>,<d>|\<rangle\>>+<frac|\<beta\><rsup|2>|2<around*|(|1-\<beta\>|)>>+<frac|\<rho\><d><rsup|\<top\>><X><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>><X><d>|2f<around*|(|<x>|)><rsup|2>>>>|<row|<cell|>|<cell|\<leq\>>|<cell|-<frac|\<lambda\>|2><around*|\<\|\|\>|<d>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<d>|\<\|\|\>><rsup|2>+<frac|\<beta\><rsup|2>|2<around*|(|1-\<beta\>|)>>+<frac|2\<rho\>\<beta\><rsup|2>|f<around*|(|<x>|)>>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|(|-<frac|\<lambda\>|2>-1+<frac|1|2<around*|(|1-\<beta\>|)>>+<frac|2\<rho\>|n>|)>\<beta\><rsup|2>,>>>>
  </eqnarray*>

  \;

  where the first inequality is by <math|<frac|\<rho\><X><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>><X>|f<around*|(|<x>|)><rsup|2>>\<preceq\><frac|4\<rho\>|n><I>>.

  \;

  <\equation*>
    <matrix|<tformat|<table|<row|<cell|a<rsub|11>>|<cell|>>|<row|<cell|a<rsub|21>>|<cell|a<rsub|22>>>>>><matrix|<tformat|<table|<row|<cell|x<rsub|1>>>|<row|<cell|dx<rsub|2>>>>>>=<matrix|<tformat|<table|<row|<cell|b<rsub|1>>>|<row|<cell|b<rsub|2>>>>>>
  </equation*>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|4|5>>
    <associate|auto-11|<tuple|5|8>>
    <associate|auto-12|<tuple|6|9>>
    <associate|auto-13|<tuple|6.1|9>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|2|2>>
    <associate|auto-4|<tuple|2.1|3>>
    <associate|auto-5|<tuple|2.2|4>>
    <associate|auto-6|<tuple|2.3|4>>
    <associate|auto-7|<tuple|2.4|5>>
    <associate|auto-8|<tuple|3|5>>
    <associate|auto-9|<tuple|1|5>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Solving NETLIB LPs
      in <with|color|<quote|red>|<with|font-series|<quote|bold>|1000>>
      iterations>|<pageref|auto-9>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Dimension-reduced
      Method for Potential Reduction> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Two directions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Potential
      Reduction for LP> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Potential Reduction for HSD
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>Acceleration by negative
      curvature <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|1tab>|2.3<space|2spc>Direct computation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|2.4<space|2spc>Scaled Hessian
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Numerical
      Experiments> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Algorithm
      Design> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Analysis>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>General
      Potential Method> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.5fn>

      <with|par-left|<quote|1tab>|6.1<space|2spc>Second-order Potential
      Reduction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>
    </associate>
  </collection>
</auxiliary>