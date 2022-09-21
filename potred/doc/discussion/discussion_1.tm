<TeXmacs|2.1>

<style|beamer>

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
  </hide-preamble>

  <screens|<\shown>
    <tit|Dimension-Reduced Interior Point Method>

    \;

    \;

    \;

    \;

    <doc-data|<doc-title|Dimension-reduced Interior Point
    Method>|<doc-author|<author-data|<\author-affiliation>
      Discussion 1

      \;

      <date|>
    </author-affiliation>>>>
  </shown>|<\hidden>
    <tit|Substituting Conjugate Gradient>

    One direction applies dimension-reduced method to replace CG to solve
    <math|<A><x>=<b>>.

    <\itemize>
      <item>The iterations are almost identical when
      <math|<around*|\<\|\|\>|<A><x><rsup|k>-<b>|\<\|\|\>>\<geq\>10<rsup|-15>>.

      <center|<image|file:///Users/gaowenzhi/Desktop/potred/matlab/conv.pdf|400px|||>>

      <item>Sometimes better than CG when the problem is ill-conditioned

      <center|<image|file:///Users/gaowenzhi/Desktop/potred/matlab/\<#622A\>\<#5C4F\>2022-08-04
      11.45.02.png|600px|||>>
    </itemize>
  </hidden>|<\hidden>
    <tit|Potential Reduction Method>

    We start from the simple case of simplex-constrained QP

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<x>>>|<cell|<frac|1|2><around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>|<cell|\<backassign\>f<around*|(|<x>|)>>>|<row|<cell|<text|subject
      to>>|<cell|<e><rsup|\<top\>><x>=1>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>>|<cell|>>>>
    </eqnarray*>

    and use the potential function

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<varphi\><around*|(|<x>|)>>|<cell|\<assign\>>|<cell|\<rho\>
      log <around*|(|f<around*|(|<x>|)>|)>-<small|<big|sum><rsub|i=1><rsup|n>>log
      x<rsub|i>>>|<row|<cell|<n>\<varphi\><around*|(|<x>|)>>|<cell|=>|<cell|<frac|\<rho\><n>f<around*|(|<x>|)>|f<around*|(|<x>|)>>-<X><rsup|-1><e>.>>>>
    </eqnarray*>

    Potential reduction solves for <math|\<Delta\>\<assign\><x><rsup|k+1>-<x><rsup|k>>
    at each iteration.

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|\<Delta\>>>|<cell|<around*|\<langle\>|<n>\<varphi\><around*|(|<x><rsup|k>|)>,\<Delta\>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<e><rsup|\<top\>>\<Delta\>=0>|<cell|>>|<row|<cell|>|<cell|<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1>\<Delta\>|\<\|\|\>>\<leq\>\<beta\>,>|<cell|>>>>
    </eqnarray*>
  </hidden>|<\hidden>
    <tit|Introducing Momentum>

    We follow the gradient projection framework and define
    <math|<P><rsub|\<Delta\>><around*|[|<x>|]>\<assign\><around*|(|<I>-<frac|<e><e><rsup|\<top\>>|<around*|\<\|\|\>|<e>|\<\|\|\>><rsup|2>>|)><x>>.
    Then we consider

    <\eqnarray*>
      <tformat|<table|<row|<cell|<x><rsup|k+<frac|1|2>>>|<cell|\<leftarrow\>>|<cell|<x><rsup|k>+\<alpha\><rsup|g>\<nabla\>\<varphi\><around*|(|<x><rsup|k>|)>+\<alpha\><rsup|m><around*|(|<x><rsup|k>-<x><rsup|k-1>|)>>>|<row|<cell|<x><rsup|k+1>>|<cell|\<leftarrow\>>|<cell|<P><rsub|\<Delta\>><around*|[|<x><rsup|k+<frac|1|2>>|]>,>>>>
    </eqnarray*>

    where <math|\<alpha\><rsup|g>,\<alpha\><rsup|d>> come through

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|\<alpha\>>>|<cell|<frac|1|2>\<alpha\><rsup|\<top\>><H>\<alpha\>+<h><rsup|\<top\>>\<alpha\>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<around*|\<\|\|\>|\<alpha\><rsup|g><around*|(|<X><rsup|k>|)><rsup|-1><g><rsup|k>+\<alpha\><rsup|m><around*|(|<X><rsup|k>|)><rsup|-1><m><rsup|k>|\<\|\|\>>\<leq\>\<beta\>,>|<cell|>>>>
    </eqnarray*>

    <\equation*>
      <H>\<assign\><matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<g><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<langle\>|<g><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>|<row|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>>>>,<space|1em><h>=<matrix|<tformat|<table|<row|<cell|<around*|\<\|\|\>|<g><rsup|k>|\<\|\|\>><rsup|2>>>|<row|<cell|<around*|\<langle\>|<g><rsup|k>,<m><rsup|k>|\<rangle\>>>>>>>
    </equation*>

    and\ 

    <\equation*>
      \<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)>=-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>>|f<around*|(|<x><rsup|k>|)><rsup|2>>+<around*|(|<X><rsup|k>|)><rsup|-2>.
    </equation*>
  </hidden>|<\hidden>
    <tit|Preliminary observations>

    <\itemize>
      <item>The adaptive trust radius replaces original
      <math|\<beta\><rsup|k>=<frac|1|2+<frac|\<rho\>\<gamma\>|f<around*|(|<x><rsup|k>|)>>>>

      <item>May not be strictly decreasing for some steps due to large
      <math|\<beta\><rsup|k>>

      <item>Some issues are observed
    </itemize>

    <big-figure|<image|file:///Users/gaowenzhi/Desktop/potred/matlab/resi.pdf|500px|||>|<math|<around*|\<\|\|\>|<A><rsup|k><x>|\<\|\|\>><rsup|2>\<sim\>k>>
  </hidden>|<\hidden>
    <tit|Issues with potential reduction>

    \;

    Recall that

    <\equation*>
      <H>=<matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<g><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<langle\>|<g><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>|<row|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><g><rsup|k>|\<rangle\>>>|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>>>>
    </equation*>

    and\ 

    <\equation*>
      \<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)>=-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>>|f<around*|(|<x><rsup|k>|)><rsup|2>>+<around*|(|<X><rsup|k>|)><rsup|-2>
    </equation*>

    <\equation*>
      <g><rsup|k>=<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)>|f<around*|(|<x><rsup|k>|)>>-<around*|(|<X><rsup|k>|)><rsup|-1><e>.
    </equation*>

    \;

    <\itemize>
      <item>The matrix <math|<H>> is almost always ill-conditioned since
      <math|<around*|\||<around*|\<langle\>|<g><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><g><rsup|k>|\<rangle\>>|\|>\<gg\>0>

      this drives <math|\<alpha\><rsup|g>> and <math|\<alpha\><rsup|d>>
      imbalanced and subproblem hard to solve

      <item>Often only <math|\<alpha\><rsup|g>> works and this looks like
      steepest descent with adaptive <math|\<beta\><rsup|k>>.
    </itemize>
  </hidden>|<\hidden>
    <tit|Future directions>

    \;

    One possible direction is to consider the projected gradient in the model
    problem

    <\equation*>
      <H><rsub|><rsup|<P>>=<matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<P><rsub|\<Delta\>><around*|[|<g><rsup|k>|]>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><P><rsub|\<Delta\>><around*|[|<g><rsup|k>|]>|\<rangle\>>>|<cell|<around*|\<langle\>|<P><rsub|\<Delta\>><around*|[|<g><rsup|k>|]>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>|<row|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><P><rsub|\<Delta\>><around*|[|<g><rsup|k>|]>|\<rangle\>>>|<cell|<around*|\<langle\>|<m><rsup|k>,\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)><m><rsup|k>|\<rangle\>>>>>>>.
    </equation*>

    Or look for other proper ways to replace steepest descent in potential
    reduction.

    \;

    The other direction is to consider a more general formulation. e.g., the
    dual potential function

    <\equation*>
      \<rho\>log<around*|(|z-<b><rsup|\<top\>><y>|)>-<big|sum><rsub|i=1><rsup|n>log<around*|(|<c><rsub|i>-<math-bf|a><rsub|i><rsup|\<top\>><y>|)>.
    </equation*>

    \;
  </hidden>>
</body>

<\initial>
  <\collection>
    <associate|math-color|black>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
  </collection>
</references>