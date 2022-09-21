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

  <screens|<\hidden>
    <tit|Dimension-Reduced Interior Point Method>

    \;

    \;

    \;

    \;

    <doc-data|<doc-title|Dimension-reduced Interior Point
    Method>|<doc-author|<author-data|<\author-affiliation>
      Discussion 2

      \;

      <date|>
    </author-affiliation>>>>
  </hidden>|<\hidden>
    <tit|Substituting Conjugate Gradient>

    \;

    More directions are chosen for CG.\ 

    <\equation*>
      <x><rsup|k>+<d><rsup|k>=<x><rsup|k>+\<alpha\><rsup|g><around*|(|<A><x><rsup|k>-<b>|)>+\<alpha\><rsup|m><m><rsup|k>
    </equation*>

    If <math|<m><rsup|k>\<approx\><x><rsup|\<ast\>>-<x><rsup|k>>, then
    <math|\<alpha\><rsup|m>=1> solves the system.

    \;

    Intuitively

    <\itemize>
      <item>We choose <math|<wide|<x>|^>> heuristically to be close to
      <math|<x>>

      <item>Then take <math|<m><rsup|k>=<wide|<x>|^>-<x><rsup|k>>

      e.g. <math|<wide|<x>|^>=diag<around*|(|<A>|)><rsup|-1><b>>,
      <math|<wide|<x>|^>=tridiag<around*|(|<A>|)><rsup|-1><b>>

      <item>Return to CG after some steps
    </itemize>

    \;

    A little better than CG wehen <math|<A>> has special structure (e.g.,
    dominant diagonal).
  </hidden>|<\hidden>
    <tit|Potential Reduction Method>

    Still solving

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<x>>>|<cell|<frac|1|2><around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>|<cell|\<backassign\>f<around*|(|<x>|)>>>|<row|<cell|<text|subject
      to>>|<cell|<e><rsup|\<top\>><x>=1>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>>|<cell|>>>>
    </eqnarray*>

    with

    <\eqnarray*>
      <tformat|<table|<row|<cell|<d><rsup|k>>|<cell|\<leftarrow\>>|<cell|\<alpha\><rsup|g><P><rsub|\<Delta\>><around*|[|\<nabla\>\<varphi\><around*|(|<x><rsup|k>|)>|]>+\<alpha\><rsup|m><around*|(|<x><rsup|k>-<x><rsup|k-1>|)>>>|<row|<cell|<x><rsup|k>>|<cell|\<leftarrow\>>|<cell|<x><rsup|k>+<d><rsup|k>>>>>
    </eqnarray*>

    where <math|<P><rsub|\<Delta\>><around*|[|\<cdummy\>|]>> is the
    orthogonal projection onto <math|<e><rsup|\<top\>><x>=0>.
    <math|\<alpha\><rsup|g>,\<alpha\><rsup|d>> come from the following model

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<d>,\<alpha\><rsup|g>,\<alpha\><rsup|m>>>|<cell|<frac|1|2><d><rsup|\<top\>><H><d>+<h><rsup|\<top\>><d>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<around*|\<\|\|\>|<X><rsup|-1><d>|\<\|\|\>>\<leq\>\<Delta\>>|<cell|>>|<row|<cell|>|<cell|<d>=\<alpha\><rsup|g><g><rsup|k>+\<alpha\><rsup|m><m><rsup|k>>|<cell|>>>>
    </eqnarray*>
  </hidden>|<\hidden>
    <tit|Implementation>

    \;

    <\itemize>
      <item>Hessian vector product is relatively cheap

      <math|\<nabla\><rsup|2><rsub|<x>,<x>>\<varphi\><around*|(|<x><rsup|k>|)>=-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>>|f<around*|(|<x><rsup|k>|)><rsup|2>>+\<rho\><frac|<A><rsup|\<top\>><A>|f<around*|(|<x><rsup|k>|)>>+<around*|(|<X><rsup|k>|)><rsup|-2>>

      <item>Trust radius <math|\<beta\>> is adjusted by\ 

      <\equation*>
        <frac|m<rsup|\<varphi\>><around*|(|\<alpha\>|)>-m<rsup|\<varphi\>><around*|(|0|)>|\<varphi\><around*|(|<x><rsup|k>+<d><rsup|k>|)>-\<varphi\><around*|(|<x><rsup|k>|)>>
      </equation*>

      and <math|\<beta\>\<leq\>1> to ensure feasibility

      <item>Scaling is imposed to enhance stability

      <\eqnarray*>
        <tformat|<table|<row|<cell|min<rsub|<d>,\<alpha\><rsup|g>,\<alpha\><rsup|m>>>|<cell|<frac|1|2><d><rsup|\<top\>><X><rsup|-1><H><X><rsup|-1><d>+<around*|(|<X><rsup|-1><h>|)><rsup|\<top\>><d>>|<cell|>>|<row|<cell|<text|subject
        to>>|<cell|<around*|\<\|\|\>|<d>|\<\|\|\>>\<leq\>\<Delta\>>|<cell|>>|<row|<cell|>|<cell|<d>=\<alpha\><rsup|g><X><g><rsup|k>+\<alpha\><rsup|m><X><m><rsup|k>>|<cell|>>>>
      </eqnarray*>

      <item>Potential function reduces much faster and more stably
    </itemize>
  </hidden>|<\shown>
    <tit|Performance>

    <big-figure|<image|file:///Users/gaowenzhi/Desktop/potred/matlab/compare.pdf||400px||><space|1em><image|file:///Users/gaowenzhi/Desktop/potred/matlab/iter.pdf||400px||>|Left:
    <math|f<around*|(|<x><rsup|k>|)>><space|1em>Right: Momentum might
    accelerate convergence>

    \;

    <\itemize>
      <item>Reaching <verbatim|1e-06> accuracy is easy on synthetic data

      Around 20 iterations are required

      <item>Momentum plays an important role as
      <math|\<alpha\><rsup|m>\<gtr\>\<alpha\><rsup|g>>.
    </itemize>
  </shown>>
</body>

<\initial>
  <\collection>
    <associate|math-color|black>
    <associate|page-height|auto>
    <associate|page-medium|paper>
    <associate|page-type|4:3>
    <associate|page-width|auto>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Left:
      <with|color|<quote|black>|font-family|<quote|rm>|<with|mode|<quote|math>|f<around*|(|<with|color|<quote|black>|font-family|<quote|rm>|<with|mode|<quote|math>|<rigid|<with|mode|<quote|text>|<with|font-family|<quote|rm>|font-series|<quote|bold>|font-shape|<quote|right>|x>>>>><rsup|k>|)>>><space|1em>Right:
      Momentum might accelerate convergence>|<pageref|auto-1>>
    </associate>
  </collection>
</auxiliary>