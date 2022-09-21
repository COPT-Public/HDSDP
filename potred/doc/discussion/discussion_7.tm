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

    <assign|v|<macro|<math-bf|v>>>

    <assign|a|<macro|<math|<math-bf|a>>>>
  </hide-preamble>

  <screens|<\hidden>
    <tit|Dimension-Reduced Interior Point Method>

    \;

    \;

    \;

    \;

    <doc-data|<doc-title|Dimension-reduced Interior Point
    Method>|<doc-author|<author-data|<\author-affiliation>
      \;

      Discussion 7

      \;

      \;

      <date|>
    </author-affiliation>>>>
  </hidden>|<\hidden>
    <tit|Current Progress>

    \;

    \;

    <strong|Theory>

    <\itemize>
      <item><math|\<cal-O\><around*|(|\<varepsilon\><rsup|-1>
      log<around*|(|<frac|1|\<varepsilon\>>|)>|)>> convergence without extra
      assumption

      <item><math|\<cal-O\><around*|(|\<varepsilon\><rsup|-3/4>log<around*|(|<frac|1|\<varepsilon\>>|)>|)>>
      convergence with assumption on the Hessian
    </itemize>

    \;

    <strong|Practice>

    <\itemize>
      <item>The Lanczos solver has been tuned preliminarily

      <item>Now a first order method with accuracy <math|10<rsup|-5>>
      (<math|10<rsup|-8>> if full eigen-decomposition is allowed)
    </itemize>

    \;

    \;

    Now transforming <verbatim|MATLAB> into <verbatim|C> implementation. (In
    one or two weeks)
  </hidden>|<\hidden>
    <tit|Overview of Dimension-reduced Potential Reduction IPM>

    \;

    We consider the model of simplex-constrained QP

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<x>,<y>>>|<cell|<frac|1|2><around*|\<\|\|\>|<A><x>+<math-bf|B><y>|\<\|\|\>><rsup|2>>|<cell|\<backassign\>f<around*|(|<x>,<y>|)>>>|<row|<cell|<text|subject
      to>>|<cell|<e><rsup|\<top\>><x>=1>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>>|<cell|>>>>
    </eqnarray*>

    in the framework of potential reduction

    <\equation*>
      \<varphi\><around*|(|<x>|)>\<assign\>\<rho\> log
      <around*|(|f<around*|(|<x>|)>|)>-<big|sum><rsub|i=1><rsup|n>log
      x<rsub|i>.
    </equation*>

    We update in a dimension-reduced fashion

    <\eqnarray*>
      <tformat|<table|<row|<cell|<d><rsup|k>>|<cell|\<leftarrow\>>|<cell|\<alpha\><rsup|g><g><rsup|k>+\<alpha\><rsup|m><m><rsup|k>>>|<row|<cell|<x><rsup|k>>|<cell|\<leftarrow\>>|<cell|<x><rsup|k>+<d><rsup|k>>>>>
    </eqnarray*>
  </hidden>|<\hidden>
    <tit|Subproblem and Directions>

    \;

    Trust-region controls the accuracy of approximation

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<d>,\<alpha\><rsup|g>,\<alpha\><rsup|m>>>|<cell|<frac|1|2><d><rsup|\<top\>><H><d>+<h><rsup|\<top\>><d>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<around*|\<\|\|\>|<X><rsup|-1><d>|\<\|\|\>>\<leq\>\<Delta\>>|<cell|>>|<row|<cell|>|<cell|<d>=\<alpha\><rsup|g><g>+\<alpha\><rsup|m><m>,>|<cell|>>>>
    </eqnarray*>

    where\ 

    <\eqnarray*>
      <tformat|<table|<row|<cell|<g>\<leftarrow\><n>\<varphi\><around*|(|<x>|)>>|<cell|=>|<cell|<frac|\<rho\><n>f<around*|(|<x>|)>|f<around*|(|<x>|)>>-<X><rsup|-1><e>>>|<row|<cell|<H>\<leftarrow\>\<nabla\><rsup|2>\<varphi\><around*|(|<x>|)>>|<cell|=>|<cell|-<frac|\<rho\><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>>|f<around*|(|<x>|)><rsup|2>>+\<rho\><frac|<A><rsup|\<top\>><A>|f<around*|(|<x>|)>>+<X><rsup|-2>>>|<row|<cell|<m><rsup|k>>|<cell|\<in\>>|<cell|<around*|{|<x><rsup|k>-<x><rsup|k-1>,<v><rsub|min><rsup|<e><rsub|\<bot\>>><around*|(|<H>|)>|}>>>>>
    </eqnarray*>

    \;

    The subproblem solvable using bisection/Newton's method
  </hidden>|<\hidden>
    <tit|Heuristics>

    \;

    <\itemize>
      <item><strong|Radius adjustment>

      The trust radius <math|\<Delta\>\<leq\>1> is adjusted by

      <\equation*>
        \<rho\>=<frac|m<rsup|k><around*|(|\<alpha\>|)>-m<rsup|k><around*|(|<0>|)>|\<varphi\><around*|(|<x><rsup|k>+<d><rsup|\<alpha\>>|)>-\<varphi\><around*|(|<x><rsup|k>|)>>
      </equation*>

      by conventional trust region rules

      \;

      <item><strong|Line-search>

      Line-search procedure is employed to compute

      <\equation*>
        \<alpha\><rsub|max>=max<around*|{|\<gamma\>:<x><rsup|k>+\<gamma\><d><rsup|\<alpha\>>\<geq\><0>|}>
      </equation*>

      and search is performed over the potential function
    </itemize>
  </hidden>|<\hidden>
    <tit|Computational Aspects>

    \;

    <\itemize>
      <item>The algorithm is sensitive to numerical due to the potential
      function

      <item>Simplex is magnified by <math|<e><rsup|\<top\>><x>=n>

      <item>Hessian <math|\<nabla\><rsup|2>\<varphi\><around*|(|<x><rsup|k>|)>>
      has special structures to exploit
    </itemize>

    \;

    <\eqnarray*>
      <tformat|<table|<row|<cell|<around*|\<langle\>|<a>,\<nabla\><rsup|2>\<varphi\><around*|(|<x><rsup|k>|)><a>|\<rangle\>>>|<cell|=>|<cell|<around*|\<langle\>|<a>,-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><a>|f<around*|(|<x><rsup|k>|)><rsup|2>>|\<rangle\>>+<frac|<around*|\<\|\|\>|<A><a>|\<\|\|\>><rsup|2>|f<around*|(|<x><rsup|k>|)>>+<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><a>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|-\<rho\><around*|(|<frac|<n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><a>|f<around*|(|<x><rsup|k>|)>>|)><rsup|2>+<frac|<around*|\<\|\|\>|<A><a>|\<\|\|\>><rsup|2>|f<around*|(|<x><rsup|k>|)>>+<around*|\<\|\|\>|<around*|(|<X><rsup|k>|)><rsup|-1><a>|\<\|\|\>><rsup|2>>>|<row|<cell|<around*|\<langle\>|<a>,\<nabla\><rsup|2>\<varphi\><around*|(|<x><rsup|k>|)><b>|\<rangle\>>>|<cell|=>|<cell|<around*|\<langle\>|<a>,-<frac|\<rho\><n>f<around*|(|<x><rsup|k>|)><n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><b>|f<around*|(|<x><rsup|k>|)><rsup|2>>|\<rangle\>>+<frac|<around*|\<langle\>|<A><a>,<A><b>|\<rangle\>>|f<around*|(|<x><rsup|k>|)>>+<around*|\<langle\>|<a>,<around*|(|<X><rsup|k>|)><rsup|-2><b>|\<rangle\>>>>|<row|<cell|>|<cell|=>|<cell|-\<rho\><around*|(|<frac|<n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><a>|f<around*|(|<x><rsup|k>|)>>|)><around*|(|<frac|<n>f<around*|(|<x><rsup|k>|)><rsup|\<top\>><b>|f<around*|(|<x><rsup|k>|)>>|)>+<frac|<around*|\<langle\>|<A><a>,<A><b>|\<rangle\>>|f<around*|(|<x><rsup|k>|)>>+<around*|\<langle\>|<a>,<around*|(|<X><rsup|k>|)><rsup|-2><b>|\<rangle\>>.>>>>
    </eqnarray*>

    \;

    The most important and ill-conditioned step: negative curvature
    computation
  </hidden>|<\hidden>
    <tit|Negative Curvature Computation>

    Negative curvature is quite efficient in accelerating convergence but is
    hard to compute.

    <\itemize>
      <item>Direct method (cheaper but less stable)

      <\eqnarray*>
        <tformat|<table|<row|<cell|min<rsub|<around*|\<\|\|\>|<v>|\<\|\|\>>=1>>|<cell|<v><rsup|\<top\>><tiny|<around*|{|<frac|2\<rho\><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><A><rsup|\<top\>><A><u><u><rsup|\<top\>><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4>>+<small|<tiny|<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X><rsup|-2>>>>>>>>|}>><v>>|<cell|>>|<row|<cell|<text|subject
        to>>|<cell|<e><rsup|\<top\>><v><rsub|<x>>=0.>|<cell|>>>>
      </eqnarray*>

      <item>Scaled Hessian (expensive but more stable)

      <\eqnarray*>
        <tformat|<table|<row|<cell|min<rsub|<around*|\<\|\|\>|<v>|\<\|\|\>>=1>>|<cell|<v><rsup|\<top\>><small|<tiny|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><small|<around*|{|<tiny|<frac|2\<rho\><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|2>>-<frac|4\<rho\><A><rsup|\<top\>><A><u><u><rsup|\<top\>><A><rsup|\<top\>><A>|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|4>>>+<small|<tiny|<matrix|<tformat|<table|<row|<cell|<0><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X><rsup|-2>>>>>>>>|}>>>><tiny|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>>><v>>|<cell|>>|<row|<cell|<text|subject
        to>>|<cell|<x><rsup|\<top\>><v><rsub|<x>>=0.>|<cell|>>>>
      </eqnarray*>

      <item>Reduced support (cheap and more stable)

      <\equation*>
        <v><rsub|k>=0,k\<leq\>\<delta\>
      </equation*>
    </itemize>

    Last, a customized Lanczos procedure is employed to solve for <math|<v>>.
  </hidden>|<\shown>
    <tit|Customized Lanczos>

    \;

    A customized Lanczos procedure

    <\eqnarray*>
      <tformat|<table|<row|<cell|<x><rprime|'>>|<cell|\<leftarrow\>>|<cell|<frac|<x>|<around*|\<\|\|\>|<x>|\<\|\|\>>>>>|<row|<cell|<v>>|<cell|\<leftarrow\>>|<cell|<matrix|<tformat|<table|<row|<cell|<v><rsub|<y>>>>|<row|<cell|<v><rsub|<x>>-<around*|(|<x><rprime|'><rsup|\<top\>><v><rsub|<x>>|)><x><rprime|'>>>>>>>>|<row|<cell|<u><rsub|1>>|<cell|\<leftarrow\>>|<cell|<matrix|<tformat|<table|<row|<cell|<0>>>|<row|<cell|<v><rsub|<x>>-<around*|(|<x><rprime|'><rsup|\<top\>><v><rsub|<x>>|)><x><rprime|'>>>>>>>>|<row|<cell|<u><rsub|2>>|<cell|\<leftarrow\>>|<cell|<small|<matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><rprime|'><x><rprime|'><rsup|\<top\>>>>>>>><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><A><rsup|\<top\>><A><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><v>>>|<row|<cell|<u><rsub|3>>|<cell|\<leftarrow\>>|<cell|<g><rsup|\<top\>><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><v><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<I><rsub|n>-<x><rprime|'><x><rprime|'><rsup|\<top\>>>>>>><matrix|<tformat|<table|<row|<cell|<I><rsub|m>>|<cell|>>|<row|<cell|>|<cell|<X>>>>>><g>>>|<row|<cell|<text|Assemble>>|<cell|\<leftarrow\>>|<cell|f<rsup|2><u><rsub|1>+\<rho\>f<u><rsub|2>-\<rho\><u><rsub|3>>>>>
    </eqnarray*>

    Re-orthogonalization is employed to improve convergence.

    Early stop adapted from SDPT3 reduces number of iterations.

    \;
  </shown>|<\hidden>
    <tit|Numerical Experiments>

    <big-table|<tiny|<block|<tformat|<cwith|1|19|1|8|cell-halign|c>|<table|<row|<cell|Instance>|<cell|Pinf.>|<cell|Dinf.>|<cell|Compl.>|<cell|Instance>|<cell|Pinf.>|<cell|Dinf.>|<cell|Compl.>>|<row|<cell|ADLITTLE>|<cell|3.326e-07>|<cell|5.581e-07>|<cell|7.617e-06>|<cell|SC105>|<cell|3.096e-06>|<cell|2.351e-06>|<cell|1.324e-05>>|<row|<cell|AFIRO>|<cell|3.597e-06>|<cell|3.487e-06>|<cell|1.467e-05>|<cell|SC205>|<cell|2.591e-05>|<cell|2.077e-05>|<cell|9.955e-05>>|<row|<cell|AGG2>|<cell|4.608e-04>|<cell|9.277e-05>|<cell|1.121e-03>|<cell|SC50A>|<cell|1.062e-05>|<cell|5.956e-06>|<cell|4.219e-05>>|<row|<cell|BANDM>|<cell|3.869e-06>|<cell|2.935e-06>|<cell|2.232e-05>|<cell|SC50B>|<cell|1.562e-06>|<cell|1.259e-06>|<cell|7.846e-06>>|<row|<cell|BEACONFD>|<cell|9.495e-07>|<cell|1.642e-06>|<cell|1.753e-05>|<cell|SCAGR25>|<cell|7.654e-06>|<cell|4.435e-06>|<cell|1.075e-04>>|<row|<cell|BLEND>|<cell|1.479e-06>|<cell|2.673e-06>|<cell|9.001e-06>|<cell|SCAGR7>|<cell|9.592e-07>|<cell|3.253e-07>|<cell|7.192e-06>>|<row|<cell|BOEING2>|<cell|1.841e-05>|<cell|1.571e-06>|<cell|3.407e-05>|<cell|SCFXM1>|<cell|3.558e-05>|<cell|2.766e-05>|<cell|7.141e-05>>|<row|<cell|BORE3D>|<cell|2.493e-05>|<cell|7.667e-05>|<cell|1.895e-04>|<cell|SCORPION>|<cell|1.174e-06>|<cell|1.328e-06>|<cell|1.249e-05>>|<row|<cell|BRANDY>|<cell|3.477e-05>|<cell|1.398e-05>|<cell|7.888e-05>|<cell|SCTAP1>|<cell|9.530e-07>|<cell|1.702e-06>|<cell|9.649e-06>>|<row|<cell|FINNIS>|<cell|3.468e-05>|<cell|3.486e-05>|<cell|3.622e-04>|<cell|SEBA>|<cell|2.459e-08>|<cell|1.014e-07>|<cell|5.075e-07>>|<row|<cell|FORPLAN>|<cell|3.323e-05>|<cell|1.922e-05>|<cell|1.717e-04>|<cell|SHARE1B>|<cell|2.614e-05>|<cell|2.470e-05>|<cell|1.034e-04>>|<row|<cell|GFRD-PNC>|<cell|4.032e-04>|<cell|1.004e-04>|<cell|3.425e-04>|<cell|SHARE2B>|<cell|8.259e-04>|<cell|1.968e-04>|<cell|5.570e-04>>|<row|<cell|GROW7>|<cell|3.069e-04>|<cell|4.716e-05>|<cell|6.817e-04>|<cell|STAIR>|<cell|4.136e-04>|<cell|7.065e-06>|<cell|2.526e-05>>|<row|<cell|ISRAEL>|<cell|2.326e-03>|<cell|1.776e-04>|<cell|9.293e-04>|<cell|STANDATA>|<cell|6.528e-06>|<cell|9.310e-06>|<cell|1.734e-04>>|<row|<cell|KB2>|<cell|4.759e-06>|<cell|3.436e-05>|<cell|4.135e-06>|<cell|STANDGUB>|<cell|1.175e-04>|<cell|3.768e-05>|<cell|6.312e-04>>|<row|<cell|LOTFI>|<cell|1.445e-06>|<cell|1.496e-06>|<cell|3.365e-05>|<cell|STOCFOR1>|<cell|2.798e-05>|<cell|2.346e-05>|<cell|1.125e-04>>|<row|<cell|MODSZK1>|<cell|5.552e-05>|<cell|5.675e-04>|<cell|2.214e-03>|<cell|VTP-BASE>|<cell|1.414e-05>|<cell|2.472e-06>|<cell|2.246e-05>>|<row|<cell|RECIPELP>|<cell|7.721e-06>|<cell|9.249e-06>|<cell|1.676e-05>|<cell|>|<cell|>|<cell|>|<cell|>>>>>>|Computation
    without full eigen-decomposition (1000 iterations)>

    The method solves LPs to

    <\itemize>
      <item><math|10<rsup|-8>\<sim\>10<rsup|-10>> if full-eigen decomposition
      is allowed

      <item><math|10<rsup|-4>\<sim\>10<rsup|-6>> if Lanczos is employed
    </itemize>

    Now optimizing the code and enhancing its implementation
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

<\auxiliary>
  <\collection>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Computation without
      full eigen-decomposition (1000 iterations)>|<pageref|auto-1>>
    </associate>
  </collection>
</auxiliary>