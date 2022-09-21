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

      Discussion 4. Part 1

      \;

      \;

      <date|>
    </author-affiliation>>>>
  </hidden>|<\hidden>
    <tit|Solving LPs>

    Simplex constrained QP formulation

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<y>,<u>=<around*|(|<x>,<s>,\<kappa\>,\<tau\>|)>>>|<cell|f<around*|(|<u>|)>\<assign\><frac|1|2><around*|\<\|\|\>|<wide|<A>|^><around*|(|<y>;<u>|)>|\<\|\|\>><rsup|2>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<e><rsup|\<top\>><u>=1>|<cell|>>|<row|<cell|>|<cell|<u>\<geq\><0>.>|<cell|>>>>
    </eqnarray*>

    Using the potential function

    <\equation*>
      \<varphi\><around*|(|<u>|)>\<assign\>\<rho\> log
      <around*|(|f<around*|(|<u>|)>|)>-B<around*|(|<x>|)>-B<around*|(|<s>|)>-log
      \<kappa\>-log \<tau\>
    </equation*>

    \;

    <\itemize>
      <item><math|<y>> is treated unconstrained
    </itemize>

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<d>,\<alpha\><rsup|g>,\<alpha\><rsup|m>>>|<cell|<frac|1|2><d><rsup|\<top\>><H><d>+<h><rsup|\<top\>><d>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<around*|\<\|\|\>|<math-bf|U><rsup|-1><u>|\<\|\|\>><rsup|2>\<leq\>\<beta\>\<leq\>1>|<cell|>>|<row|<cell|>|<cell|<d>=\<alpha\><rsup|g><g><rsup|k>+\<alpha\><rsup|m><m><rsup|k>.>|<cell|>>>>
    </eqnarray*>
  </hidden>|<\hidden>
    <tit|Preliminary results \U Gradient + Momentum>

    The new method

    <\itemize>
      <item>is able to solving NETLIB LPs to
      <math|10<rsup|-3>\<sim\>10<rsup|-4>> relative accuracy
      (<with|color|red|20000> iterations)

      <item>slows down in reducing potential half way
    </itemize>

    <big-figure|<image|file:///Users/gaowenzhi/Desktop/potred/doc/fig1.jpg|1000px|||>|Reduction
    in potential>

    Main observation: steps may be trapped at local solutions.

    So we use one more direction to <with|color|red|escape> it.
  </hidden>|<\hidden>
    <tit|Adding Negative Curvatire>

    Note that

    <\itemize>
      <item><math|<math-bf|v><rsub|min>> that corresponds to
      <math|\<lambda\><rsub|min><around*|(|<n><rsub|<x><x>>\<varphi\>|)>> is
      an ideal direction to escape a local solution.
    </itemize>

    We only need to consider <math|<tabular|<tformat|<table|<row|<cell|min<rsub|<v>\<neq\><0>>>|<cell|<frac|<around*|\<langle\>|<v>,<n><rsub|<x><x>>\<varphi\><v>|\<rangle\>>|<around*|\<\|\|\>|<v>|\<\|\|\>><rsup|2>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsup|\<top\>><v>=0,>|<cell|>>>>>> or simply

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<v>\<neq\><0>>>|<cell|<frac|<around*|\<langle\>|<v>,<P><rsub|\<Delta\>><n><rsub|<x><x>>\<varphi\><P><rsub|\<Delta\>><v>|\<rangle\>>|<around*|\<\|\|\>|<v>|\<\|\|\>><rsup|2>>>|<cell|>>>>
    </eqnarray*>

    \;

    <\itemize>
      <item>Finding the miminum (negative) eigenvalue of
      <math|<P><rsub|\<Delta\>><n><rsub|<x><x>>\<varphi\><P><rsub|\<Delta\>>>

      <item>Efficiently solvable using Lanczos exploiting\ 

      <\equation*>
        <n><rsub|<x><x>>\<varphi\><around*|(|<x>|)>=\<rho\><around*|(|-<frac|<n>f<around*|(|<x>|)><n>f<around*|(|<x>|)>|f<around*|(|<x>|)><rsup|2>><rsup|\<top\>>+<frac|<A><rsup|\<top\>><A>|f<around*|(|<x>|)>>|)>+<X><rsup|-2>
      </equation*>
    </itemize>

    and structure of <math|<A>>.
  </hidden>|<\hidden>
    <tit|Efficiency of Negative Curvature>

    \;

    In practice we can add curvature

    <\itemize>
      <item>if <math|\<varphi\><around*|(|<x><rsup|k>|)>-\<varphi\><around*|(|<x><rsup|k+1>|)>\<leq\>\<varepsilon\>>

      <item>every <math|K> iterations
    </itemize>

    <center|<big-figure|<image|file:///Users/gaowenzhi/Desktop/potred/doc/IMG_3825.PNG|800px|||>|Adding
    negative curvature>>
  </hidden>|<\shown>
    <tit|Efficiency of Negative Curvature>

    <\small>
      Original: solving NETLIB LPs to <math|<with|color|red|10<rsup|-3>\<sim\>10<rsup|-4>>>
      relative accuracy in <with|color|red|20000> iterations

      Now: solving NETLIB LPs to <math|<with|color|red|10<rsup|-5>\<sim\>10<rsup|-8>>>
      relative accuracy in <with|color|red|1000> iterations
    </small>

    <big-table|<small|<block|<tformat|<cwith|2|13|1|4|cell-valign|c>|<cwith|1|-1|1|-1|cell-halign|c>|<cwith|1|13|5|8|cell-halign|c>|<cwith|2|7|6|8|cell-valign|c>|<cwith|2|7|5|5|cell-valign|c>|<cwith|14|18|1|4|cell-halign|c>|<cwith|14|18|1|4|cell-halign|c>|<cwith|14|18|2|4|cell-valign|c>|<cwith|14|18|1|1|cell-valign|c>|<cwith|2|3|5|8|cell-halign|c>|<cwith|2|3|5|8|cell-halign|c>|<cwith|2|2|6|8|cell-valign|c>|<cwith|2|2|5|5|cell-valign|c>|<cwith|4|13|5|5|cell-halign|c>|<cwith|4|8|5|5|cell-halign|c>|<cwith|4|13|5|5|cell-valign|c>|<cwith|4|13|6|8|cell-valign|c>|<cwith|14|21|5|5|cell-rborder|0ln>|<cwith|14|21|5|5|cell-valign|c>|<cwith|14|18|6|8|cell-valign|c>|<cwith|19|20|1|4|cell-valign|c>|<cwith|19|19|5|8|cell-halign|c>|<cwith|19|19|5|8|cell-valign|c>|<cwith|20|20|1|4|cell-valign|c>|<cwith|20|20|5|8|cell-valign|c>|<cwith|1|-1|1|-1|cell-tborder|1ln>|<cwith|1|-1|1|-1|cell-bborder|1ln>|<cwith|1|-1|1|-1|cell-lborder|0ln>|<cwith|1|-1|1|-1|cell-rborder|0ln>|<cwith|1|-1|5|5|cell-tborder|1ln>|<cwith|1|-1|5|5|cell-bborder|1ln>|<cwith|1|-1|5|5|cell-lborder|1ln>|<cwith|1|-1|5|5|cell-rborder|1ln>|<cwith|1|-1|4|4|cell-rborder|1ln>|<cwith|1|-1|6|6|cell-lborder|1ln>|<table|<row|<cell|Problem>|<cell|PInfeas>|<cell|DInfeas.>|<cell|Compl.>|<cell|Problem>|<cell|PInfeas>|<cell|DInfeas.>|<cell|Compl.>>|<row|<cell|ADLITTLE>|<cell|1.347e-10>|<cell|2.308e-10>|<cell|2.960e-09>|<cell|KB2>|<cell|5.455e-11>|<cell|6.417e-10>|<cell|7.562e-11>>|<row|<cell|AFIRO>|<cell|7.641e-11>|<cell|7.375e-11>|<cell|3.130e-10>|<cell|LOTFI>|<cell|2.164e-09>|<cell|4.155e-09>|<cell|8.663e-08>>|<row|<cell|AGG2>|<cell|3.374e-08>|<cell|4.859e-08>|<cell|6.286e-07>|<cell|MODSZK1>|<cell|1.527e-06>|<cell|5.415e-05>|<cell|2.597e-04>>|<row|<cell|AGG3>|<cell|2.248e-05>|<cell|1.151e-06>|<cell|1.518e-05>|<cell|RECIPELP>|<cell|5.868e-08>|<cell|6.300e-08>|<cell|1.285e-07>>|<row|<cell|BANDM>|<cell|2.444e-09>|<cell|4.886e-09>|<cell|3.769e-08>|<cell|SC105>|<cell|7.315e-11>|<cell|5.970e-11>|<cell|2.435e-10>>|<row|<cell|BEACONFD>|<cell|5.765e-12>|<cell|9.853e-12>|<cell|1.022e-10>|<cell|SC205>|<cell|6.392e-11>|<cell|5.710e-11>|<cell|2.650e-10>>|<row|<cell|BLEND>|<cell|2.018e-10>|<cell|3.729e-10>|<cell|1.179e-09>|<cell|SC50A>|<cell|1.078e-05>|<cell|6.098e-06>|<cell|4.279e-05>>|<row|<cell|BOEING2>|<cell|1.144e-07>|<cell|1.110e-08>|<cell|2.307e-07>|<cell|SC50B>|<cell|4.647e-11>|<cell|3.269e-11>|<cell|1.747e-10>>|<row|<cell|BORE3D>|<cell|2.389e-08>|<cell|5.013e-08>|<cell|1.165e-07>|<cell|SCAGR25>|<cell|1.048e-07>|<cell|5.298e-08>|<cell|1.289e-06>>|<row|<cell|BRANDY>|<cell|2.702e-05>|<cell|7.818e-06>|<cell|1.849e-05>|<cell|SCAGR7>|<cell|1.087e-07>|<cell|1.173e-08>|<cell|2.601e-07>>|<row|<cell|CAPRI>|<cell|7.575e-05>|<cell|4.488e-05>|<cell|4.880e-05>|<cell|SCFXM1>|<cell|4.323e-06>|<cell|5.244e-06>|<cell|8.681e-06>>|<row|<cell|E226>|<cell|2.656e-06>|<cell|4.742e-06>|<cell|2.512e-05>|<cell|SCORPION>|<cell|1.674e-09>|<cell|1.892e-09>|<cell|1.737e-08>>|<row|<cell|FINNIS>|<cell|8.577e-07>|<cell|8.367e-07>|<cell|1.001e-05>|<cell|SCTAP1>|<cell|5.567e-07>|<cell|8.430e-07>|<cell|5.081e-06>>|<row|<cell|FORPLAN>|<cell|5.874e-07>|<cell|2.084e-07>|<cell|4.979e-06>|<cell|SEBA>|<cell|2.919e-11>|<cell|5.729e-11>|<cell|1.448e-10>>|<row|<cell|GFRD-PNC>|<cell|4.558e-05>|<cell|1.052e-05>|<cell|4.363e-05>|<cell|SHARE1B>|<cell|3.367e-07>|<cell|1.339e-06>|<cell|3.578e-06>>|<row|<cell|GROW7>|<cell|1.276e-04>|<cell|4.906e-06>|<cell|1.024e-04>|<cell|SHARE2B>|<cell|2.142e-04>|<cell|2.014e-05>|<cell|6.146e-05>>|<row|<cell|ISRAEL>|<cell|1.422e-06>|<cell|1.336e-06>|<cell|1.404e-05>|<cell|STAIR>|<cell|5.549e-04>|<cell|8.566e-06>|<cell|2.861e-05>>|<row|<cell|STANDATA>|<cell|5.645e-08>|<cell|2.735e-07>|<cell|5.130e-06>|<cell|STANDGUB>|<cell|2.934e-08>|<cell|1.467e-07>|<cell|2.753e-06>>|<row|<cell|STOCFOR1>|<cell|6.633e-09>|<cell|9.701e-09>|<cell|4.811e-08>|<cell|VTP-BASE>|<cell|1.349e-10>|<cell|5.098e-11>|<cell|2.342e-10>>>>>>|Solving
    NETLIB LPs in <with|color|red|<with|font-series|bold|1000>> iterations>
  </shown>|<\hidden>
    <tit|More extensions>

    \;

    <strong|One observation:>

    Pre-conditioning by <math|<A><rsup|\<top\>><around*|(|<A><A><rsup|\<top\>>|)><rsup|-1><A>\<in\>\<bbb-R\><rsup|n\<times\>n>>
    sometimes unfriendly to HSD since <math|\<tau\>\<rightarrow\>0>.

    \;

    Following aspects may accelerate the algorithm

    <\itemize>
      <item>Adjust <math|\<rho\>> adaptively

      <item>More careful matrix scaling to enhance conditioning

      <item>Faster eigen-routine (may be randomized)

      Solving the eigen-problem below

      <\equation*>
        \<lambda\><rsub|min><around*|{|<n><rsub|<x><x>>\<varphi\><around*|(|<x>|)>=\<rho\><around*|(|-<frac|<n>f<around*|(|<x>|)><n>f<around*|(|<x>|)>|f<around*|(|<x>|)><rsup|2>><rsup|\<top\>>+<frac|<A><rsup|\<top\>><A>|f<around*|(|<x>|)>>|)>+<X><rsup|-2>|}>
      </equation*>

      with randomized technique for negative curvature

      <item>More directions
    </itemize>
  </hidden>|<\hidden>
    <tit|>

    \;

    \;

    \;

    \;

    <doc-data|<doc-title|Efficient Eigen-decomposition for Sparse
    Matrices>|<doc-author|<author-data|<\author-affiliation>
      \;

      Discussion 4. Part 2

      \;

      \;

      <date|>
    </author-affiliation>>>>
  </hidden>|<\hidden>
    <tit|Efficient Eigen-decomposition>

    \;

    The eigen routine in <verbatim|HDSDP> is now re-written and presented as
    a stand-alone library <verbatim|SPEIGS>

    It targets the problem

    <\equation*>
      <A><math-bf|V>=\<Lambda\><math-bf|V>
    </equation*>

    for <strong|extremely> sparse <math|<A>\<in\>\<bbb-S\><rsup|n\<times\>n>>
    and needs the full spectrum.

    Recall that in dual-scaling

    <\equation*>
      <around*|\<langle\>|<C>,<bs><rsup|-1><A><bs><rsup|-1>|\<rangle\>>=<big|sum><rsub|i=1><rsup|r<around*|(|<A>|)>>\<lambda\><rsub|i><around*|\<langle\>|<a><rsub|i>,<C><a><rsub|i>|\<rangle\>>.
    </equation*>

    \;

    <\itemize>
      <item>Efficient eigen-routine is critical as a pre-processing step

      <item>Generally hard to exploit sparsity in full-eigen decomposition
      problem

      <item>Benson implements a method that targets decomposition of SDPs
    </itemize>
  </hidden>|<\hidden>
    <tit|Efficient Eigen-decomposition>

    Many SDP coefficient matrices are so sparse that very few entries exist

    Benson's idea: permute the useful entries to a smaller dense block

    \;

    <\equation*>
      <small|<matrix|<tformat|<table|<row|<cell|>|<cell|\<ast\>>|<cell|\<ast\>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|\<ast\>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|\<ast\>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|<with|color|red|\<ast\>>>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|<with|color|red|\<ast\>>>|<cell|>>>>>><large|\<Longrightarrow\>><small|<matrix|<tformat|<cwith|1|1|1|4|cell-tborder|1ln>|<cwith|4|4|1|4|cell-bborder|1ln>|<cwith|5|5|1|4|cell-tborder|1ln>|<cwith|1|4|1|1|cell-lborder|1ln>|<cwith|1|4|4|4|cell-rborder|1ln>|<cwith|1|4|5|5|cell-lborder|1ln>|<table|<row|<cell|>|<cell|\<ast\>>|<cell|\<ast\>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|\<ast\>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|\<ast\>>|<cell|>|<cell|>|<cell|<with|color|red|\<ast\>>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|<with|color|red|\<ast\>>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>>>><large|>>
    </equation*>

    <\equation*>
      <tabular|<tformat|<cwith|1|1|1|4|cell-tborder|1ln>|<cwith|4|4|1|4|cell-bborder|1ln>|<cwith|1|4|1|1|cell-lborder|1ln>|<cwith|1|4|4|4|cell-rborder|1ln>|<table|<row|<cell|>|<cell|\<ast\>>|<cell|\<ast\>>|<cell|>>|<row|<cell|\<ast\>>|<cell|>|<cell|>|<cell|>>|<row|<cell|\<ast\>>|<cell|>|<cell|>|<cell|\<ast\>>>|<row|<cell|>|<cell|>|<cell|\<ast\>>|<cell|>>>>>=<math-bf|V>\<Lambda\><math-bf|V><rsup|\<top\>>
    </equation*>

    Then Lapack is invoked on the smaller system.\ 

    1000x faster than running direct factorization. Simple but useful.
  </hidden>|<\hidden>
    <tit|<verbatim|SPEIGS>>

    \;

    <verbatim|SPEIGS> now implments a standard-alone library for factorizing
    extremely sparse matrices.

    \;

    It implements

    \;

    <\itemize>
      <item>Submatrix detection and permutation

      <item>Diagonal detection

      <item>Givens' rotation

      <item>Rank-one fast detection
    </itemize>

    \;

    for both sparse and dense matrices and works efficiently for SDP matrices\ 

    (and probably real-life matrices that are sparse)

    Available in <verbatim|Matlab> and C.
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
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|1|1>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Reduction in
      potential>|<pageref|auto-1>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|2>||Adding negative
      curvature>|<pageref|auto-2>>
    </associate>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Solving NETLIB LPs
      in <with|color|<quote|red>|<with|font-series|<quote|bold>|1000>>
      iterations>|<pageref|auto-3>>
    </associate>
  </collection>
</auxiliary>