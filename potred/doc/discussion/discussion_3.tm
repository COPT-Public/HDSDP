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
      \;

      Discussion 3

      \;

      \;

      <date|>
    </author-affiliation>>>>
  </shown>|<\hidden>
    <tit|Another direction for LP>

    The original direction composes of <strong|projected gradient> and
    <strong|momentum>

    <\eqnarray*>
      <tformat|<table|<row|<cell|<d><rsup|k>>|<cell|\<leftarrow\>>|<cell|\<alpha\><rsup|g><P><rsub|\<Delta\>><around*|[|\<nabla\>\<varphi\><around*|(|<x><rsup|k>|)>|]>+\<alpha\><rsup|m><around*|(|<x><rsup|k>-<x><rsup|k-1>|)>>>|<row|<cell|<x><rsup|k>>|<cell|\<leftarrow\>>|<cell|<x><rsup|k>+<d><rsup|k>>>>>
    </eqnarray*>

    We can also consider <strong|scaled projected gradient>

    <\equation*>
      <small|<p><around*|(|<x><rsup|k>|)>\<assign\><frac|<X><rsup|k><around*|(|<I>-<frac|<X><rsup|k><e><e><rsup|\<top\>><X><rsup|k>|<around*|\<\|\|\>|<x><rsup|k>|\<\|\|\>><rsup|2>>|)><X><rsup|k>\<nabla\>\<varphi\><around*|(|<x><rsup|k>|)>|<around*|\<\|\|\>|<around*|(|<I>-<frac|<X><rsup|k><e><e><rsup|\<top\>><X><rsup|k>|<around*|\<\|\|\>|<x><rsup|k>|\<\|\|\>><rsup|2>>|)><X><rsup|k>\<nabla\>\<varphi\><around*|(|<x><rsup|k>|)>|\<\|\|\>>>>
    </equation*>

    and build up <math|<d><rsup|k>>.

    <\itemize>
      <item>Empirically better than projected gradient\ 

      <center|<block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<table|<row|<cell|<math|*<around*|(|m,n|)>>/Iteration>|<cell|Projected
      gradient>|<cell|Scaled Projected Gradient>>|<row|<cell|<math|<around*|(|50,100|)>>>|<cell|37>|<cell|33>>|<row|<cell|<math|<around*|(|200,1000|)>>>|<cell|45>|<cell|39>>|<row|<cell|<math|<around*|(|500,2000|)>>>|<cell|46>|<cell|34>>>>>>

      <item>Need more careful tuning for higher accuracy
    </itemize>
  </hidden>|<\hidden>
    <tit|Solving LPs using First-order Potential Reduction>

    A preliminary attempt by simplifying the dual problem

    <\eqnarray*>
      <tformat|<cwith|5|7|3|3|cell-halign|c>|<table|<row|<cell|min<rsub|<x>\<in\>\<bbb-R\><rsup|n>>>|<cell|<c><rsup|\<top\>><x>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<A><x>=<b>>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|max<rsub|<y>>>|<cell|<b><rsup|\<top\>><y>>|<cell|<b><rsup|\<top\>><y><rsup|+>-<b><rsup|\<top\>><y><rsup|->>>|<row|<cell|<text|subject
      to>>|<cell|<A><rsup|\<top\>><y>+<s>=<c>>|<cell|<A><rsup|\<top\>><y><rsup|+>-<rsup|><A><rsup|\<top\>><y><rsup|->+<s>=<c>>>|<row|<cell|>|<cell|<s>\<geq\><0>>|<cell|<y><rsup|+>,<y><rsup|->,<s>\<geq\><0>>>>>
    </eqnarray*>

    and

    <\equation*>
      <tabular|<tformat|<table|<row|<cell|<matrix|<tformat|<table|<row|<cell|<0><rsub|m\<times\>m>>|<cell|<0><rsub|m\<times\>m>>|<cell|<A>>|<cell|<0><rsub|m\<times\>n>>|<cell|<0><rsub|m\<times\>1>>|<cell|-<b>>>|<row|<cell|-<A><rsup|\<top\>>>|<cell|<A><rsup|\<top\>>>|<cell|<0><rsub|n\<times\>n>>|<cell|-<I><rsub|n\<times\>n>>|<cell|<0><rsub|n\<times\>1>>|<cell|<c>>>|<row|<cell|<b><rsup|\<top\>>>|<cell|-<b><rsup|\<top\>>>|<cell|-<c><rsup|\<top\>>>|<cell|<0><rsub|1\<times\>n>>|<cell|-1>|<cell|0>>>>><matrix|<tformat|<table|<row|<cell|<y><rsup|+>>>|<row|<cell|<y><rsup|->>>|<row|<cell|<x>>>|<row|<cell|<s>>>|<row|<cell|\<kappa\>>>|<row|<cell|\<tau\>>>>>>>|<cell|=>|<cell|<0><rsub|m+n+1>,>>>>>
    </equation*>
  </hidden>|<\hidden>
    <tit|Solving LPs>

    \;

    Simplex constrained QP formulation

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<u>=<around*|(|<y>,<x>,<s>,\<kappa\>,\<tau\>|)>>>|<cell|f<around*|(|<u>|)>\<assign\><frac|1|2><around*|\<\|\|\>|<wide|<A>|^><u>|\<\|\|\>><rsup|2>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<e><rsup|\<top\>><u>=1>|<cell|>>|<row|<cell|>|<cell|<u>\<geq\><0>.>|<cell|>>>>
    </eqnarray*>

    Using the potential function

    \;

    <\equation*>
      \<varphi\><around*|(|<u>|)>\<assign\>\<rho\> log
      <around*|(|f<around*|(|<u>|)>|)>-B<around*|(|<y><rsup|+>|)>-B<around*|(|<y><rsup|->|)>-B<around*|(|<x>|)>-B<around*|(|<s>|)>-log
      \<kappa\>-log \<tau\>
    </equation*>

    \;

    and apply the dimension-reduced method.
  </hidden>|<\hidden>
    <tit|Preliminary results>

    \;

    \;

    Solving synthetic LPs, 20000 iterations

    <\big-table>
      <block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<table|<row|<cell|<math|m>>|<cell|<math|n>>|<cell|<math|<around*|\<\|\|\>|<A><x>-<b>|\<\|\|\>>>>|<cell|<math|<around*|\<\|\|\>|<A><rsup|\<top\>><y>+<s>-<c>|\<\|\|\>>>>|<cell|<math|<c><rsup|\<top\>><x>-<b><rsup|\<top\>><y>>>>|<row|<cell|10>|<cell|100>|<cell|<verbatim|3e-06>>|<cell|<verbatim|4e-09>>|<cell|<verbatim|2e-05>>>|<row|<cell|50>|<cell|200>|<cell|<verbatim|6e-04>>|<cell|<verbatim|9e-07>>|<cell|<verbatim|2e-03>>>|<row|<cell|100>|<cell|500>|<cell|<verbatim|3e-05>>|<cell|<verbatim|3e-07>>|<cell|<verbatim|4e-04>>>|<row|<cell|500>|<cell|1000>|<cell|<verbatim|2e-03>>|<cell|<verbatim|5e-06>>|<cell|<verbatim|1e-03>>>>>>
    </big-table|Synthetic tests>

    \;

    <\itemize>
      <item>Possible to solve LPs to low accuracy

      <item>Needs much tuning for real-life LPs (e.g., <verbatim|Netlib>)

      <item>Or just using HSD without <math|<y>=<y><rsup|+>-<y><rsup|->> /
      primal-dual potential function
    </itemize>
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
    <associate|auto-1|<tuple|1|1>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Synthetic
      tests>|<pageref|auto-1>>
    </associate>
  </collection>
</auxiliary>