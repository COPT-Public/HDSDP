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

    \;

    \;

    <doc-data|<doc-title|Dimension-reduced Interior Point
    Method>|<doc-author|<author-data|<\author-affiliation>
      \;

      Discussion 8

      \;

      \;

      <date|>
    </author-affiliation>>>>
  </hidden>|<\shown>
    <tit|Current Progress>

    \;

    <strong|Potential reduction>

    \;

    <\itemize>
      <item>Working on C transformation

      <item>The warm-start Lanczos improves by 20% in speed
    </itemize>

    \;

    The solver is designed for solving general problem

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<x>>>|<cell|f<around*|(|<x>|)>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<A><x>=<b>>|<cell|>>|<row|<cell|>|<cell|<x>\<geq\><0>>|<cell|>>>>
    </eqnarray*>

    with smooth convex <math|f<around*|(|<x>|)>> via potential reduction

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<phi\><around*|(|<x>|)>>|<cell|\<assign\>>|<cell|\<rho\>
      log<around*|(|f<around*|(|<x>|)>-z|)>+<big|sum><rsub|i=1><rsup|n>log<around*|(|x<rsub|i>|)>>>>>
    </eqnarray*>

    <\itemize>
      <item>A general framework exploiting curvature in potential reduction

      <item>HSD embedding stands for <math|<A>=<e><rsup|\<top\>>,<b>=1> and
      <math|f<around*|(|<x>|)>=<frac|1|2><around*|\<\|\|\>|<wide|<A>|^><x>|\<\|\|\>><rsup|2>>
    </itemize>

    \;
  </shown>|<\hidden>
    <tit|One useful observation>

    \;

    Recall that\ 

    \;

    <\equation*>
      \<lambda\><rsub|min><around*|(|<X>\<nabla\><rsup|2>\<phi\><around*|(|<x>|)><X>|)>\<leq\><frac|-2\<rho\>|<around*|\<\|\|\>|<X><rsup|-1><around*|(|<x><rsup|\<ast\>>-<x>|)>|\<\|\|\>><rsup|2>>+1.
    </equation*>

    \;

    \;

    <\itemize>
      <item>If <math|x<rsub|i>\<rightarrow\>0> while
      <math|x<rsup|\<ast\>><rsub|i>-x<rsub|i>\<neq\>0>, then the curvature is
      harder to find

      In other words, <with|color|red|centrality> matters when exploiting the
      curvature

      <math|x<rsub|i>+\<alpha\>d<rsub|i>\<rightarrow\>0<rsup|+>> makes next
      curvature hard to detect

      \;

      <item>In practice, we now let line-search go less aggressively to
      ensure centrality

      <math|\<rho\>> is also adjusted to balance centrality and optimality
    </itemize>

    \;

    <strong|HDSDP>

    Now integrating HDSDP into the next COPT release
  </hidden>>
</body>

<\initial>
  <\collection>
    <associate|font-base-size|8>
    <associate|math-color|black>
    <associate|page-medium|paper>
  </collection>
</initial>