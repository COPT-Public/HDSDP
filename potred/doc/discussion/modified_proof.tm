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
  </hide-preamble>

  <doc-data|<doc-title|Modification of proof>|<doc-author|<author-data|<\author-affiliation>
    <date|>
  </author-affiliation>>>>

  In this note we fix a little typo involving the proof for the first-order
  potential reduction.

  Recall that in the note we claim the following relation

  <\equation*>
    <H>\<succeq\>-<frac|\<rho\>|f<around*|(|<x>|)><rsup|2>><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>>=-<frac|2\<rho\>|<with|color|red|<around*|\<\|\|\>|<math-bf|Ax>|\<\|\|\>><rsup|2>>><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A>\<succeq\>-2\<rho\>\<gamma\>
  </equation*>

  but since <math|f<around*|(|<x>|)>=<frac|1|2><around*|\<\|\|\>|<A><x>|\<\|\|\>><rsup|2>>,
  we actually have\ 

  <\equation*>
    -<frac|\<rho\>|f<around*|(|<x>|)><rsup|2>><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>>=-<frac|4\<rho\>|<around*|\<\|\|\>|<math-bf|Ax>|\<\|\|\>><rsup|4>><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A>
  </equation*>

  and the eigen-value bound becomes <math|<x>>-dependent since\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|max<rsub|<u>\<neq\>0><around*|\<langle\>|<frac|<u>|<around*|\<\|\|\>|<u>|\<\|\|\>><rsup|2>>,<frac|4\<rho\>|<around*|\<\|\|\>|<math-bf|Ax>|\<\|\|\>><rsup|4>><A><rsup|\<top\>><A><x><x><rsup|\<top\>><A><rsup|\<top\>><A><u>|\<rangle\>>>>|<row|<cell|>|<cell|=>|<cell|max<rsub|<u>\<neq\>0><frac|4\<rho\><around*|(|<x><A><rsup|\<top\>><A><u>|)><rsup|2>|<around*|\<\|\|\>|<math-bf|Ax>|\<\|\|\>><rsup|4><around*|\<\|\|\>|<u>|\<\|\|\>><rsup|2>>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|2\<rho\>|<frac|1|2><around*|\<\|\|\>|<math-bf|Ax>|\<\|\|\>><rsup|2>>max<rsub|<u>\<neq\>0><frac|<around*|\<\|\|\>|<A><u>|\<\|\|\>><rsup|2>|<around*|\<\|\|\>|<u>|\<\|\|\>><rsup|2>>\<leq\><frac|2\<rho\>\<gamma\>|f<around*|(|<x>|)>>.>>>>
  </eqnarray*>

  By the modified bound and appealing to <math|<around*|(|3|)>> in the
  original proof, we have

  <\equation*>
    <frac|1|2><around*|\<langle\>|<d>,<H><d>|\<rangle\>>\<geq\>-<frac|\<rho\>\<gamma\>|f<around*|(|<x>|)>><around*|\<\|\|\>|<d>|\<\|\|\>><rsup|2>\<geq\>-<frac|\<rho\>\<gamma\>\<beta\><rsup|2>|f<around*|(|<x>|)>>.
  </equation*>

  Though the bound is worse, it does not matter since it is a high-order term
  of <math|\<beta\>>.\ 

  Therefore,

  <\equation*>
    \<phi\><around*|(|<x>+<d>|)>-\<phi\><around*|(|<x>|)>\<leq\>-\<beta\>+<frac|\<rho\>\<gamma\>\<beta\><rsup|2>|2f<around*|(|<x>|)>>+<frac|\<beta\><rsup|2>|2>+<frac|\<rho\>\<gamma\>\<beta\><rsup|2>|f<around*|(|<x>|)>>\<leq\>-\<beta\>+<frac|3\<rho\>\<gamma\>\<beta\><rsup|2>|2f<around*|(|<x>|)>>+<frac|\<beta\><rsup|2>|2>
  </equation*>

  and the rest of the results still apply to guarantee the
  <math|\<cal-O\><around*|(|<frac|1|\<varepsilon\>>log<around*|(|<frac|1|\<varepsilon\>>|)>|)>>
  result.

  <new-page>

  <subsection|Others>

  Then we consider the case of two directions. Recall that

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<phi\><around*|(|<x>|)>>|<cell|\<assign\>>|<cell|\<rho\>
    log<around*|(|f<around*|(|<x>|)>|)>-<big|sum><rsub|i=1><rsup|n>log
    x<rsub|i>>>|<row|<cell|<n>\<phi\><around*|(|<x>|)>>|<cell|\<assign\>>|<cell|<frac|\<rho\><n>f<around*|(|<x>|)>|f<around*|(|<x>|)>>-<X><rsup|-1><e>>>|<row|<cell|<n><rsup|2>\<phi\><around*|(|<x>|)>>|<cell|\<assign\>>|<cell|-<frac|\<rho\><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>>|f<around*|(|<x>|)><rsup|2>>+<frac|\<rho\><n><rsup|2>f<around*|(|<x>|)>|f<around*|(|<x>|)>>+<X><rsup|-2>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<d>>>|<cell|<frac|1|2><around*|\<langle\>|<d>,<n><rsup|2>\<phi\><around*|(|<x>|)><d>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<e><rsup|\<top\>><d>=0>|<cell|>>|<row|<cell|>|<cell|<around*|\<\|\|\>|<X><rsup|-1><d>|\<\|\|\>><rsup|2>\<leq\>\<beta\>.>|<cell|>>>>
  </eqnarray*>

  Taking <math|<d><rprime|'>=<X><rsup|-1><d>>, we have\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<d>>>|<cell|<frac|1|2><around*|\<langle\>|<d><rprime|'>,<X><n><rsup|2>\<phi\><around*|(|<x>|)><X><d><rprime|'>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<x><rsup|\<top\>><d><rprime|'>=0>|<cell|\<nu\>>>|<row|<cell|>|<cell|<around*|\<\|\|\>|<d><rprime|'>|\<\|\|\>><rsup|2>\<leq\>\<beta\>,>|<cell|\<lambda\>>>>>
  </eqnarray*>

  whose solution is given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|<H><d><rprime|'>+\<nu\><x>+2\<lambda\><d><rprime|'>>|<cell|=>|<cell|<0>>>|<row|<cell|<x><rsup|\<top\>><d><rprime|'>>|<cell|=>|<cell|0>>|<row|<cell|<around*|\<\|\|\>|<d><rprime|'>|\<\|\|\>><rsup|2>>|<cell|=>|<cell|\<beta\>,>>>>
  </eqnarray*>

  where

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nu\>>|<cell|=>|<cell|-<frac|<x><rsup|\<top\>><H><d><rprime|'>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>>>|<row|<cell|\<lambda\>>|<cell|=>|<cell|-<frac|<d><rprime|'><rsup|\<top\>><H><d><rprime|'>|2\<beta\>>>>|<row|<cell|\<lambda\>>|<cell|=>|<cell|<frac|<around*|\<\|\|\>|<around*|(|<I>-<frac|<x><rsup|\<top\>><x>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>|)><H><d><rprime|'>|\<\|\|\>>|<around*|\<\|\|\>|<d><rprime|'>|\<\|\|\>><rsup|2>>>>>>
  </eqnarray*>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
  </collection>
</references>