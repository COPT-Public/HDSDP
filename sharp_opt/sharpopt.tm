<TeXmacs|2.1>

<style|<tuple|generic|centered-program|number-europe>>

<\body>
  <doc-data|<doc-title|Stochastic Model-based Algorithm can be Accelerated by
  Minibatching for Sharp Functions>|<\doc-author>
    <date|>
  </doc-author>>

  <section|Literature Review>

  <big-table|<block|<tformat|<cwith|1|-1|1|-1|cell-valign|c>|<cwith|6|6|2|2|cell-row-span|4>|<cwith|6|6|2|2|cell-col-span|1>|<cwith|2|2|2|2|cell-row-span|4>|<cwith|2|2|2|2|cell-col-span|1>|<cwith|2|2|3|3|cell-row-span|2>|<cwith|2|2|3|3|cell-col-span|1>|<cwith|4|4|3|3|cell-row-span|2>|<cwith|4|4|3|3|cell-col-span|1>|<cwith|6|6|3|3|cell-row-span|2>|<cwith|6|6|3|3|cell-col-span|1>|<cwith|8|8|3|3|cell-row-span|2>|<cwith|8|8|3|3|cell-col-span|1>|<cwith|14|14|2|2|cell-row-span|4>|<cwith|14|14|2|2|cell-col-span|1>|<cwith|10|10|2|2|cell-row-span|4>|<cwith|10|10|2|2|cell-col-span|1>|<cwith|10|10|3|3|cell-row-span|2>|<cwith|10|10|3|3|cell-col-span|1>|<cwith|12|12|3|3|cell-row-span|2>|<cwith|12|12|3|3|cell-col-span|1>|<cwith|14|14|3|3|cell-row-span|2>|<cwith|14|14|3|3|cell-col-span|1>|<cwith|16|16|3|3|cell-row-span|2>|<cwith|16|16|3|3|cell-col-span|1>|<cwith|2|2|1|1|cell-row-span|8>|<cwith|2|2|1|1|cell-col-span|1>|<cwith|10|10|1|1|cell-row-span|8>|<cwith|10|10|1|1|cell-col-span|1>|<cwith|2|17|1|5|cell-valign|c>|<cwith|1|-1|1|-1|cell-halign|c>|<cwith|6|9|3|3|cell-valign|c>|<cwith|6|6|3|3|cell-row-span|2>|<cwith|6|6|3|3|cell-col-span|1>|<cwith|8|8|3|3|cell-row-span|2>|<cwith|8|8|3|3|cell-col-span|1>|<cwith|6|9|3|3|cell-valign|c>|<cwith|6|9|3|3|cell-halign|c>|<cwith|10|13|3|3|cell-valign|c>|<cwith|10|10|3|3|cell-row-span|2>|<cwith|10|10|3|3|cell-col-span|1>|<cwith|12|12|3|3|cell-row-span|2>|<cwith|12|12|3|3|cell-col-span|1>|<cwith|10|13|3|3|cell-valign|c>|<cwith|10|13|3|3|cell-halign|c>|<cwith|14|17|3|3|cell-valign|c>|<cwith|14|14|3|3|cell-row-span|2>|<cwith|14|14|3|3|cell-col-span|1>|<cwith|16|16|3|3|cell-row-span|2>|<cwith|16|16|3|3|cell-col-span|1>|<cwith|14|17|3|3|cell-valign|c>|<cwith|14|17|3|3|cell-halign|c>|<cwith|10|17|2|2|cell-valign|c>|<cwith|14|14|2|2|cell-row-span|4>|<cwith|14|14|2|2|cell-col-span|1>|<cwith|10|10|2|2|cell-row-span|4>|<cwith|10|10|2|2|cell-col-span|1>|<cwith|10|17|2|2|cell-valign|c>|<cwith|10|17|2|2|cell-halign|c>|<cwith|6|7|4|4|cell-valign|c>|<cwith|6|7|4|4|cell-valign|c>|<cwith|6|7|4|4|cell-halign|c>|<cwith|8|9|4|4|cell-valign|c>|<cwith|8|9|4|4|cell-valign|c>|<cwith|8|9|4|4|cell-halign|c>|<cwith|10|11|4|4|cell-valign|c>|<cwith|10|11|4|4|cell-valign|c>|<cwith|10|11|4|4|cell-halign|c>|<cwith|12|13|4|4|cell-valign|c>|<cwith|12|13|4|4|cell-valign|c>|<cwith|12|13|4|4|cell-halign|c>|<cwith|14|15|4|4|cell-valign|c>|<cwith|14|15|4|4|cell-valign|c>|<cwith|14|15|4|4|cell-halign|c>|<cwith|16|17|4|4|cell-valign|c>|<cwith|16|17|4|4|cell-valign|c>|<cwith|16|17|4|4|cell-halign|c>|<table|<row|<cell|Algorithm>|<cell|Convexity>|<cell|Randomness>|<cell|Stepsize>|<cell|Complexity>>|<row|<cell|SGD>|<cell|Convex>|<cell|Determinic>|<cell|Constant>|<cell|<math|
  log<around*|(|1/\<varepsilon\>|)>>>>|<row|<cell|>|<cell|>|<cell|>|<cell|Geometrically>|<cell|<math|
  log<around*|(|1/\<varepsilon\>|)>>>>|<row|<cell|>|<cell|>|<cell|Stochastic>|<cell|Constant>|<cell|\U>>|<row|<cell|>|<cell|>|<cell|>|<cell|Geometrically
  >|<cell|<math| log<around*|(|1/\<varepsilon\>|)>>>>|<row|<cell|>|<cell|Weakly>|<cell|Deterministic>|<cell|Constant>|<cell|<math|
  log<around*|(|1/\<varepsilon\>|)>>>>|<row|<cell|>|<cell|>|<cell|>|<cell|Geometrically>|<cell|<math|
  log<around*|(|1/\<varepsilon\>|)>>>>|<row|<cell|>|<cell|>|<cell|Stochastic>|<cell|Constant>|<cell|\U>>|<row|<cell|>|<cell|>|<cell|>|<cell|Geometrically>|<cell|<math|
  log<around*|(|1/\<varepsilon\>|)>>>>|<row|<cell|SPL/SPP>|<cell|Convex>|<cell|Deterministic>|<cell|Constant>|<cell|<math|log
  log<around*|(|1/\<varepsilon\>|)>>>>|<row|<cell|>|<cell|>|<cell|>|<cell|Geometrically>|<cell|<with|color|red|Needed>>>|<row|<cell|>|<cell|>|<cell|Stochastic>|<cell|Constant>|<cell|<with|color|red|<math|log<around*|(|1/\<varepsilon\>|)><rsup|\<dag\>>>>>>|<row|<cell|>|<cell|>|<cell|>|<cell|Geometrically>|<cell|<math|
  log<around*|(|1/\<varepsilon\>|)>>>>|<row|<cell|>|<cell|Weakly>|<cell|Deterministic>|<cell|Constant>|<cell|<math|log
  log<around*|(|1/\<varepsilon\>|)>>>>|<row|<cell|>|<cell|>|<cell|>|<cell|Geometrically>|<cell|<with|color|red|Needed>>>|<row|<cell|>|<cell|>|<cell|Stochastic>|<cell|Constant>|<cell|<with|color|red|Needed>>>|<row|<cell|>|<cell|>|<cell|>|<cell|Geometrically>|<cell|<math|
  log<around*|(|1/\<varepsilon\>|)>>>>>>>|Literature over optimization with
  sharpness>

  <math|\<dag\>>: minibatch acceleration is already proven for easy problems
  (<math|arg min<rsub|x> f<around*|(|x,\<xi\>|)>=x<rsup|\<ast\>>,\<forall\>\<xi\>>).

  <section|Preliminaries>

  Consider the following optimization problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|x\<in\>\<cal-X\>>>|<cell|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|x,\<xi\>|)>|]>>|<cell|>>>>
  </eqnarray*>

  <strong|Assumption 1>. It is possible to sample i.i.d.
  <math|<around*|{|\<xi\><rsub|1>,\<ldots\>,\<xi\><rsub|n>|}>>.

  <strong|Assumption 2>. <math|f> is <math|\<lambda\>>-weakly convex.

  We assume that <math|f+<frac|\<lambda\>|2><around*|\<\|\|\>|x|\<\|\|\>><rsup|2>>
  is convex.

  <strong|Assumption 3>. <math|f> is sharp. In other words,

  <\equation*>
    \<mu\>\<cdummy\>dist<around*|(|x,\<cal-X\><rsup|\<ast\>>|)>\<leq\>f<around*|(|x|)>-f<rsup|\<ast\>>,\<forall\>x\<in\>\<cal-X\><rsup|\<ast\>>,
  </equation*>

  where <math|\<cal-X\><rsup|\<ast\>>> is the set of optimal solutions to the
  problem.

  <strong|Assumption 4>. <math|f> is locally Lipschitz-continuous.\ 

  Define the tube <math|\<cal-T\><rsub|\<gamma\>>\<assign\><around*|{|x\<in\>\<cal-X\>:dist<around*|(|x,\<cal-X\><rsup|\<ast\>>|)>\<leq\><frac|\<gamma\>\<mu\>|\<tau\>>|}>>
  and we have

  <\equation*>
    min<rsub|g\<in\>\<partial\>f<rsub|x><around*|(|x,\<xi\>|)>><around*|\<\|\|\>|g|\<\|\|\>>\<leq\>L,\<forall\>x\<in\>\<cal-T\><rsub|2>,\<xi\>.
  </equation*>

  <strong|Assumption 5>. Two-sided accuracy is available. i.e.,

  <\equation*>
    <around*|\||f<around*|(|y|)>-f<rsub|x><around*|(|y,\<xi\>|)>|\|>\<leq\><frac|\<tau\>|2><around*|\<\|\|\>|x-y|\<\|\|\>><rsup|*2>.
  </equation*>

  It is already known that in the convex case, the proximal point method
  converges quadratically [1] and its stochastic variant has linear
  convergence when using a geometrically decaying stepsize [2]. Hence there
  is space for acceleration.

  <section|Convex Optimization>

  To analyze the case of convex optimization, we specially let
  <math|\<lambda\>=0> and further assume that global Lipschitzness of the
  model <math|f<rsub|x><around*|(|\<cdummy\>,\<xi\>|)>> holds.\ 

  <subsection|Restarting Strategy with Decaying Stepsize>

  <\lemma>
    The algorithm in [SMOD] initialized with <math|y<rsub|0>> and satisfies\ 

    <\equation*>
      \<bbb-E\><around*|[|f<around*|(|x<rsup|K+1>|)>-f<rsup|\<ast\>>|]>\<leq\><frac|2\<tau\>dist<rsup|2><around*|(|y<rsub|0>,\<cal-X\><rsup|\<ast\>>|)>|*<around*|(|K+1|)><around*|(|K+2|)>>+<frac|4<sqrt|2>L
      dist<around*|(|y<rsub|0>,\<cal-X\><rsup|\<ast\>>|)>|<sqrt|3m<rsub|t><around*|(|K+1|)>>>.
    </equation*>
  </lemma>

  <\lemma>
    For some growth function <math|g\<gtr\>0>, denote
    <math|E<rsub|t>\<assign\><around*|{|dist<around*|(|x<rsub|t>,\<cal-X\><rsup|\<ast\>>|)>\<leq\><frac|R<rsub|0>|g<around*|(|t|)>>|}>>
    and we have the following relation holds

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<bbb-P\><around*|(|E<rsub|T>|)>>|<cell|\<geq\>>|<cell|1-<big|sum><rsub|t=0><rsup|T-1><around*|[|<frac|2\<tau\>R<rsub|0>|\<mu\>K<rsup|2>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)><rsup|2>>+<frac|4<sqrt|6>L
      |3<sqrt|m<rsub|t> <around*|(|K+1|)>>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)>>|]>.>>>>
    </eqnarray*>
  </lemma>

  <\proof>
    Without loss of generality we have

    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|>|<cell|\<bbb-P\><around*|(|E<rsub|t+1>|)>>>|<row|<cell|>|<cell|=>|<cell|\<bbb-P\><around*|(|E<rsub|t+1>\|<wide|E<rsub|t>|\<bar\>>|)>\<bbb-P\><around*|(|<wide|E<rsub|t>|\<bar\>>|)>+\<bbb-P\><around*|(|E<rsub|t>|)>\<bbb-P\><around*|(|E<rsub|t+1>\|E<rsub|t>|)>\<bbb-P\><around*|(|E<rsub|t>|)>>>|<row|<cell|>|<cell|\<geq\>>|<cell|\<bbb-P\><around*|(|E<rsub|t>|)>\<bbb-P\><around*|(|E<rsub|t+1>\|E<rsub|t>|)>>>>>
    </eqnarray*>

    and that

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<bbb-P\><around*|(|E<rsub|t+1>\|E<rsub|t>|)>>|<cell|=>|<cell|1-\<bbb-P\><around*|(|<wide|E<rsub|t+1>|\<bar\>>\|E<rsub|t>|)>>>|<row|<cell|>|<cell|=>|<cell|1-\<bbb-P\><around*|(|dist<around*|(|x<rsub|t+1>,\<cal-X\><rsup|\<ast\>>|)>\<geq\><frac|R<rsub|0>|g<around*|(|t+1|)>>\|E<rsub|t>|)>>>|<row|<cell|>|<cell|\<geq\>>|<cell|1-<frac|\<bbb-E\><around*|[|dist<around*|(|x<rsub|t+1>,\<cal-X\><rsup|\<ast\>>|)><mid|\|>E<rsub|t>|]>|R<rsub|0>/g<around*|(|t+1|)>>>>|<row|<cell|>|<cell|=>|<cell|1-<frac|\<bbb-E\><around*|[|dist<around*|(|x<rsub|t+1>,\<cal-X\><rsup|\<ast\>>|)>\<bbb-I\><around*|{|E<rsub|t>|}>|]>|R<rsub|0>/g<around*|(|t+1|)>><frac|1|\<bbb-P\><around*|(|E<rsub|t>|)>>,>>>>
    </eqnarray*>

    where the inequality is by Markov's inequality.\ 

    Then we consider\ 

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<bbb-E\><around*|[|dist<around*|(|x<rsub|t+1>,\<cal-X\><rsup|\<ast\>>|)>\<bbb-I\><around*|{|E<rsub|t>|}>|]>>|<cell|\<leq\>>|<cell|<frac|1|\<mu\>>\<bbb-E\><around*|[|<around*|(|f<around*|(|x<rsub|t+1>|)>-f<rsup|\<ast\>>|)>\<bbb-I\><around*|{|E<rsub|t>|}>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|1|\<mu\>><around*|{|<frac|2\<tau\>\<bbb-E\><around*|[|dist<rsup|2><around*|(|x<rsub|t>,\<cal-X\><rsup|\<ast\>>|)>\<bbb-I\><around*|{|E<rsub|t>|}>|]>|<around*|(|K+1|)><around*|(|K+2|)>>+<frac|4<sqrt|2>L\<bbb-E\><around*|[|dist<around*|(|x<rsub|t>,\<cal-X\><rsup|\<ast\>>|)>\<bbb-I\><around*|{|E<rsub|t>|}>|]>|<sqrt|3m<rsub|t><around*|(|K+1|)>>>|}>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|2\<tau\>R<rsub|0><rsup|2>|\<mu\>K<rsup|2>>\<cdummy\><frac|1|g<around*|(|t|)><rsup|2>>+<frac|4<sqrt|6>L
      R<rsub|0>|3\<mu\><sqrt|m<rsub|t> <around*|(|K+1|)>>>\<cdummy\><frac|1|g<around*|(|t|)>>.>>>>
    </eqnarray*>

    Next we combine the above and obtain that

    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|>|<cell|\<bbb-P\><around*|(|E<rsub|t+1>|)>>>|<row|<cell|>|<cell|\<geq\>>|<cell|\<bbb-P\><around*|(|E<rsub|t>|)><around*|{|1-<frac|\<bbb-E\><around*|[|dist<around*|(|x<rsub|t+1>,\<cal-X\><rsup|\<ast\>>|)>\<bbb-I\><around*|{|E<rsub|t>|}>|]>|R<rsub|0>/g<around*|(|t+1|)>><frac|1|\<bbb-P\><around*|(|E<rsub|t>|)>>|}>>>|<row|<cell|>|<cell|=>|<cell|\<bbb-P\><around*|(|E<rsub|t>|)>-<frac|\<bbb-E\><around*|[|dist<around*|(|x<rsub|t+1>,\<cal-X\><rsup|\<ast\>>|)>\<bbb-I\><around*|{|E<rsub|t>|}>|]>|R<rsub|0>/g<around*|(|t+1|)>>>>|<row|<cell|>|<cell|\<geq\>>|<cell|\<bbb-P\><around*|(|E<rsub|t>|)>-<around*|[|<frac|2\<tau\>R<rsub|0>|\<mu\>K<rsup|2>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)><rsup|2>>+<frac|4<sqrt|6>L
      |3\<mu\><sqrt|m<rsub|t> <around*|(|K+1|)>>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)>>|]>.>>>>
    </eqnarray*>

    Summing over <math|t=0,\<ldots\>,T-1> gives

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<bbb-P\><around*|(|E<rsub|T>|)>>|<cell|\<geq\>>|<cell|1-<big|sum><rsub|t=0><rsup|T-1><around*|[|<wide*|<frac|2\<tau\>R<rsub|0>|\<mu\>K<rsup|2>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)><rsup|2>>|\<wide-underbrace\>><rsub|<text|Quadratic>>+<wide*|<frac|4<sqrt|6>L
      |3\<mu\><sqrt|m<rsub|t><around*|(|K+1|)>>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)>>|\<wide-underbrace\>><rsub|<text|Linear>>|]>>>>>
    </eqnarray*>

    \;
  </proof>

  <\remark>
    For <verbatim|SPP> algorithm we have <math|\<tau\>=0> and the quadratic
    acceleration term is not present and we hence have

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<bbb-P\><around*|(|E<rsub|T>|)>>|<cell|\<geq\>>|<cell|1-<frac|4<sqrt|6>L
      |3\<mu\><sqrt|m <around*|(|K+1|)>>><big|sum><rsub|t=0><rsup|T-1><frac|g<around*|(|t+1|)>|g<around*|(|t|)>>>>>>
    </eqnarray*>
  </remark>

  <\remark>
    To recover the deterministic quadratic convergence, we let
    <math|m\<rightarrow\>\<infty\>> and get

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<bbb-P\><around*|(|E<rsub|T>|)>>|<cell|\<geq\>>|<cell|1-<frac|2\<tau\>R<rsub|0>|\<mu\>K<rsup|2>><big|sum><rsub|t=0><rsup|T-1><frac|g<around*|(|t+1|)>|g<around*|(|t|)><rsup|2>>>>>>
    </eqnarray*>

    and this allows us to take growth function to
    <math|g<around*|(|t|)>=2<rsup|2<rsup|t>>> such that
    <math|<frac|g<around*|(|t+1|)>|g<around*|(|t|)><rsup|2>>=2=\<cal-O\><around*|(|1|)>>.
    Then we can follow [<em|Dmitri>] to recover the quadratic convergence.
  </remark>

  Now we analyze the way to choose <math|<around*|(|g,<around*|{|m<rsub|t>|}>|)>>
  for faster convergence.

  Consider taking <math|m<rsub|t>=m<around*|(|t|)>> and we get\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-P\><around*|(|E<rsub|T>|)>>|<cell|\<geq\>>|<cell|1-<big|sum><rsub|t=0><rsup|T-1><around*|(|<frac|2\<tau\>R<rsub|0>|\<mu\>K<rsup|2>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)><rsup|2>>+<frac|4<sqrt|6>L
    |3\<mu\><sqrt|K+1>>\<cdummy\><frac|g<around*|(|t+1|)>|<sqrt|m<rsub|t>>g<around*|(|t|)>>|)>.>>>>
  </eqnarray*>

  For brevity we first consider the proximal point method with
  <math|\<tau\>=0> and we get the bound

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-P\><around*|(|E<rsub|T>|)>>|<cell|\<geq\>>|<cell|1-<big|sum><rsub|t=0><rsup|T-1><around*|(|<frac|4<sqrt|6>L
    |3\<mu\><sqrt|K<rsub|t>+1>>\<cdummy\><frac|g<around*|(|t+1|)>|<sqrt|m<rsub|t>>g<around*|(|t|)>>|)>.>>>>
  </eqnarray*>

  <\strong>
    Super-linear Batchsize
  </strong>

  Take <math|g<around*|(|t|)>=2<rsup|t<rsup|2>>> and we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-P\><around*|(|E<rsub|T>|)>>|<cell|\<geq\>>|<cell|1-<big|sum><rsub|t=0><rsup|T-1><around*|(|<frac|8<sqrt|6>L
    |3\<mu\><sqrt|K<rsub|t>+1>>\<cdummy\><frac|4<rsup|t>|<sqrt|m<rsub|t>>>|)>.>>>>
  </eqnarray*>

  Take <math|m<rsub|t>=16<rsup|t>,T=<around*|\<lceil\>|<sqrt|log<rsub|2><around*|(|<frac|R<rsub|0>|\<varepsilon\>>|)>>|\<rceil\>>>
  and <math|K<rsub|t>\<equiv\><around*|\<lfloor\>|<frac|128T<rsup|2>|3>\<cdummy\><around*|(|<frac|L|\<delta\>\<mu\>>|)><rsup|2>|\<rfloor\>>>,
  we have the total sample complexity of\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsub|t=0><rsup|T-1>m<rsub|t>K<rsub|t>>|<cell|=>|<cell|<frac|128T<rsup|2>|3><around*|(|<frac|L|\<delta\>\<mu\>>|)><rsup|2><big|sum><rsub|t=0><rsup|T-1>16<rsup|t>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|128T<rsup|2>|45><around*|(|<frac|L|\<delta\>\<mu\>>|)><rsup|2>exp<around*|(|4<sqrt|
    log<rsub|2><around*|(|<frac|R<rsub|0>|\<varepsilon\>>|)>>|)><rsup|>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|128
    log<rsub|2><around*|(|<frac|R<rsub|0>|\<varepsilon\>>|)>|45><around*|(|<frac|L|\<delta\>\<mu\>>|)><rsup|2>exp<around*|(|4<sqrt|
    log<rsub|2><around*|(|<frac|R<rsub|0>|\<varepsilon\>>|)>>|)>>>>>
  </eqnarray*>

  <strong|Optimal Choice for Parameters>

  Last we consider the general choice of <math|g<around*|(|t|)>,m<rsub|t>>
  and <math|K<rsub|t>>. For brevity we use <math|m<around*|(|t|)>> and
  <math|K<around*|(|t|)>> as functions of discrete values <math|t>. Then due
  to monotonicity of <math|g> we have <math|T=g<rsup|-1><around*|(|t|)>> and
  <math|>that

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-P\><around*|(|E<rsub|T>|)>>|<cell|\<geq\>>|<cell|1-<big|sum><rsub|t=0><rsup|g<rsup|-1><around*|(|R<rsub|0>/\<varepsilon\>|)>-1><around*|(|<frac|8<sqrt|6>L
    |3\<mu\><sqrt|K<around*|(|t|)>+1>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)><sqrt|m<around*|(|t|)>>>|)>.>>>>
  </eqnarray*>

  Also, we have the total sample complexity given by

  <\equation*>
    <big|sum><rsub|t=0><rsup|g<rsup|-1><around*|(|R<rsub|0>/\<varepsilon\>|)>-1>m<around*|(|t|)>K<around*|(|t|)>.
  </equation*>

  Then we use <math|K<around*|(|t|)>+1> to replace <math|K<around*|(|t|)>>
  and get an abstract optimization problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|g,m,K>>|<cell|<big|sum><rsub|t=0><rsup|g<rsup|-1><around*|(|R<rsub|0>/\<varepsilon\>|)>-1>m<around*|(|t|)>K<around*|(|t|)>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<big|sum><rsub|t=0><rsup|g<rsup|-1><around*|(|R<rsub|0>/\<varepsilon\>|)>-1><around*|(|<frac|8<sqrt|6>L
    |3\<mu\>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)><sqrt|m<around*|(|t|)>K<around*|(|t|)>>>|)>\<leq\>\<delta\>>|<cell|.>>>>
  </eqnarray*>

  To solve the problem, we first denote <math|\<alpha\>\<assign\>R<rsub|0>/\<varepsilon\>,\<theta\>\<assign\><frac|<sqrt|6>\<mu\>\<delta\>|16L>,u<around*|(|t|)>\<assign\>m<around*|(|t|)>K<around*|(|t|)>>
  and get\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|g,u>>|<cell|<big|sum><rsub|t=0><rsup|g<rsup|-1><around*|(|\<alpha\>|)>-1>u<around*|(|t|)>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<big|sum><rsub|t=0><rsup|g<rsup|-1><around*|(|\<alpha\>|)>-1><frac|1|<sqrt|u<around*|(|t|)>>>\<cdummy\><frac|g<around*|(|t+1|)>|g<around*|(|t|)>>\<leq\>\<theta\>>|<cell|.>>>>
  </eqnarray*>

  Now we consider the following cases.

  <strong|Linear Convergence>

  In this case we have <math|<frac|g<around*|(|t+1|)>|g<around*|(|t|)>>=\<beta\>>
  and by optimality condition we know that it is optimal to let
  <math|u<around*|(|t<rsub|1>|)>=u<around*|(|t<rsub|2>|)>,\<forall\>t<rsub|1>,t<rsub|2>>
  and the constraint is transformed into

  <\equation*>
    <frac|log<rsub|\<beta\>><around*|(|\<alpha\>|)>|<sqrt|u<around*|(|0|)>>>\<leq\>\<theta\>/\<beta\>\<Rightarrow\>u<around*|(|0|)>\<geq\><frac|\<beta\><rsup|2>
    log<rsup|2><rsub|\<beta\>><around*|(|\<alpha\>|)>|\<theta\><rsup|2>>=<frac|128L<rsup|2>\<beta\><rsup|2>
    log<rsup|2><rsub|\<beta\>><around*|(|\<alpha\>|)>|3\<mu\><rsup|2>\<delta\><rsup|2>>.
  </equation*>

  Also the objective is into

  <\equation*>
    <big|sum><rsub|t=0><rsup|g<rsup|-1><around*|(|\<alpha\>|)>-1>u<around*|(|t|)>=log<rsub|\<beta\>><around*|(|\<alpha\>|)>u<around*|(|0|)>\<geq\><around*|(|<frac|\<beta\>
    |log<rsup|3><around*|(|\<beta\>|)>>|)><around*|(|<frac|128L<rsup|2>|3\<mu\><rsup|2>\<delta\><rsup|2>>|)>log<rsup|3><around*|(|\<alpha\>|)>.
  </equation*>

  Hence the best bound in terms of linear convergence is attained by
  <math|\<beta\>=e<rsup|3>\<Rightarrow\><frac|\<beta\>
  |log<rsup|3><around*|(|\<beta\>|)>>=<frac|e<rsup|3>|27>> with constant
  batchsize and this gives the best sample complexity

  <\equation*>
    <frac|128e<rsup|3>|81><around*|(|<frac|L<rsup|2>|\<mu\><rsup|2>\<delta\><rsup|2>>|)>log<rsup|3><around*|(|<frac|R<rsub|0>|\<varepsilon\>>|)>.
  </equation*>

  <strong|Constant Sample per Iteration>

  In this case we assume that <math|u<around*|(|t|)>\<equiv\>u> and we have\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|g,u>>|<cell|g<rsup|-1><around*|(|\<alpha\>|)>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<big|sum><rsub|t=0><rsup|g<rsup|-1><around*|(|\<alpha\>|)>-1><frac|g<around*|(|t+1|)>|g<around*|(|t|)>>\<leq\>\<theta\><sqrt|u>>|<cell|.>>>>
  </eqnarray*>

  Or more abstractly, we have to solve

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|f>>|<cell|f<rsup|-1><around*|(|\<alpha\>|)>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<big|int><rsub|0><rsup|f<rsup|-1><around*|(|\<alpha\>|)>><frac|f<around*|(|x+1|)>|f<around*|(|x|)>>
    d x\<leq\>1>|<cell|>>>>
  </eqnarray*>

  \;

  <\strong>
    Super-linear <math|exp<around*|(|t log <around*|(|t+1|)>|)>>
  </strong>

  In this case we have <math|<frac|g<around*|(|t+1|)>|g<around*|(|t|)>>=<around*|(|1+<frac|1|t+1>|)><rsup|t><around*|(|t+2|)>>
  and in this case we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|g,u>>|<cell|<big|sum><rsub|t=0><rsup|W<around*|(|R<rsub|0>/\<varepsilon\>|)>-1>u<around*|(|t|)>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<big|sum><rsub|t=0><rsup|W<around*|(|e
    R<rsub|0>/\<varepsilon\>|)>-2><frac|1|<sqrt|u<around*|(|t|)>>>\<cdummy\><around*|(|1+<frac|1|t>|)><rsup|t><around*|(|t+1|)>\<leq\>\<theta\>>|<cell|,>>>>
  </eqnarray*>

  where <math|W<around*|(|x|)>> is the Lambert-W function. By taking
  <math|m<around*|(|t|)>\<equiv\>m,K<around*|(|t|)>=<frac|512L<rsup|2>e<rsup|2>|3m\<mu\><rsup|2>\<delta\><rsup|2>>log<rsup|4><around*|(|<frac|R<rsub|0>|\<varepsilon\>>|)>>
  we have the sample complexity of <math|o<around*|(|<frac|512L<rsup|2>|3\<mu\><rsup|2>\<delta\><rsup|2>>log<rsup|5><around*|(|<frac|R<rsub|0>|\<varepsilon\>>|)>|)>>.
  Hence we achieve super-linear convergence.

  \;

  <strong|Super-linear <math|exp<around*|(|\<cal-P\><around*|(|t|)>|)>>>

  In this case we consider a special case of super-linear convergence with
  <math|g<around*|(|t|)>=e<rsup|\<beta\>t<rsup|p>>>. In this case we have
  <math|<frac|g<around*|(|t+1|)>|g<around*|(|t|)>>=exp<around*|(|\<beta\><around*|(|t+1|)><rsup|p>-\<beta\>t<rsup|p>|)>>
  and <math|T=g<rsup|-1><around*|(|\<alpha\>|)>=<around*|[|log<around*|(|\<alpha\>|)>|]><rsup|1/p>>.
  Hence we have the optimization problem given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|p,u>>|<cell|<big|sum><rsub|t=0><rsup|<around*|[|log<around*|(|\<alpha\>|)>|]><rsup|1/p>-1>u<around*|(|t|)>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<big|sum><rsub|t=0><rsup|<around*|[|log<around*|(|\<alpha\>|)>|]><rsup|1/p>-1><frac|1|<sqrt|u<around*|(|t|)>>>\<cdummy\>exp<around*|(|<around*|(|t+1|)><rsup|p>-t<rsup|p>|)>\<leq\>\<theta\>.>|<cell|>>>>
  </eqnarray*>

  A trivial selection is <math|p=2> and <math|<frac|g<around*|(|t+1|)>|g<around*|(|t|)>>=exp<around*|(|2\<beta\>t+\<beta\>|)>>.
  Then we have <math|><math|<around*|[|log<around*|(|\<alpha\>|)>|]><rsup|1/p>-1=<sqrt|log<around*|(|\<alpha\>|)>>-1>,
  giving

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|u>>|<cell|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1>u<around*|(|t|)>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1><frac|exp<around*|(|2\<beta\>t|)>|<sqrt|u<around*|(|t|)>>>\<leq\>\<theta\>e<rsup|-\<beta\>>.>|<cell|>>>>
  </eqnarray*>

  Then by writing the Lagrangian function

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-L\><around*|(|<around*|{|u<around*|(|t|)>|}>,\<lambda\>|)>>|<cell|\<assign\>>|<cell|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1>u<around*|(|t|)>-\<lambda\><around*|(|\<theta\>e<rsup|-\<beta\>>-<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1><frac|exp<around*|(|2\<beta\>t|)>|<sqrt|u<around*|(|t|)>>>|)>>>>>
  </eqnarray*>

  we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|u<around*|(|t|)>>\<cal-L\>>|<cell|=>|<cell|1-<frac|\<lambda\>exp<around*|(|2\<beta\>t|)>|2>u<around*|(|t|)><rsup|-3/2>>>>>
  </eqnarray*>

  and that\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|u<around*|(|t|)>>|<cell|=>|<cell|\<lambda\><rsup|2/3><around*|(|<frac|exp<around*|(|2\<beta\>t|)>|2>|)><rsup|2/3>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1><frac|exp<around*|(|2\<beta\>t|)>|<sqrt|u<around*|(|t|)>>>>|<cell|=>|<cell|<around*|[|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1><frac|exp<around*|(|2\<beta\>t|)><rsup|2/3>|2<rsup|-1/3>>|]>\<lambda\><rsup|-1/3>>>|<row|<cell|>|<cell|=>|<cell|\<theta\>e<rsup|-\<beta\>>.>>>>
  </eqnarray*>

  Hence we have <math|\<lambda\><rsup|\<ast\>2/3>=<frac|<around*|(|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1><frac|exp<around*|(|2\<beta\>t|)><rsup|2/3>|2<rsup|-1/3>>|)><rsup|2>|\<theta\><rsup|2>e<rsup|-2\<beta\>>>>
  and\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|u<around*|(|t|)>>|<cell|=>|<cell|<frac|<around*|(|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1>exp<around*|(|2\<beta\>t|)><rsup|2/3>|)><rsup|2>|\<theta\><rsup|2>e<rsup|-2\<beta\>>><around*|(|exp<around*|(|2\<beta\>t|)>|)><rsup|2/3>,>>>>
  </eqnarray*>

  giving

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1>u<around*|(|t|)>>|<cell|=>|<cell|<frac|1|\<theta\><rsup|2>e<rsup|-2\<beta\>>><around*|(|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1>exp<around*|(|4\<beta\>t/3|)>|)><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<theta\><rsup|2>e<rsup|-2\<beta\>>><around*|(|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1>exp<around*|(|4\<beta\>/3|)><rsup|t>|)><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<theta\><rsup|2>e<rsup|-2\<beta\>>><around*|(|<big|sum><rsub|t=0><rsup|<sqrt|log<around*|(|\<alpha\>|)>>-1>exp<around*|(|4\<beta\>/3|)><rsup|t>|)><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<theta\><rsup|2>e<rsup|-2\<beta\>>><around*|(|<frac|exp<around*|(|4\<beta\>/3|)><rsup|<sqrt|log<around*|(|\<alpha\>|)>>>-1<rsup|>|exp<around*|(|4\<beta\>/3|)>-1>|)><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<theta\><rsup|2>><frac|exp<around*|(|2\<beta\>|)><around*|(|exp<around*|(|4\<beta\>/3|)><rsup|<sqrt|log<around*|(|\<alpha\>|)>>>-1<rsup|>|)><rsup|2>|<around*|[|exp<around*|(|4\<beta\>/3|)>-1|]><rsup|2>>.>>>>
  </eqnarray*>

  For some given <math|\<beta\>>, we get the total complexity of\ 

  <\equation*>
    <frac|128L<rsup|2>|3\<mu\><rsup|2>\<delta\><rsup|2>>\<cdummy\><frac|exp<around*|(|2\<beta\>|)><around*|(|exp<around*|(|4\<beta\>/3|)><rsup|<sqrt|log<around*|(|\<alpha\>|)>>>-1<rsup|>|)><rsup|2>|<around*|[|exp<around*|(|4\<beta\>/3|)>-1|]><rsup|2>>=\<cal-O\><around*|(|<frac|128L<rsup|2>|3\<mu\><rsup|2>\<delta\><rsup|2>>e<rsup|<sqrt|log<around*|(|R<rsub|0>/\<varepsilon\>|)>>>|)>
  </equation*>

  <subsection|Improved Analysis for Stability with Sharpness>

  For brevity we assume that <math|f<rsub|z><around*|(|x,\<xi\>|)>=f<around*|(|x,\<xi\>|)>>
  and that the sharpness condition holds with parameter <math|\<mu\>>. Now we
  consider the stability analysis where we are given\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|y|^>>|<cell|=>|<cell|arg
    min<rsub|x\<in\>\<cal-X\>> <around*|{|f<around*|(|x,B|)>+<frac|\<gamma\>|2><around*|\<\|\|\>|x-y|\<\|\|\>><rsup|2>|}>>>|<row|<cell|<wide|y|^><rsub|i>>|<cell|=>|<cell|arg
    min<rsub|x\<in\>\<cal-X\>> <around*|{|f<around*|(|x,B<rsub|<around*|(|i|)>>|)>+<frac|\<gamma\>|2><around*|\<\|\|\>|x-y|\<\|\|\>><rsup|2>|}>.>>>>
  </eqnarray*>

  By sharpness condition we immediately have

  <\equation*>
    <around*|\<\|\|\>|y-x<rsup|\<ast\>>|\<\|\|\>>\<leq\>\<mu\><rsup|-1><around*|[|f<around*|(|y|)>-f<rsup|\<ast\>>|]>
  </equation*>

  and by letting <math|x<rsup|\<ast\>>=Proj<rsub|\<cal-X\>><around*|(|y|)>>,
  we have

  <\equation*>
    f<around*|(|<wide|y|^>,B|)>+<frac|\<gamma\>|2><around*|\<\|\|\>|<wide|y|^>-y|\<\|\|\>><rsup|2>\<leq\>f<around*|(|x<rsup|\<ast\>>,B|)>+<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|\<ast\>>-y|\<\|\|\>><rsup|2>-<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|\<ast\>>-<wide|y|^>|\<\|\|\>><rsup|2>
  </equation*>

  <\equation*>
    f<around*|(|<wide|y|^><rsub|i>,B<rsub|<around*|(|i|)>>|)>+<frac|\<gamma\>|2><around*|\<\|\|\>|<wide|y|^><rsub|i>-y|\<\|\|\>><rsup|2>\<leq\>f<around*|(|x<rsup|\<ast\>>,B<rsub|<around*|(|i|)>>|)>+<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|\<ast\>>-y|\<\|\|\>><rsup|2>-<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|\<ast\>>-<wide|y|^><rsub|i>|\<\|\|\>><rsup|2>.
  </equation*>

  Rearranging the terms, we have

  <\equation*>
    <frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|\<ast\>>-<wide|y|^>|\<\|\|\>><rsup|2>\<leq\>f<around*|(|x<rsup|\<ast\>>,B|)>-f<around*|(|<wide|y|^>,B|)>+<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|\<ast\>>-y|\<\|\|\>><rsup|2>-<frac|\<gamma\>|2><around*|\<\|\|\>|<wide|y|^>-y|\<\|\|\>><rsup|2>
  </equation*>

  and taking expectation gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|\<gamma\>|2>\<bbb-E\><around*|[|<around*|\<\|\|\>|x<rsup|\<ast\>>-<wide|y|^>|\<\|\|\>><rsup|2>|]>>|<cell|\<leq\>>|<cell|\<bbb-E\><around*|[|f<around*|(|x<rsup|\<ast\>>,B|)>-f<around*|(|<wide|y|^>,B|)>|]>+<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|\<ast\>>-y|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|f<around*|(|x<rsup|\<ast\>>|)>-f<around*|(|<wide|y|^>|)>+<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|\<ast\>>-y|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|\<ast\>>-y|\<\|\|\>><rsup|2>>>>>
  </eqnarray*>

  Now we have

  <\equation*>
    \<bbb-E\><around*|[|<around*|\<\|\|\>|x<rsup|\<ast\>>-<wide|y|^>|\<\|\|\>>|]>\<leq\><sqrt|\<bbb-E\><around*|[|<around*|\<\|\|\>|x<rsup|\<ast\>>-<wide|y|^>|\<\|\|\>><rsup|2>|]>>\<leq\><around*|\<\|\|\>|x<rsup|\<ast\>>-y|\<\|\|\>>.
  </equation*>

  On the other hand we have

  <\equation*>
    \<bbb-E\><around*|[|<around*|\<\|\|\>|x<rsup|\<ast\>>-<wide|y|^><rsub|i>|\<\|\|\>>|]>\<leq\><around*|\<\|\|\>|x<rsup|\<ast\>>-y|\<\|\|\>>
  </equation*>

  and by combining the two terms above we have

  <\equation*>
    \<bbb-E\><around*|[|<around*|\<\|\|\>|y<rsub|i>-<wide|y|^><rsub|i>|\<\|\|\>>|]>\<leq\>\<bbb-E\><around*|[|<around*|\<\|\|\>|y<rsub|i>-x<rsup|\<ast\>>|\<\|\|\>>+<around*|\<\|\|\>|<wide|y|^><rsub|i>-x<rsup|\<ast\>>|\<\|\|\>>|]>\<leq\>2<around*|\<\|\|\>|x<rsup|\<ast\>>-y|\<\|\|\>>\<leq\>2\<mu\><rsup|-1><around*|[|f<around*|(|y|)>-f<rsup|\<ast\>>|]>
  </equation*>

  Thus we can conclude that given <math|y=x<rsup|k>>, the stability holds
  with\ 

  <\equation*>
    \<varepsilon\><rsub|k>=min<around*|{|<frac|2L|m\<gamma\>>,2\<mu\><rsup|-1><around*|[|f<around*|(|x<rsup|k>|)>-f<rsup|\<ast\>>|]>|}>.
  </equation*>

  Now we take <math|\<varepsilon\>=2\<mu\><rsup|-1><around*|[|f<around*|(|y|)>-f<rsup|\<ast\>>|]>>
  and return to the original analysis by

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<mu\>\<bbb-E\><rsub|k><around*|[|<around*|\<\|\|\>|x<rsup|k+1>-x<rsup|\<ast\>>|\<\|\|\>>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|\<bbb-E\><rsub|k><around*|[|f<around*|(|x<rsup|k+1>|)>-f<rsup|\<ast\>>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|2<around*|\<\|\|\>|x<rsup|k>-x<rsup|\<ast\>>|\<\|\|\>>+<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|k>-x<rsup|\<ast\>>|\<\|\|\>><rsup|2>-<frac|\<gamma\>|2>\<bbb-E\><rsub|k><around*|[|<around*|\<\|\|\>|x<rsup|k+1>-x<rsup|\<ast\>>|\<\|\|\>><rsup|2>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|2\<mu\><rsup|-1><around*|[|f<around*|(|x<rsup|k>|)>-f<rsup|\<ast\>>|]>+<frac|\<gamma\>|2><around*|\<\|\|\>|x<rsup|k>-x<rsup|\<ast\>>|\<\|\|\>><rsup|2>-<frac|\<gamma\>|2>\<bbb-E\><rsub|k><around*|[|<around*|\<\|\|\>|x<rsup|k+1>-x<rsup|\<ast\>>|\<\|\|\>><rsup|2>|]>>>>>
  </eqnarray*>

  which gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|<frac|1|2><around*|\<\|\|\>|x<rsup|k+1>-x<rsup|\<ast\>>|\<\|\|\>><rsup|2>+<frac|\<mu\>|\<gamma\>><around*|\<\|\|\>|x<rsup|k+1>-x<rsup|\<ast\>>|\<\|\|\>>|]>>|<cell|\<leq\>>|<cell|<frac|1|2><around*|\<\|\|\>|x<rsup|k>-x<rsup|\<ast\>>|\<\|\|\>><rsup|2>+<frac|2|\<gamma\>><around*|\<\|\|\>|x<rsup|k>-x<rsup|\<ast\>>|\<\|\|\>>>>>>
  </eqnarray*>

  or\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|\<mu\>|\<gamma\>>\<bbb-E\><around*|[|<around*|\<\|\|\>|x<rsup|k+1>-x<rsup|\<ast\>>|\<\|\|\>>|]>>|<cell|\<leq\>>|<cell|<wide*|<frac|1|2><around*|\<\|\|\>|x<rsup|k>-x<rsup|\<ast\>>|\<\|\|\>><rsup|2>|\<wide-underbrace\>><rsub|<text|Quadratic>>+<wide*|<frac|2|\<gamma\>><around*|\<\|\|\>|x<rsup|k>-x<rsup|\<ast\>>|\<\|\|\>>|\<wide-underbrace\>><rsub|<text|Linear
    if <math|\<mu\>\<gtr\>2>>>.>>>>
  </eqnarray*>

  \;

  \;

  [1] Bertsekas, Dimitri.<nbsp><with|font-shape|italic|Convex optimization
  algorithms>. Athena Scientific, 2015.

  [2] Davis, D. , D. Drusvyatskiy , and V Charisopoulos. ``Stochastic
  algorithms with geometric step decay converge linearly on sharp
  functions."<nbsp><with|font-shape|italic|arXiv><nbsp>(2019).

  [3] Davis, Damek, et al. ``Subgradient methods for sharp weakly convex
  functions."<nbsp><with|font-shape|italic|Journal of Optimization Theory and
  Applications><nbsp>179.3 (2018): 962-982.

  \;

  \;

  \;
</body>

<initial|<\collection>
</collection>>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|1|?>>
    <associate|auto-3|<tuple|2|?>>
    <associate|auto-4|<tuple|3|?>>
    <associate|auto-5|<tuple|3.1|?>>
    <associate|auto-6|<tuple|3.2|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Literature over
      optimization with sharpness>|<pageref|auto-2>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Literature
      Review> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Preliminaries>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Convex
      Optimization> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Restarting Strategy with
      Decaying Stepsize <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Improved Analysis for
      Stability with Sharpness <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>
    </associate>
  </collection>
</auxiliary>