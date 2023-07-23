<TeXmacs|2.1>

<style|generic>

<\body>
  <\hide-preamble>
    \;

    <assign|X|<macro|<math|<math-bf|X>>>>

    <assign|A|<macro|<math|<math-bf|A>>>>

    <assign|C|<macro|<math|<math-bf|C>>>>

    <assign|b|<macro|<math-bf|b>>>

    <assign|0|<macro|<math-bf|0>>>

    <assign|x|<macro|<math-bf|x>>>

    <assign|y|<macro|<math-bf|y>>>

    <assign|bs|<macro|<math-bf|S>>>

    <assign|c|<macro|<math-bf|c>>>

    <assign|z|<macro|<math|<math-bf|z>>>>

    <assign|I|<macro|<math|<math-bf|I>>>>

    <assign|fa|<macro|\<cal-A\>>>

    <assign|r|<macro|<math-bf|r>>>

    <assign|p|<macro|<math-bf|p>>>

    <assign|M|<macro|<math|<math-bf|M>>>>

    <assign|d|<macro|<math|<math-bf|d>>>>

    <assign|R|<macro|<math|<math-bf|R>>>>

    <assign|m|<macro|<math|<math-bf|m>>>>

    <assign|u|<macro|<math-bf|u>>>

    <assign|s|<macro|<math-bf|s>>>

    <assign|e|<macro|<math-bf|e>>>

    <assign|a|<macro|<math|<math-bf|a>>>>
  </hide-preamble>

  <doc-data|<doc-title|HDSDP: Software for Semidefinite
  Programming>|<doc-author|<author-data|<\author-affiliation>
    Wenzhi Gao<space|1em>Dongdong Ge

    \;

    Shanghai University of Finance and Economics

    \;

    \ Yinyu Ye

    \;

    Stanford University

    \;

    <date|>
  </author-affiliation>>>>

  <abstract-data|<\abstract>
    HDSDP is a numerical software designed for solving the semidefinite
    programming problems. The main framework of HDSDP resembles the
    dual-scaling interior point solver DSDP<cite|benson2008algorithm> and
    several new features, especially a dual algorithm based on the simplified
    homogeneous self-dual embedding, have been implemented. The self-dual
    embedding has several desirable features which enhance numerical
    stability and several new heuristics and computational techniques are
    designed to accelerate its convergence. HDSDP aims to show how
    dual-scaling algorithms benefit from the self-dual embedding and it is
    developed parallel to <verbatim|DSDP5.8>. Numerical experiments over
    several SDP benchmark datasets prove the robustness and efficiency of
    HDSDP.
  </abstract>>

  <section|Introduction>

  Semi-definite programming (SDP) is a mathematical programming problem
  defined by

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<X>>>|<cell|<around*|\<langle\>|<C>,<X>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<fa><X>=<b>>|<cell|>>|<row|<cell|>|<cell|<X>\<in\>\<bbb-S\><rsub|+><rsup|n>,>|<cell|>>>>
  </eqnarray*>

  where we do linear optimization subject to affine constraints over the cone
  of positive-semidefinite matrices. Due to its extensive modeling
  capability, SDP has been employed in various communities including
  combinatorial optimization<cite|goemans1995improved><cite|laurent2005semidefinite>,
  dynamic systems <cite|vandenberghe1996semidefinite>, sums of squares
  optimization<cite|laurent2009sums>, quantum
  information<cite|hayashi2016quantum>, and distance geometry
  <cite|biswas2004semidefinite><cite|so2007theory>. We refer the interested
  readers to <cite|wolkowicz2005semidefinite> for a more comprehensive review
  of SDP applications.

  While SDP is useful in many applications, a fundamental issue is how to
  numerically solve them. Theoretically SDP is a convex conic problem which
  admits efficient polynomial-time algorithms and for general SDPs, the
  interior point method (IPM) is known as a most robust and efficient
  approach. Since the 1990s, high-performance SDPs softwares based on the IPM
  have been developed, including <verbatim|DSDP> <cite|benson2008algorithm>,
  <verbatim|COPT> <cite|copt>, <verbatim|Mosek> <cite|aps2019mosek>,
  <verbatim|Sedumi> <cite|polik2007sedumi>, <verbatim|SDPT3>
  <cite|toh2012implementation>, <verbatim|CSDP> <cite|borchers2006csdp> and
  <verbatim|SDPA> <cite|yamashita2012latest>. While most SDP codes implement
  the IPM, there is also a number of successful attempts adopting other
  algorithms including <cite|kocvara2006pensdp><cite|kwasniewiczimplementation><cite|yang2015sdpnal>
  and a list of such SDP softwares is available in <cite|majumdar2020recent>.

  SDP solvers based on different IPM variants enjoy nice convergence behavior
  both theoretically and practically. Most SDP solvers implement a
  path-following primal-dual approach either using infeasible start
  <cite|potra1998superlinearly> or the homogeneous self-dual embedding
  <cite|potra1998homogeneous> with <verbatim|DSDP> being an exception.
  <verbatim|DSDP>, as its name suggests, implements a dual IPM algorithm
  based on the potential reduction framework proposed in
  <cite|benson1999mixed>. Since the initial release <cite|benson2000solving>,
  <verbatim|DSDP> has gone a long way through several major updates and has
  evolved into an efficient and robust solver for general SDPs
  <cite|benson2008algorithm>. And to further enhance the efficiency and
  robustness of <verbatim|DSDP>, we make another extension by incorporating
  the well-known homogeneous self-dual embedding (HSD) into the dual
  algorithm. This new implementation, named <verbatim|HDSDP(v1.0)>, is
  presented in this manuscript.

  The rest of the manuscript is organized as follows. Section <strong|2>
  describes the SDP formulation of interest and some basic notations. In
  <strong|Section <reference|sec2>>, we review the dual-scaling algorithm for
  SDP and describe how it is combined with the simplified HSD. In
  <strong|Section <reference|sec3>>, we introduce the new infeasibility
  detection and initalization scheme for <verbatim|HDSDP>. Last we present
  the computational results of <verbatim|HDSDP v1.0> on several SDP benchmark
  datasets.

  <section|Formulation and Notations>

  <verbatim|HDSDP> is interested in the standard form primal SDP and its dual
  form

  <center|<math|<tabular|<tformat|<cwith|1|3|2|3|cell-halign|c>|<table|<row|<cell|<around*|(|P|)>>|<cell|min<rsub|<X>>
  >|<cell|<around*|\<langle\>|<C>,<X>|\<rangle\>>>|<cell|>>|<row|<cell|>|<cell|<text|subject
  to>>|<cell|<around*|\<langle\>|<A><rsub|i>,<X>|\<rangle\>>=<b>,>|<cell|i=1,\<ldots\>,m>>|<row|<cell|>|<cell|>|<cell|<X>\<in\>\<bbb-S\><rsub|+><rsup|n>>|<cell|>>>>>><space|2em><math|<tabular|<tformat|<cwith|1|3|2|3|cell-halign|c>|<table|<row|<cell|<around*|(|D|)>>|<cell|max<rsub|<y>,<bs>>>|<cell|<b><rsup|\<top\>><y>>>|<row|<cell|>|<cell|<text|subject
  to>>|<cell|<big|sum><rsub|i=1><rsup|m><A><rsub|i>y<rsub|i>+<bs>=<C>>>|<row|<cell|>|<cell|>|<cell|<bs>\<in\>\<bbb-S\><rsub|+><rsup|n>,>>>>>>>

  where the problem data <math|<A><rsub|i>,<C>> are <math|n\<times\>n>
  symmetric matrices (<math|\<bbb-S\><rsup|n>>) and
  <math|<b>\<in\>\<bbb-R\><rsup|m>> is a real vector. Matrix inner product is
  defined by <math|<around*|\<langle\>|<C>,<X>|\<rangle\>>\<assign\><big|sum><rsub|k,l>c<rsub|i
  j>x<rsub|i j>> and <math|\<bbb-S\><rsub|+><rsup|n>> denotes the cone of
  positive semi-definite matrices. For brevity we use <math|<X>\<succeq\><0>>
  to denote the relation <math|<X>\<in\>\<bbb-S\><rsub|+><rsup|n>> and the
  linear map <math|<fa>:\<bbb-S\><rsup|n>\<rightarrow\>\<bbb-R\><rsup|m>> and
  its adjoint <math|<fa><rsup|\<ast\>>:\<bbb-R\><rsup|m>\<rightarrow\>\<bbb-S\><rsup|n>>
  are respectively defined by <math|<fa><X>\<assign\><around*|(|<around*|\<langle\>|<A>,<X><rsub|1>|\<rangle\>>,\<ldots\>,<around*|\<langle\>|<A><rsub|m>,<X><rsub|m>|\<rangle\>>|)><rsup|\<top\>>>
  and <math|<fa><rsup|\<ast\>><y>\<assign\><big|sum><rsub|i=1><rsup|m><A><rsub|i>y<rsub|i>>.
  <math|<around*|\<\|\|\>|<A>|\<\|\|\>><rsub|F>\<assign\><sqrt|<big|sum><rsub|i
  j>a<rsub|i j><rsup|2>>> denotes matrix Frobenius norm. With the above
  notations we can rewrite the primal and dual problems by

  <center|<math|<tabular|<tformat|<cwith|1|3|2|3|cell-halign|c>|<table|<row|<cell|<around*|(|P|)>>|<cell|min<rsub|<X>>
  >|<cell|<around*|\<langle\>|<C>,<X>|\<rangle\>>>|<cell|>>|<row|<cell|>|<cell|<text|subject
  to>>|<cell|<fa><X>=<b>,>|<cell|i=1,\<ldots\>,m>>|<row|<cell|>|<cell|>|<cell|<X>\<succeq\><0>>|<cell|>>>>>><space|2em><math|<tabular|<tformat|<cwith|1|3|2|3|cell-halign|c>|<table|<row|<cell|<around*|(|D|)>>|<cell|max<rsub|<y>,<bs>>>|<cell|<b><rsup|\<top\>><y>>>|<row|<cell|>|<cell|<text|subject
  to>>|<cell|<fa><rsup|\<ast\>><y>+<bs>=<C>>>|<row|<cell|>|<cell|>|<cell|<bs>\<succeq\><0>.>>>>>>>

  and the feasible region for <math|<around*|(|P|)>> and
  <math|<around*|(|D|)>> are respectively denoted by
  <math|\<cal-F\><around*|(|P|)>\<assign\><around*|{|<X>:<fa><X>=<b>,<X>\<succeq\><0>|}>>
  and <math|\<cal-F\><around*|(|D|)>\<assign\><around*|{|<around*|(|<y>,<bs>|)>:<fa><rsup|\<ast\>><y>+<bs>=<C>,<bs>\<succeq\><0>|}>>.
  Any <math|<X>\<in\>\<cal-F\><around*|(|P|)>> is primal feasible and
  <math|<around*|(|<y>,<bs>|)>\<in\>\<cal-F\><around*|(|D|)>> is dual
  feasible. The interior of the feasible regions are denoted by
  <math|\<cal-F\><rsup|0><around*|(|P|)>,\<cal-F\><rsup|0><around*|(|D|)>>
  and any point in <math|\<cal-F\><rsup|0>> is called an interior point
  solution.

  <verbatim|HDSDP> is an optimization software designed to solve both
  (<em|P>) and (<em|D>) through the simplified HSD and in the next section we
  get down to the underlying theoretical intuitions.

  <\remark>
    Although in this manuscript we focus on the SDP of single-block for ease
    of exposition, a natural generalization to multi-block SDPs applies
    rather straightforward.
  </remark>

  <section|Homogeneous Dual Scaling Algorithm ><label|sec2>

  In this section we first briefly review the dual-scaling algorithm in
  <verbatim|DSDP> and one of its interpretations through Newton's method.
  Then we show how dual-scaling can leverage this interpretation to get
  applied to the simplified HSD embedding.\ 

  <subsection|Dual-scaling Algorithm>

  The dual-scaling algorithm for SDP is initially proposed in
  <cite|benson1999mixed> it works under three conditions. <strong|1)> the
  constriant data <math|<around*|{|<A><rsub|i>|}>> are linearly independent.
  <strong|2)> (<em|P>) and (<em|D>) both admit an interior point solution.
  <strong|3)> an interior dual feasible solution
  <math|<around*|(|<y><rsup|0>,<bs><rsup|0>|)>\<in\>\<cal-F\><rsup|0><around*|(|D|)>>
  is known. The first two conditions imply strong duality and the existence
  of primal-dual optimal pair <math|<around*|(|<X><rsup|\<ast\>>,<y><rsup|\<ast\>>,<bs><rsup|\<ast\>>|)>>
  satisfying complementarity <math|<around*|\<langle\>|<X><rsup|\<ast\>>,<bs><rsup|\<ast\>>|\<rangle\>>=0>
  and <math|<X><bs>=<0>>. Also, we have the existence of the central path
  <math|\<cal-C\><around*|(|\<mu\>|)>\<assign\><around*|{|<around*|(|<X>,<y>,<bs>|)>\<in\>\<cal-F\><rsup|0><around*|(|P|)>\<times\>\<cal-F\><rsup|0><around*|(|D|)>:<X><bs>=\<mu\><I>|}>>
  which serves as the foundation of the path-following IPMs.\ 

  By the last condition, dual-scaling starts from a dual-feasible solution
  <math|<around*|(|<y>,<bs>|)>> and takes Newton's step towards
  <math|<fa><X>=<b>,<fa><rsup|\<ast\>><y>+<bs>=<C>> and
  <math|<X><bs>=\<mu\><I>> by solving

  <\eqnarray*>
    <tformat|<table|<row|<cell|<fa><around*|(|<X>+\<Delta\><X>|)>>|<cell|=>|<cell|<b>>>|<row|<cell|\<cal-A\><rsup|\<ast\>>\<Delta\><y>+\<Delta\><bs>>|<cell|=>|<cell|<0>>>|<row|<cell|\<mu\><bs><rsup|-1>\<Delta\><bs><bs><rsup|-1>+\<Delta\><X>>|<cell|=>|<cell|\<mu\><bs><rsup|-1>-<X>,>>>>
  </eqnarray*>

  where the last relation linearizes <math|<X>=\<mu\><bs><rsup|-1>> instead
  of <math|<X><bs>=\<mu\><I>> and <math|<bs><rsup|-1>> is called a scaling
  matrix. Then <math|<X>> and <math|\<Delta\><X>> vanish in the Schur
  complement system

  <\equation*>
    <matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<A><rsub|1>,<bs><rsup|-1><A><rsub|1><bs><rsup|-1>|\<rangle\>>>|<cell|\<cdots\>>|<cell|<around*|\<langle\>|<A><rsub|1>,<bs><rsup|-1><A><rsub|m><bs><rsup|-1>|\<rangle\>>>>|<row|<cell|\<vdots\>>|<cell|\<ddots\>>|<cell|\<vdots\>>>|<row|<cell|<around*|\<langle\>|<A><rsub|m>,<bs><rsup|-1><A><rsub|1><bs><rsup|-1>|\<rangle\>>>|<cell|\<cdots\>>|<cell|<around*|\<langle\>|<A><rsub|m>,<bs><rsup|-1><A><rsub|m><bs><rsup|-1>|\<rangle\>>>>>>>\<Delta\><y>=<frac|1|\<mu\>><b>-<fa><bs><rsup|-1>
  </equation*>

  and the dual algorithm thereby avoids explicit reference to the primal
  variable <math|<X>>. For brevity we will denote the left-hand side matrix
  by <math|<M>>. When <math|<around*|{|<A><rsub|i>|}>,<C>> are sparse, the
  dual variable <math|<bs>=<C>-\<cal-A\><rsup|\<ast\>><y>> inherits the
  sparsity pattern from the data and it is sometimes much cheaper to iterate
  in the dual space to exploit sparsity. Another desirable feature of
  dual-scaling is the availablity of primal solution by solving a projection
  subproblem at the end of the algorithm. The above desirable properties
  characterize the behavior of dual-scaling.

  However, the nice theoretical properties of the dual algorithm does come
  for free. First, an initial dual feasible solution is needed but obtaining
  such a solution is often as difficult as solving the original problem.
  Second, due to a lack of information from the primal space, the
  dual-scaling has to be property guided to avoid deviating too far away from
  the central path. Last, dual-scaling linearizes the highly nonlinear
  relation <math|<X>=\<mu\><bs><rsup|-1>> and this imposes strict constraint
  on the damping factor towards the Newton's direction. To overcome the
  aforementioned difficulties, <verbatim|DSDP> introduces slack variables
  with big-<math|M> penalty to ensure nonempty interior and a trivial dual
  feasible solution. Morever, <verbatim|DSDP> elegantly adopts a potential
  function to guide the dual iterations. This attempt works well in practice
  and makes <verbatim|DSDP> an efficient general SDP solver. For a complete
  description of the theoretical aspects of <verbatim|DSDP> and its delicate
  implementation, we refer the interested readers to
  <cite|benson2000solving><cite|benson2008algorithm>.\ 

  Although <verbatim|DSDP> proves efficient in real practice, the
  big-<math|M> methods needs a prior estimation of the optimal solution to
  avoid affecting optimality. Also, large penalty might lead to numerical
  difficulties and misclassification of infeasible problems when the problem
  is not well-conditioned. It is natural to ask if there exists an
  alternative to big-<math|M> method to address these issues and in this
  work, we propose to leverage the well-known HSD embedding.

  <subsection|Homogeneous Dual-scaling Algorithm>

  In this section, we introduce the theoretical idea behind <verbatim|HDSDP>.
  HSD embedding is a skew symmetric system whose non-trivial interior point
  solution certificates primal-dual feasibility. Here we adopt the simplified
  HSD embedding <cite|xu1996simplified> and its extension to SDP
  <cite|potra1998homogeneous>.

  <strong|HSD embedding for SDP>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<fa><X>-<b>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|-<fa><rsup|\<ast\>><y>+<C>\<tau\>-<bs>>|<cell|=>|<cell|<0>>>|<row|<cell|<b><rsup|\<top\>><y>-<around*|\<langle\>|<C>,<X>|\<rangle\>>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|<X>,<bs>\<succeq\>0,>|<cell|>|<cell|\<kappa\>,\<tau\>\<geq\>0,>>>>
  </eqnarray*>

  where <math|\<kappa\>,\<tau\>\<geq\>0> are two homogenizing variables for
  infeasibility detection. Given parameter <math|\<mu\>>, we can define
  central path for the embedding by

  <\eqnarray*>
    <tformat|<table|<row|<cell|<fa><X>-<b>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|-<fa><rsup|\<ast\>><y>+<C>\<tau\>-<bs>>|<cell|=>|<cell|<0>>>|<row|<cell|<b><rsup|\<top\>><y>-<around*|\<langle\>|<C>,<X>|\<rangle\>>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|<X><bs>=\<mu\><I>,>|<cell|>|<cell|\<kappa\>\<tau\>=\<mu\>>>|<row|<cell|<X>,<bs>\<succeq\>0,>|<cell|>|<cell|\<kappa\>,\<tau\>\<geq\>0.>>>>
  </eqnarray*>

  Here <math|<around*|(|<y>,<bs>,\<tau\>|)>> are jointly considered as dual
  variables. Given a dual point <math|<around*|(|<y>,<bs>,\<tau\>|)>> such
  that <math|-<fa><rsup|\<ast\>><y>+<C>\<tau\>-<bs>=<R>>, <verbatim|HDSDP>
  chooses a damping factor <math|\<gamma\>\<in\><around*|(|0,1|]>> and takes
  Newton's step towards\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<fa><around*|(|<X>+\<Delta\><X>|)>-<b><around*|(|\<tau\>+\<Delta\>\<tau\>|)>>|<cell|=>|<cell|<0>>>|<row|<cell|-\<cal-A\><rsup|\<ast\>><around*|(|<y>+\<Delta\><y>|)>+<C><around*|(|\<tau\>+\<Delta\>\<tau\>|)>-<around*|(|<bs>+\<Delta\><bs>|)>>|<cell|=>|<cell|-\<gamma\><R>>>|<row|<cell|\<mu\><bs><rsup|-1>\<Delta\><bs><bs><rsup|-1>+\<Delta\><X>>|<cell|=>|<cell|\<mu\><bs><rsup|-1>-<X>,>>|<row|<cell|\<mu\>\<tau\><rsup|-2>\<Delta\>\<tau\>+\<Delta\>\<kappa\>>|<cell|=>|<cell|\<mu\>\<tau\><rsup|-1>-\<kappa\>,>>>>
  </eqnarray*>

  where, as in <verbatim|DSDP>, we linearize <math|<X>=\<mu\><bs><rsup|-1>>
  and <math|\<kappa\>=\<mu\>\<tau\><rsup|-1>>. We note that the damping
  factor can be chosen after the directions are computed and for simplicity
  we set can <math|\<gamma\>=0>. Then a Newton direction
  <math|<around*|(|\<Delta\><y>,\<Delta\>\<tau\>|)>> is computed through the
  following system.\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<matrix|<tformat|<table|<row|<cell|\<mu\><M>>|<cell|-<b>-\<mu\>\<cal-A\><bs><rsup|-1><C><bs><rsup|-1>>>|<row|<cell|-<b>+\<mu\>\<cal-A\><bs><rsup|-1><C><bs><rsup|-1>>|<cell|-\<mu\><around*|(|<around*|\<langle\>|<C>,<bs><rsup|-1><C><bs><rsup|-1>|\<rangle\>>+\<tau\><rsup|-2>|)>>>>>><matrix|<tformat|<table|<row|<cell|\<Delta\><y>>>|<row|<cell|\<Delta\>\<tau\>>>>>>>>|<row|<cell|>|<cell|=>|<cell|<matrix|<tformat|<table|<row|<cell|<b>\<tau\>>>|<row|<cell|<b><rsup|\<top\>><y>-\<mu\>\<tau\><rsup|-1>>>>>>-\<mu\><matrix|<tformat|<table|<row|<cell|<fa><bs><rsup|-1>>>|<row|<cell|<around*|\<langle\>|<C>,<bs><rsup|-1>|\<rangle\>>>>>>>+\<mu\><matrix|<tformat|<table|<row|<cell|<fa><bs><rsup|-1><R><bs><rsup|-1>>>|<row|<cell|<around*|\<langle\>|<C>,<bs><rsup|-1><R><bs><rsup|-1>|\<rangle\>>>>>>>>>>>
  </eqnarray*>

  In practice, <verbatim|HDSDP> solves <math|\<Delta\><y><rsub|1>\<assign\><M><rsup|-1><b>,\<Delta\><y><rsub|2>\<assign\><M><rsup|-1><fa><bs><rsup|-1>,\<Delta\><y><rsub|3>\<assign\><M><rsup|-1><fa><bs><rsup|-1><R><bs><rsup|-1>,\<Delta\><y><rsub|4>=<M><rsup|-1>\<cal-A\><bs><rsup|-1><C><bs><rsup|-1>>,
  plugs <math|\<Delta\><y>=<frac|\<tau\>+\<Delta\>\<tau\>|\<mu\>>\<Delta\><y><rsub|1>-\<Delta\><y><rsub|2>+\<Delta\><y><rsub|3>+\<Delta\><y><rsub|4>\<Delta\>\<tau\>>
  into the second equation to solve for <math|\<Delta\>\<tau\>> and finally
  assembles <math|\<Delta\><y>> using <math|\<tau\>,\<Delta\>\<tau\>> and
  <math|\<mu\>>. Given the above relations,
  <math|<X><around*|(|\<mu\>|)>\<assign\>\<mu\><bs><rsup|-1><around*|(|<bs>-<R>+\<cal-A\><rsup|\<ast\>>\<Delta\><y>-<C>\<Delta\>\<tau\>|)><bs><rsup|-1>>
  satisfies <math|<fa><X><around*|(|\<mu\>|)>-<b><around*|(|\<tau\>+\<Delta\>\<tau\>|)>=<0>>
  and <math|<X><around*|(|\<mu\>|)>\<succeq\><0>> iff.
  <math|-\<cal-A\><rsup|\<ast\>><around*|(|<y>-\<Delta\><y>|)>+<C><around*|(|\<tau\>-\<Delta\>\<tau\>|)>-2<R>\<succeq\><0>>.
  When <math|-\<cal-A\><rsup|\<ast\>><around*|(|<y>-\<Delta\><y>|)>+<C><around*|(|\<tau\>-\<Delta\>\<tau\>|)>-2<R>>
  is positive definite, an objective bound follows by\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|z|\<bar\>>>|<cell|=>|<cell|<around*|\<langle\>|<C>\<tau\>,<X><around*|(|\<mu\>|)>|\<rangle\>>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<langle\>|<R>+<bs>+<fa><rsup|\<ast\>><y>,\<mu\><bs><rsup|-1><around*|(|<bs>-<R>+\<cal-A\><rsup|\<ast\>>\<Delta\><y>-<C>\<Delta\>\<tau\>|)><bs><rsup|-1>|\<rangle\>>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|\<tau\>+\<Delta\>\<tau\>|)><b><rsup|\<top\>><y>+n\<mu\>+<around*|(|<fa><bs><rsup|-1>+<fa><bs><rsup|-1><R><bs><rsup|-1>|)><rsup|\<top\>>\<Delta\><y>+\<mu\><around*|(|<around*|\<langle\>|<C>,<bs><rsup|-1>|\<rangle\>>+<around*|\<langle\>|<C>,<bs><rsup|-1><C><bs><rsup|-1>|\<rangle\>>|)>\<Delta\>\<tau\>>>>>
  </eqnarray*>

  and dividing both sides by <math|\<tau\>>. Alternatively <verbatim|HDSDP>
  extracts a bound from the projection problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<X>>>|<cell|<around*|\<\|\|\>|<bs><rsup|1/2><X><bs><rsup|1/2>-\<mu\><I>|\<\|\|\>><rsub|F><rsup|2>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<fa><X>=<b>\<tau\>,>|<cell|>>>>
  </eqnarray*>

  whose optimal solution <math|<X><rprime|'><around*|(|\<mu\>|)>\<assign\>\<mu\><bs><rsup|-1><around*|(|<C>\<tau\>-<fa><rsup|\<ast\>><around*|(|<y>-\<Delta\><rprime|'><y>|)>-<R>|)><bs><rsup|-1>>
  is given by <math|\<Delta\><y><rprime|'>=<frac|\<tau\>|\<mu\>>\<Delta\><y><rsub|1>-\<Delta\><y><rsub|2>>
  and\ 

  <\equation*>
    z=<around*|\<langle\>|<C>\<tau\>,<X><rprime|'><around*|(|\<mu\>|)>|\<rangle\>>=\<mu\><around*|{|<around*|\<langle\>|<R>,<bs><rsup|-1>|\<rangle\>>+<around*|(|<fa><bs><rsup|-1><R><rsub|<y>><bs><rsup|-1>+<fa><bs><rsup|-1>|)><rsup|\<top\>>\<Delta\><rprime|'><y>+n|}>+\<tau\><b><rsup|\<top\>><y>.
  </equation*>

  In practice, when the explicit solution <math|<X><around*|(|\<mu\>|)>> or
  <math|<X><rprime|'><around*|(|\<mu\>|)>> is needed, we can decompose
  <math|<bs>> and solve two sets of linear systems to obtain them.

  One major computationally intensive part in <verbatim|HDSDP> is the setup
  of the Schur complement matrix <math|<M>>. Since the objective <math|<C>>
  often does not enjoy the same structure as
  <math|<around*|{|<A><rsub|i>|}>>, the additional cost from
  <math|<fa><bs><rsup|-1><C><bs><rsup|-1>> and
  <math|<around*|\<langle\>|<C>,<bs><rsup|-1><C><bs><rsup|-1>|\<rangle\>>>
  calls for more efficient techniques to set up <math|<M>>. In
  <verbatim|HDSDP>, a permutation <math|\<sigma\><around*|(|m|)>> of the rows
  of <math|<M>> is generated heuristically and each row is set up using one
  of the following five techniques.

  <\equation*>
    <matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<A><rsub|\<sigma\><around*|(|1|)>>,<bs><rsup|-1><A><rsub|\<sigma\><around*|(|1|)>><bs><rsup|-1>|\<rangle\>>>|<cell|\<cdots\>>|<cell|<around*|\<langle\>|<A><rsub|\<sigma\><around*|(|1|)>>,<bs><rsup|-1><A><rsub|\<sigma\><around*|(|m|)>><bs><rsup|-1>|\<rangle\>>>>|<row|<cell|>|<cell|\<ddots\>>|<cell|\<vdots\>>>|<row|<cell|>|<cell|>|<cell|<around*|\<langle\>|<A><rsub|\<sigma\><around*|(|m|)>>,<bs><rsup|-1><A><rsub|\<sigma\><around*|(|m|)>><bs><rsup|-1>|\<rangle\>>>>>>>
  </equation*>

  Technique <strong|M1> and <strong|M2> are inherited from <verbatim|DSDP>.
  They exploit the low-rank structure of the problem data and an
  eigen-decomposition of the problem data\ 

  <\equation*>
    <A><rsub|i>=<big|sum><rsub|r=1><rsup|rank<around*|(|<A><rsub|i>|)>>\<lambda\><rsub|i
    r><a><rsub|i,r><a><rsup|\<top\>><rsub|i,r>
  </equation*>

  has to be computed at the beginning of the algorithm.\ 

  <strong|Technique M1>

  <\enumerate>
    <item><strong|Setup> <math|><math|<math-bf|B><rsub|\<sigma\><around*|(|i|)>>=<big|sum><rsub|r=1><rsup|rank<around*|(|<A><rsub|\<sigma\><around*|(|i|)>>|)>>\<lambda\><rsub|r><around*|(|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,r>|)><around*|(|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,r>|)><rsup|\<top\>>>.

    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>
    \<sigma\><around*|(|j|)>>=<around*|\<langle\>|<math-bf|B><rsub|\<sigma\><around*|(|i|)>>,<A><rsub|\<sigma\><around*|(|j|)>>|\<rangle\>>,\<forall\>j\<geq\>i>.
  </enumerate>

  <strong|Technique M2>

  <\enumerate>
    <item><strong|Setup> <math|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,r>,r=1,\<ldots\>,r<rsub|\<sigma\><around*|(|i|)>>>.

    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>\<sigma\><around*|(|j|)>>=<big|sum><rsub|r=1><rsup|rank<around*|(|<A><rsub|\<sigma\><around*|(|i|)>>|)>><around*|(|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,r>|)><rsup|\<top\>><A><rsub|\<sigma\><around*|(|j|)>><around*|(|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,r>|)>>.
  </enumerate>

  Technique <strong|M3>, <strong|M4> and <strong|M5> follow the idea of
  <verbatim|SDPA> and exploit the sparsity of the constraints
  <cite|fujisawa1997exploiting>.

  <strong|Technique M3>

  <\enumerate>
    <item><strong|Setup> <math|<math-bf|B><rsub|\<sigma\><around*|(|i|)>>=<bs><rsup|-1><A><rsub|\<sigma\><around*|(|i|)>><bs><rsup|-1>>.

    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>\<sigma\><around*|(|j|)>>=<around*|\<langle\>|<math-bf|B><rsub|\<sigma\><around*|(|i|)>>,<A><rsub|\<sigma\><around*|(|j|)>>|\<rangle\>>,\<forall\>j\<geq\>i>.
  </enumerate>

  <strong|Technique M4>

  <\enumerate>
    <item><strong|Setup> <math|<math-bf|B><rsub|\<sigma\><around*|(|i|)>>=<bs><rsup|-1><A><rsub|\<sigma\><around*|(|i|)>>>.

    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>\<sigma\><around*|(|j|)>>=<around*|\<langle\>|<math-bf|B><rsub|\<sigma\><around*|(|i|)>><bs><rsup|-1>,<A><rsub|\<sigma\><around*|(|j|)>>|\<rangle\>>,\<forall\>j\<geq\>i>.
  </enumerate>

  <strong|Technique M5>

  <\enumerate>
    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>\<sigma\><around*|(|j|)>>=<around*|\<langle\>|<bs><rsup|-1><A><rsub|\<sigma\><around*|(|i|)>><bs><rsup|-1>,<A><rsub|\<sigma\><around*|(|j|)>>|\<rangle\>>,\<forall\>j\<geq\>i>
    directly.
  </enumerate>

  At the beginning of the algorithm, <verbatim|HDSDP> first estimates the
  number of flops of each technique and each row is associated with the
  cheapest technique that minimizes the overall computation cost. This is a
  major improvement over <verbatim|DSDP> in the computational aspect and we
  give a more detailed description of the implementation in the later
  sections.

  After setting up the Schur matrix, a conjugate gradient (CG) method is
  employed to solve the linear systems. The maximum number of iterations is
  chosen around <math|50/m> and is heuristically adjusted. Either the
  diagonal of <math|<M>> or its Cholesky decomposition is chosen to
  precondition the system and after each Cholesky preconditioner is computed,
  <verbatim|HDSDP> reuses in the consecutive solves till a heuristic
  determines that the current preconditioner needs update. When the algorithm
  approaches optimality, <math|<M>> might become ill-conditioned and an
  <samp|LDLT> decomposition is applied once Cholesky fails.

  Using the Newton's direction, <verbatim|HDSDP> computes the maximum
  stepsize <math|\<alpha\>=max<around*|{|\<alpha\>\<in\><around*|[|0,1|]>:<bs>+\<alpha\>\<Delta\><bs>\<succeq\><0>,\<tau\>+\<alpha\>\<Delta\>\<tau\>\<geq\>0|}>>
  via a Lanczos procedure <cite|toh2002note>. Then to determine a proper
  stepsize, one linesearch determines <math|\<alpha\><rsub|c>> such that the
  barrier satisfies <math|-log det <around*|(|<bs>+\<alpha\><rsub|c>\<Delta\><bs>|)>-log
  det <around*|(|\<tau\>+\<alpha\><rsub|c>\<Delta\>\<tau\>|)>\<leq\>-log det
  <around*|(|<bs>+\<alpha\><rsub|c>\<Delta\><bs>|)>-log det
  <around*|(|\<tau\>+\<alpha\><rsub|c>\<Delta\>\<tau\>|)>> and a second
  linesearch chooses <math|\<gamma\>> such that
  <math|<bs>+\<alpha\><rsub|c>\<cal-A\><rsup|\<ast\>>\<Delta\><y><rsub|2>+\<alpha\><rsub|c>\<gamma\><around*|(|<R>-\<cal-A\><rsup|\<ast\>>\<Delta\><y><rsub|3>|)>\<succeq\><0>>.
  Then a full direction <math|<around*|(|\<Delta\><y><rsup|\<gamma\>>,\<Delta\><bs><rsup|\<gamma\>>,\<Delta\>\<tau\><rsup|\<gamma\>>|)>>
  is finally determined from the Newton system with damping factor
  <math|\<gamma\>> and a third Lanczos procedure computes
  <math|\<alpha\><rprime|'>=max<around*|{|\<alpha\>\<in\><around*|[|0,1|]>:<bs>+\<alpha\>\<Delta\><bs><rsup|r>\<succeq\><0>,\<tau\>+\<alpha\>\<Delta\>\<tau\><rsup|r>\<geq\>0|}>>.
  Last <verbatim|HDSDP> updates <math|<y>\<leftarrow\><y>+0.95\<alpha\><rprime|'>\<Delta\><y><rsup|r>>
  and <math|\<tau\>\<leftarrow\>\<tau\>+0.95\<alpha\>\<Delta\>\<tau\><rsup|r>>.
  To make further use of the Schur matrix and maintain centrality of the
  iterates, the above procedure is repeated several times in an iteration.\ 

  In <verbatim|HDSDP>, the update of the barrier parameter <math|\<mu\>> is
  another critical factor. At the end of each iteration, <verbatim|HDSDP>
  updates the barrier parameter <math|\<mu\>> by
  <math|<around*|(|<wide|z|\<bar\>>-<b><rsup|\<top\>><y>+\<theta\><around*|\<\|\|\>|<R>|\<\|\|\>><rsub|F>|)>/\<rho\>n>,
  where <math|\<rho\>> and <math|\<theta\>> are pre-defined parameters.
  Heuristics also adjust <math|\<mu\>> by
  <math|\<alpha\><rsub|c>,\<alpha\><rprime|'>> and <math|\<gamma\>> computed
  in each iteration.\ 

  To get the best of the two worlds, <verbatim|HDSDP> implements the same
  dual-scaling algorithm as in <verbatim|DSDP5.8>. When the dual
  infeasibility <math|<around*|\<\|\|\>|<fa><rsup|*\<ast\>><y>+<bs>-<C>|\<\|\|\>><rsub|F>\<leq\>\<varepsilon\>\<tau\>>
  and <math|\<mu\>> are sufficiently small, <verbatim|HDSDP> fixs
  <math|\<tau\>=1> and re-starts with <math|<around*|(|<y>/\<tau\>,<bs>/\<tau\>|)>>
  and uses dual-scaling to guide convergence. For the detailed implementation
  of dual-scaling in <verbatim|DSDP5.8> refer to <cite|benson2008algorithm>.
  To sum up, <verbatim|HDSDP> implements strategies and computational tricks
  tailored for the embedding and several of them are designed to mimic the
  behavior of dual-scaling from <verbatim|DSDP5.8>.

  <section|Initalization and Feasibility Certificate>

  Internally <verbatim|HDSDP> deals with the following problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<X>>>|<cell|<around*|\<langle\>|<C>,<X>|\<rangle\>>+<u><rsup|\<top\>><x><rsub|u>+<math-bf|l><rsup|\<top\>><x><rsub|l>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<fa><X>+<x><rsup|u>-<x><rsup|l>=<b>>|<cell|>>|<row|<cell|>|<cell|<X>\<succeq\><0>,<x><rsub|u>\<geq\><0>,<x><rsub|l>\<geq\><0>>|<cell|>>>>
  </eqnarray*>

  and together with its dual, the HSD embedding is\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<fa><X>+<x><rsup|u>-<x><rsup|l>-<b>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|-<fa><rsup|\<ast\>><y>+<C>\<tau\>-<bs>>|<cell|=>|<cell|<0>>>|<row|<cell|<b><rsup|\<top\>><y>-<around*|\<langle\>|<C>,<X>|\<rangle\>>-<u><rsup|\<top\>><x><rsub|u>-<math-bf|l><rsup|\<top\>><x><rsub|l>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|<X>,<bs>\<succeq\>0,>|<cell|>|<cell|\<kappa\>,\<tau\>\<geq\>0,>>>>
  </eqnarray*>

  where the primal problem is relaxed by two slack variables with penalty
  <math|<math-bf|<strong|l>>,<u>> to prevents <math|<y>> goes too large. Due
  to the adoption of HSD embedding, <verbatim|HDSDP> needs no big-<math|M>
  initialization in the dual variable and starts from arbitrary
  <math|<bs>\<succ\><0>,\<tau\>\<gtr\>0> and by default <math|<bs>> is
  initialized by a multiple of the identity matrix.\ 

  One major feature of the <verbatim|HSD> embedding is its capability to
  detect infeasibility and <verbatim|HDSDP> detects dual infeasibility via
  the emedding. Given infeasibility tolerance <math|\<varepsilon\><rsub|f>>,
  <verbatim|HDSDP> classifies the problem as primal unbounded, dual
  infeasible if <math|<around*|\<\|\|\>|<R>|\<\|\|\>><rsub|F>\<gtr\>\<varepsilon\><rsub|f>\<tau\>,\<tau\>/\<kappa\>\<less\>\<varepsilon\><rsub|f>>
  and <math|\<mu\>/\<mu\><rsub|0>\<leq\>\<varepsilon\><rsub|f><rsup|2>>. If
  <math|<R>> is eliminated by <verbatim|HDSDP>, then the embedding
  certificates dual feasibility. If <math|<b><rsup|\<top\>><y>\<gtr\>\<varepsilon\><rsub|f><rsup|-1>>,
  then <verbatim|HDSDP> begins checking the Newton's step
  <math|<around*|(|\<Delta\><y>,\<Delta\><bs>|)>> is a dual improving ray
  every five iterations and once the ray is detected, the problem is
  classified as primal infeasible dual unbounded.

  <section|HDSDP Software><label|sec3>

  <verbatim|HDSDP> is written in <verbatim|ANSI C> and provides a
  self-contained user interface. After the user inputs the data and invokes
  the optimization routine, <verbatim|HDSDP> goes through several modules
  including input check, pre-solving, two-phase optimization, solution
  recovery and post-solving. The modules are implemented independently and
  are sequentially organized in a pipeline by <verbatim|HDSDP>.

  <subsection|Pre-solver>

  One special feature of <verbatim|HDSDP> is the special pre-solving module
  specially designed for SDPs. It inherits several strategies from
  <verbatim|DSDP5.8> and adds several new tricks to work jointly with the
  optimization procedure.\ 

  When the pre-solver is invoked, it first goes through the problem data
  <math|<around*|{|<A><rsub|i>|}>> two rounds to detect the possible low-rank
  structure. The first round uses Gaussian elimination for rank-one structure
  and the second round applies eigenvalue decompostion from
  <verbatim|Lapack>. Two exceptions are when the data matrix is too dense or
  sparse. If a matrix looks too dense to be decomposed efficiently, it is
  skipped and marked as full-rank; if it has very few entries, an elegant
  solution from <verbatim|DSDP5.8> is applied: <strong|1)> a permutation
  gathers the non-zeros to a much smaller dense block. <strong|2)> dense
  eigen routines from <verbatim|Lapack> applies. <strong|3)> the inverse
  permutation recovers the decomposition.\ 

  After detecting the hidden low-rank structures, the pre-solver moves on to
  the analysis of the Schur matrix <math|<M>>: <strong|1)>. the sparsity and
  rank information of the matrices are collected. <strong|2)>. a permutation
  of <math|<around*|{|<A><rsub|i>|}>> is generated in descending order of
  sparsity. <strong|3)>. for each row of <math|<M>>, the flops using each of
  <strong|M1> to <strong|M5> technique is computed and the cheapest technique
  is recorded. The recorded techniques reduces the flops to set up <math|<M>>
  and accelerates the convergence.

  Last the pre-solver scales the coefficients and goes on to look for the
  following structures. <strong|1)>. implied trace: constraints implying
  <math|tr<around*|(|<X>|)>=\<theta\>>. <strong|2)>. implied dual bound:
  constraints implying <math|<math-bf|l>\<leq\><y>\<leq\><u>>. <strong|3)>.
  Empty primal interior: constraints <math|>implying
  <math|tr<around*|(|<X><a><a><rsup|\<top\>>|)>\<approx\>0>. <strong|4)>.
  empty dual interior: constraints implying <math|<A><rsup|\<top\>><y>=<C>>.
  <strong|5)>. feasibility problem. <math|<C>=<0>>. For each of the cases
  above the solver adjusts its internal strategies to enhance the numerical
  stability and convergence behavior.

  <subsection|Two-phase Optimization>

  <verbatim|HDSDP> implements two phase-algorithm which integrates HSD
  embedding (<verbatim|Phase A>) and dual-scaling (<verbatim|Phase B>).
  \ <verbatim|Phase A> targets feasibility certificate and numerical
  stability, while <verbatim|Phase B> aims to drive a dual-feasible solution
  to optimality using potential functition.

  <big-figure|<image|fig.pdf|1000px|||>|Pipeline of <verbatim|HDSDP>>

  \;

  In a word, two phases share the same backend computation routines but are
  associated with different goals and strategies.

  <\subsection>
    Iteration Monitor
  </subsection>

  To help the users capture the progress of the solver, <verbatim|HDSDP>
  prints related information to the screen at different stages of
  optimization.\ 

  <\small>
    <\code>
      --------------------------------------------------------------------------------------------------

      \| Start presolving\ 

      \| - XXX completes in 0.001 seconds\ 

      \| - Matrix statistics ready in 0.000 seconds\ 

      \| \ \ \ Schur Re-ordering: M1: 0 \ M2: 3000 \ M3: 0 \ M4: 1 \ M5: 0\ 

      \| - Special structures found\ 

      \| \ \ \ tr(X) = 3.00e+03 : Bound of X fixed\ 

      \| Presolve Ends. Elapsed Time: 0.371 seconds\ 

      --------------------------------------------------------------------------------------------------

      \| Matrix statistics [Including C]:\ 

      --------------------------------------------------------------------------------------------------

      \| \ \ \ \ \ \ Zero \| \ \ \ \ Sparse \| \ \ \ \ \ Dense \|
      \ \ \ \ Rank-1 \| \ \ \ \ \ \|A\| \ \ \ \ \| \ \ \ \ \ \|b\| \ \ \ \ \|
      \ \ \ \ \ \|C\| \ \ \ \ 

      ----------------------------------------------------\|---------------------------------------------

      \| \ \ \ \ \ \ \ \ \ 0 \| \ \ \ \ \ \ \ \ \ 1 \| \ \ \ \ \ \ \ \ \ 0 \|
      \ \ \ \ \ \ 3000 \| \ 3.00000e+03 \| \ 3.00000e+03 \| \ 1.20000e+04\ 

      --------------------------------------------------------------------------------------------------

      \| Parameter Summary:\ 

      --------------------------------------------------------------------------------------------------

      \| Rhon [1.0, 10.0]: 5\ 

      \| Golden linesearch {0, 1}: 0\ 

      \| Primal relaxation penalty (0.0, inf): 1e+07\ 

      \| Time limit (0.0, inf): 15000\ 

      \| Corrector A: 4 \ Corrector B: 0\ 

      --------------------------------------------------------------------------------------------------

      \| DSDP is initialized with Ry = -1.225e+05 * I
      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      \| DSDP Phase A starts. Eliminating dual infeasibility
      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      --------------------------------------------------------------------------------------------------
    </code>
  </small>

  While pre-solving, <verbatim|HDSDP> prints time statistics when an
  operation is done. Speically, the row

  <\code>
    \ \ \ \ \ \ Schur Re-ordering: M1: 0 \ M2: 3000 \ M3: 0 \ M4: 1 \ M5: 0
  </code>

  indicates how many times each Schur complement technique is applied and if
  special SDP structures are detected, they are also printed to the screen.

  <\code>
    \ \ \ \ \ \ tr(X) = 3.00e+03 : Bound of X fixed.
  </code>

  When the pre-solving ends, <verbatim|HDSDP> prints matrix statistics
  including type and norm. Also, the solver internally adjusts its parameters
  based on the collected information and display them to the user. After
  pre-solving, the solver enters the two-phase optimization procedure and
  invokes HSD embedding (<verbatim|Phase A>).\ 

  <\small>
    <\code>
      Phase A Log: 'P': Primal solution found. '*': Phase A ends. 'F': Error
      happens. 'M': Max iteration

      --------------------------------------------------------------------------------------------------

      \| Iter \| \ \ \ \ \ \ \ \ pObj \| \ \ \ \ \ \ \ \ dObj \|
      \ \ \ \ \ dInf \| \ \ \ \ \ k/t \| \ \ \ \ \ \ mu \| \ \ step \|
      \ \ \ Pnorm \| \ \ E \|

      --------------------------------------------------------------------------------------------------

      \| \ \ \ 1 \| \ 1.00000e+05 \| \ 0.00000e+00 \| \ 6.71e+06 \| 1.00e+00
      \| 4.08e+04 \| \ \ 0.00 \| \ 1.0e+20 \| \ \ \ \ \|

      \| \ \ \ 2 \| \ 2.45010e+08 \| -7.34621e+08 \| \ 0.00e+00 \| 1.00e+00
      \| 3.88e+04 \| \ \ 1.00 \| \ 1.1e+02 \| \ \ * \|

      --------------------------------------------------------------------------------------------------

      Phase A Log Ends.\ 
    </code>
  </small>

  <big-table|<block|<tformat|<table|<row|<cell|<verbatim|Iter>>|<cell|the
  iteration number>>|<row|<cell|<verbatim|pObj>>|<cell|the primal objective
  bound>>|<row|<cell|<verbatim|dObj>>|<cell|the dual objective
  >>|<row|<cell|<verbatim|k/t>>|<cell|<math|\<kappa\>/\<tau\>> from the
  embedding>>|<row|<cell|<verbatim|mu>>|<cell|the current barrier
  parameter>>|<row|<cell|<verbatim|step>>|<cell|stepsize <math|\<alpha\>>
  taken>>|<row|<cell|<verbatim|Pnorm>>|<cell|the proximity to the central
  path>>|<row|<cell|<verbatim|E>>|<cell|event monitor>>>>>|Iteration monitor
  of Phase <verbatim|A>>

  When <verbatim|Phase A> finds a dual-feasible solution, <verbatim|HDSDP>
  utilizes the solution statistic to further adjust the parameters for
  dual-scaling in <verbatim|Phase B>, prints related information and gets
  down to <verbatim|Phase B>.

  <\small>
    <\code>
      --------------------------------------------------------------------------------------------------

      \| DSDP Phase A ends with status: DSDP_PRIMAL_DUAL_FEASIBLE
      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      \| Elapsed Time: 0.535 seconds \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      --------------------------------------------------------------------------------------------------

      \| DSDP Phase A certificates primal-dual feasibility
      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      \| Primal relaxation penalty is set to \ 2.449e+06\ 

      \| Perturbing dual iterations by \ 0.000e+00\ 

      \| DSDP Phase B starts. Restarting dual-scaling
      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      \| Heuristic start: mu: \ 2.177e+04 pObj: \ 2.450e+08 \ dObj:
      -7.346e+08 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      --------------------------------------------------------------------------------------------------
    </code>
  </small>

  The log from <verbatim|Phase B> is similar to <verbatim|Phase A> but
  <verbatim|dInf> is now replaced by <verbatim|pInf> to characterize primal
  infeasibility. Also, <verbatim|k/t> is dropped from log since the embedding
  is not applied.

  <\small>
    <\code>
      Phase B Log: 'P': Primal solution found. '*': Optimal. 'F': Error
      happens. 'M': Max iteration

      --------------------------------------------------------------------------------------------------

      \| Iter \| \ \ \ \ \ \ \ \ \ \ \ \ \ pObj \|
      \ \ \ \ \ \ \ \ \ \ \ \ \ dObj \| \ \ \ \ \ \ pInf \| \ \ \ \ \ \ mu \|
      \ \ step \| \ \ \ Pnorm \| \ \ E \|

      --------------------------------------------------------------------------------------------------

      \| \ \ \ 1 \| \ 2.4500964210e+08 \| -7.3462054382e+08 \| \ 3.001e+03 \|
      2.18e+04 \| \ \ 1.00 \| \ 1.1e+02 \| \ \ \ \ \|

      \| \ \ \ 2 \| \ 2.4500964210e+08 \| -3.6742615954e+07 \| \ 3.001e+03 \|
      2.18e+04 \| \ \ 0.09 \| \ 5.6e+02 \| \ \ \ \ \|

      \| \ \ \ 3 \| \ 1.3061748309e+08 \| -1.8487193898e+06 \| \ 8.915e-03 \|
      3.72e+03 \| \ \ 0.41 \| \ 1.3e+02 \| \ \ P \|

      <text-dots>

      \| \ \ 14 \| -1.1997904329e+04 \| -1.2000000049e+04 \| \ 9.510e-11 \|
      4.66e-06 \| \ \ 0.00 \| \ 9.4e+00 \| \ \ P \|

      \| \ \ 15 \| -1.1999958090e+04 \| -1.2000000002e+04 \| \ 1.902e-12 \|
      9.32e-07 \| \ \ 0.02 \| \ 5.1e+01 \| \ \ F \|

      --------------------------------------------------------------------------------------------------

      Phase B Log Ends.\ 
    </code>
  </small>

  When <verbatim|Phase B> converges, the solver extracts the solution status,
  recovers primal feasible solution (if available) and briefly prints
  solution and time statistics.

  <\small>
    <\code>
      --------------------------------------------------------------------------------------------------

      \| DSDP Phase B ends with status: DSDP_INTERNAL_ERROR
      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      \| Elapsed Time: 5.960 seconds \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      --------------------------------------------------------------------------------------------------

      \| DSDP Ends. \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      --------------------------------------------------------------------------------------------------

      \| Primal solution is extracted. \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 

      \| Final pObj: -1.19993e+04 \ \ dObj: -1.20000e+04\ 

      --------------------------------------------------------------------------------------------------

      \| DSDP Time Summary:\ 

      --------------------------------------------------------------------------------------------------

      \| \ \ \ \ \ \ \ \ \ \ Event \| \ \ \ Time(s) \|\ 

      --------------------------------------------------------------------------------------------------

      \| \ \ \ \ \ \ \ Presolve \| \ \ \ \ \ 0.371 \|\ 

      \| \ Phase A (iter) \| \ \ \ \ \ 0.535 \| (2)\ 

      \| \ Phase B (iter) \| \ \ \ \ \ 5.960 \| (15)\ 

      \| \ \ \ \ \ \ \ \ \ \ Get X \| \ \ \ \ \ 0.809 \|\ 

      \| \ \ \ \ \ \ Postsolve \| \ \ \ \ \ 0.000 \|\ 

      \| \ \ \ \ \ \ \ \ \ \ \ \ All \| \ \ \ \ \ 7.675 \| (17)\ 

      --------------------------------------------------------------------------------------------------
    </code>
  </small>

  <section|Computational results>

  The efficiency of robustness of <verbatim|DSDP5.8> has been proven through
  years of computational experience and <verbatim|HDSDP> aims to achieves
  further improvement on a special class of SDPs where dual-scaling has
  advantage over the primal-dual solvers. In this section, we introduce
  several classes of SDPs suitable for the dual algorithm and we compare the
  performance of <verbatim|HDSDP>, <verbatim|DSDP5.8> and <verbatim|COPT 5.0>
  on several benchmark datasets to verify the performance improvement of
  <verbatim|HDSDP>. For each problem, we only describe the mathematical model
  and refer the readers to <cite|mittelmann2003independent><cite|borchers1999sdplib>
  for the detailed background and formulation. The time statistics in this
  section are obtained using <verbatim|Intel(R) Xeon(R) CPU E5-2640 v4 @
  2.40GHz> with <verbatim|64GB> memory.

  <subsection|Maximum-cut>

  The SDP relaxation of the max-cut problem is represented by

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<X>>>|<cell|<around*|\<langle\>|<C>,<X>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|diag<around*|(|<X>|)>=<math-bf|1>>|<cell|>>|<row|<cell|>|<cell|<X>\<succeq\><0>.>|<cell|>>>>
  </eqnarray*>

  Let <math|<e><rsub|i>> be the <math|i>-th column of the identity matrix and
  the constraint <math|diag<around*|(|<X>|)>=<e>> is decomposed into
  <math|<around*|\<langle\>|<X>,<e><rsub|i><e><rsub|i><rsup|\<top\>>|\<rangle\>>=1,i=1,\<ldots\>,n>.
  Note that <math|<e><rsub|i><e><rsub|i><rsup|\<top\>>> is rank-one and has
  only one non-zero entry, <strong|M2> and <strong|M5> can greatly reduce the
  computation of the Schur matrix.

  <big-table|<block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<cwith|4|16|5|5|cell-halign|c>|<table|<row|<cell|Instance>|<cell|<verbatim|HDSDP>>|<cell|<verbatim|DSDP5.8>>|<cell|<verbatim|COPT
  v5.0>>|<cell|Instance>|<cell|<verbatim|HDSDP>>|<cell|<verbatim|DSDP5.8>>|<cell|<verbatim|COPT
  v5.0>>>|<row|<cell|<verbatim|mcp100>>|<cell|0.03>|<cell|0.03>|<cell|0.12>|<cell|<verbatim|maxG51>>|<cell|1.46>|<cell|4.65>|<cell|5.97>>|<row|<cell|<verbatim|mcp124-1>>|<cell|0.05>|<cell|0.03>|<cell|0.15>|<cell|<verbatim|maxG55>>|<cell|146.92>|<cell|606.66>|<cell|520.05>>|<row|<cell|<verbatim|mcp124-2>>|<cell|0.05>|<cell|0.02>|<cell|0.17>|<cell|<verbatim|maxG60>>|<cell|323.62>|<cell|1485.00>|<cell|1269.83>>|<row|<cell|<verbatim|mcp124-3>>|<cell|0.05>|<cell|0.04>|<cell|0.16>|<cell|<verbatim|G40_mb>>|<cell|12.76>|<cell|17.08>|<cell|25.76>>|<row|<cell|<verbatim|mcp124-4>>|<cell|0.06>|<cell|0.05>|<cell|0.15>|<cell|<verbatim|G40_mc>>|<cell|8.09>|<cell|38.54>|<cell|49.07>>|<row|<cell|<verbatim|mcp250-1>>|<cell|0.15>|<cell|0.10>|<cell|0.66>|<cell|<verbatim|G48_mb>>|<cell|16.70>|<cell|29.71>|<cell|50.49>>|<row|<cell|<verbatim|mcp250-2>>|<cell|0.11>|<cell|0.15>|<cell|0.66>|<cell|<verbatim|G48mc>>|<cell|4.42>|<cell|18.10>|<cell|33.19>>|<row|<cell|<verbatim|mcp250-3>>|<cell|0.17>|<cell|0.21>|<cell|0.62>|<cell|<verbatim|G55mc>>|<cell|118.11>|<cell|344.80>|<cell|505.18>>|<row|<cell|<verbatim|mcp250-4>>|<cell|0.19>|<cell|0.29>|<cell|0.69>|<cell|<verbatim|G59mc>>|<cell|181.20>|<cell|774.70>|<cell|727.89>>|<row|<cell|<verbatim|mcp500-1>>|<cell|0.60>|<cell|0.32>|<cell|0.48>|<cell|<verbatim|G60_mb>>|<cell|261.50>|<cell|650.00>|<cell|964.20>>|<row|<cell|<verbatim|mcp500-2>>|<cell|0.69>|<cell|0.74>|<cell|0.72>|<cell|<verbatim|G60mc>>|<cell|257.10>|<cell|600.08>|<cell|962.79>>|<row|<cell|<verbatim|mcp500-3>>|<cell|0.82>|<cell|1.12>|<cell|1.11>|<cell|<verbatim|torusg3-8>>|<cell|0.85>|<cell|1.42>|<cell|1.04>>|<row|<cell|<verbatim|mcp500-4>>|<cell|1.11>|<cell|1.84>|<cell|2.35>|<cell|<verbatim|torusg3-15>>|<cell|23.77>|<cell|178.8>|<cell|137.60>>|<row|<cell|<verbatim|maxG11>>|<cell|1.27>|<cell|0.82>|<cell|1.07>|<cell|<verbatim|toruspm3-8-50>>|<cell|0.76>|<cell|0.93>|<cell|0.76>>|<row|<cell|<verbatim|maxG32>>|<cell|5.13>|<cell|8.57>|<cell|10.77>|<cell|<verbatim|toruspm3-15-50>>|<cell|19.27>|<cell|91.67>|<cell|117.92>>>>>|Max-cut
  problems>

  <subsection|Graph Partitioning>

  The SDP relaxation of the graph partitioning problem is given by\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<X>>>|<cell|<around*|\<langle\>|<C>,<X>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|diag<around*|(|<X>|)>=<math-bf|1>>|<cell|>>|<row|<cell|>|<cell|<around*|\<langle\>|<math-bf|<strong|1>1><rsup|\<top\>>,<X>|\<rangle\>>=\<beta\>>|<cell|>>|<row|<cell|>|<cell|k<X>-<math-bf|<strong|1>1><rsup|\<top\>>\<succeq\><0>>|<cell|>>|<row|<cell|>|<cell|<X>\<geq\><0>,>|<cell|>>>>
  </eqnarray*>

  where <math|<math-bf|1>> denotes the all-one vector and <math|k,\<beta\>>
  are problem parameters. Although the dual <math|<bs>> no longer enjoys
  sparsity, the low-rank structure is still available to accelerate
  convergence.\ 

  <big-table|<block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<table|<row|<cell|Instance>|<cell|<verbatim|HDSDP>>|<cell|<verbatim|DSDP5.8>>|<cell|<verbatim|COPT
  v5.0>>|<cell|Instance>|<cell|<verbatim|HDSDP>>|<cell|<verbatim|DSDP5.8>>|<cell|<verbatim|COPT
  v5.0>>>|<row|<cell|<verbatim|gpp100>>|<cell|0.03>|<cell|0.04>|<cell|0.19>|<cell|<verbatim|gpp250-4>>|<cell|0.19>|<cell|0.16>|<cell|1.21>>|<row|<cell|<verbatim|gpp124-1>>|<cell|0.04>|<cell|0.09>|<cell|0.28>|<cell|<verbatim|gpp500-1>>|<cell|0.60>|<cell|0.63>|<cell|0.64>>|<row|<cell|<verbatim|gpp124-2>>|<cell|0.05>|<cell|0.05>|<cell|0.28>|<cell|<verbatim|gpp500-2>>|<cell|0.69>|<cell|0.60>|<cell|0.56>>|<row|<cell|<verbatim|gpp124-3>>|<cell|0.05>|<cell|0.05>|<cell|0.34>|<cell|<verbatim|gpp500-3>>|<cell|0.82>|<cell|0.55>|<cell|0.56>>|<row|<cell|<verbatim|gpp124-4>>|<cell|0.06>|<cell|0.05>|<cell|0.36>|<cell|<verbatim|gpp500-4>>|<cell|1.11>|<cell|0.58>|<cell|0.51>>|<row|<cell|<verbatim|gpp250-1>>|<cell|0.15>|<cell|0.16>|<cell|1.58>|<cell|<verbatim|bm1>>|<cell|2.48>|<cell|2.36>|<cell|1.74>>|<row|<cell|<verbatim|gpp250-2>>|<cell|0.11>|<cell|0.15>|<cell|1.33>|<cell|<verbatim|biomedP>>|<cell|221.2>|<cell|<verbatim|Failed>>|<cell|<verbatim|Failed>>>|<row|<cell|<verbatim|gpp250-3>>|<cell|0.17>|<cell|0.14>|<cell|1.46>|<cell|<verbatim|industry2>>|<cell|<verbatim|Failed>>|<cell|<verbatim|Failed>>|<cell|<verbatim|Failed>>>>>>|Graph
  partitioning problems>

  <subsection|Optimal Diagonal Pre-conditioning>

  The optimal diagonal pre-conditioning problem originates from
  <cite|qu2020diagonal>, where given a matrix <math|<math-bf|B>\<succ\><0>>,
  finding a diagonal matrix <math|<math-bf|D>> to minimize the condition
  number <math|\<kappa\><around*|(|<math-bf|D><rsup|-1/2><M><math-bf|D><rsup|-1/2>|)>>
  can be modeled as an SDP. The formulation of optimal diagonal
  pre-conditioning is given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|max<rsub|\<tau\>,<math-bf|D>>>|<cell|\<tau\>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<math-bf|D>\<preceq\><math-bf|B>>|<cell|>>|<row|<cell|>|<cell|\<tau\><math-bf|B>-<math-bf|D>\<preceq\><0>.>|<cell|>>>>
  </eqnarray*>

  Expressing <math|<math-bf|D>=<big|sum><rsub|i><e><rsub|i><e><rsub|i><rsup|\<top\>>d<rsub|i>>,
  the problem is exactly in the SDP dual form. If <math|<math-bf|B>> is also
  sparse, the problem can be efficiently solved using the dual method.

  <big-table|<block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<table|<row|<cell|Instance>|<cell|<verbatim|HDSDP>>|<cell|<verbatim|DSDP5.8>>|<cell|<verbatim|COPT
  v5.0>>>|<row|<cell|<verbatim|diag-bench-500-0.1>>|<cell|6.70>|<cell|7.50>|<cell|6.80>>|<row|<cell|<verbatim|diag-bench-1000-0.01>>|<cell|44.50>|<cell|55.20>|<cell|53.90>>|<row|<cell|<verbatim|diag-bench-2000-0.05>>|<cell|34.30>|<cell|340.70>|<cell|307.10>>|<row|<cell|<verbatim|diag-bench-west0989>>|<cell|6.72>|<cell|113.23>|<cell|108.20>>|<row|<cell|<verbatim|diag-bench-DK01R>>|<cell|13.18>|<cell|<verbatim|Failed>>|<cell|<verbatim|Failed>>>|<row|<cell|<verbatim|diag-bench-micromass_10NN>>|<cell|9.35>|<cell|60.127>|<cell|51.71>>>>>|Optimal
  diagonal pre-conditioning problems>

  <\remark>
    In the optimal pre-conditioning experiment, we start <verbatim|HDSDP>
    from a non-default trivial dual feasible solution
    <math|\<tau\>=-10<rsup|6>,<math-bf|D>=<math-bf|0>>.
  </remark>

  <subsection|Other Problems>

  In this section, we present some other benchmark datasets that are suitable
  for <verbatim|HDSDP>. They enjoy at least one of sparsity and low-rank
  structure.

  <big-table|<block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<table|<row|<cell|Instance>|<cell|Background>|<cell|Feature>|<cell|<verbatim|HDSDP>>|<cell|<verbatim|DSDP5.8>>|<cell|<verbatim|COPT
  v5.0>>>|<row|<cell|<verbatim|checker_1.5>>|<cell|unknown>|<cell|sparse,
  low-rank>|<cell|55.38>|<cell|146.80>|<cell|137.15>>|<row|<cell|<verbatim|foot>>|<cell|unknown>|<cell|sparse,
  low-rank>|<cell|25.76>|<cell|23.83>|<cell|262.54>>|<row|<cell|<verbatim|hand>>|<cell|unknown>|<cell|low-rank>|<cell|5.60>|<cell|5.57>|<cell|50.70>>|<row|<cell|<verbatim|ice_2.0>>|<cell|unknown>|<cell|low-rank>|<cell|706.16>|<cell|1106.00>|<cell|1542.96>>|<row|<cell|<verbatim|p_auss2_3.0>>|<cell|unknown>|<cell|sparse,
  low-rank>|<cell|739.40>|<cell|1066.00>|<cell|1111.72>>|<row|<cell|<verbatim|rendl1_2000_1e-6>>|<cell|unknown>|<cell|low-rank>|<cell|15.74>|<cell|22.31>|<cell|231.80>>|<row|<cell|<verbatim|trto3>>|<cell|topology
  design>|<cell|sparse, low-rank>|<cell|1.67>|<cell|1.73>|<cell|3.04>>|<row|<cell|<verbatim|trto4>>|<cell|topology
  design>|<cell|sparse, low-rank>|<cell|12.28>|<cell|12.40>|<cell|30.09>>|<row|<cell|<verbatim|trto5>>|<cell|topology
  design>|<cell|sparse, low-rank>|<cell|111.59>|<cell|151.0>|<cell|233.16>>|<row|<cell|<verbatim|sensor_500b>>|<cell|sensor
  network localization>|<cell|sparse, low-rank>|<cell|86.73>|<cell|37.87>|<cell|7.01>>|<row|<cell|<verbatim|sensor_1000b>>|<cell|sensor
  network localization>|<cell|sparse, low-rank>|<cell|232.05>|<cell|143.32>|<cell|32.27>>>>>|Feature
  of several benchmark problems>

  <subsection|Mittelmann's Benchmark Test>

  So far <verbatim|HDSDP> is tested and tuned over a large set of benchmarks
  including <verbatim|SDPLIB> <cite|borchers1999sdplib> and Hans Mittelmann's
  sparse SDP benchmark <cite|mittelmann2003independent>. We report the
  solution accuracy and CPU time of <verbatim|HDSDP> on Mittlelmann's
  benchmark in the appendix. Readers can refer to
  <cite|mittelmann2003independent> for a detailed explanation of the error
  measures and the criterion of a successful solve. Among all of 75 problems,
  71 are successfully solved; 3 problems fail due to insufficient memory and
  1 one fails due to error in primal solution recovery. The benchmark test
  proves the efficiency and robustness of <verbatim|HDSDP> as a general
  purpose SDP solver.

  <subsection|When to use DSDP/HDSDP>

  While <verbatim|HDSDP> is designed for general SDPs, it targets the
  problems more tractable in the dual form than by the primal-dual methods.
  Here are some rules in mind when deciding whether to use the dual-scaling
  method (or <verbatim|HDSDP>).

  <\enumerate>
    <item>Does the problem enjoys nice dual structure?

    Many combinatorical problems have formulations friendly to the dual
    methods. Some typical features include (aggregated) sparsity, low-rank
    structure. Dual methods can exploit these features by iterating in the
    dual space and using efficient computational tricks. If the problem is
    dense and most constraints are full-rank, dual method has no advantage of
    the primal-dual solvers due to <strong|1)> comparable iteration cost to
    primal-dual methods. <strong|2)> more iterations for convergence.

    <item>Do we need the primal optimal solution or just the optimal value?

    For some applications dual method fails to recover a correct primal
    solution due to numerical difficulties. If the optimal value is
    sufficient, there is no problem. But if an accurate primal optimal
    solution is always necessary, it is better to be more careful and to test
    the recovery procedure in case of failure at the last step.

    <item>Do we need to certificate infeasibility strictly?

    One weakness of the dual method is the difficulty in infeasibility
    certificate. Although on the dual side this issue is addressed by
    <verbatim|HDSDP> using the embedding, dual methods still suffer from
    failure to identify primal infeasibility.\ 

    <item>Is dual-feasibility hard to attain?

    The first phase of <verbatim|HDSDP> adopts the infeasible Newton's method
    and focuses on eliminating the dual infeasibility. This principle works
    well if the dual constraints are relatively easy to satisfy, but if this
    condition fails to hold (e.g., empty dual interior), the embedding will
    spend a long time deciding feasibility. In this case it is suggested
    using <verbatim|DSDP5.8> or supply an initial dual solution.
  </enumerate>

  To conclude, <verbatim|HDSDP> solves SDPs but it solves <em|certain> SDPs
  <em|efficiently>.

  <section|Conclusions>

  In this paper we propose an extension of the dual-scaling algorithm based
  on the HSD embedding. The resultant solver, <verbatim|HDSDP>, is presented
  to demonstrate how dual method can be effectively integrated with the
  embedding idea. <verbatim|HDSDP> is developed in parallel to
  <verbatim|DSDP5.8> and is entailed with several newly added computational
  techniques. The solver exhibits promising performance on several benchmark
  datasets and is under active development.

  <section|Acknowledgement>

  We thank Dr. Huangfu Qi from <verbatim|COPT> development team for his
  constructive ideas in the solver design and implementation. We also
  appreciate Hans Mittelmann's efforts in benchmarking the solver. Finally,
  we sincerely repect the developers of <verbatim|DSDP> for their precious
  suggestions <cite|benson2008algorithm> and their invaluable efforts getting
  <verbatim|DSDP> through all along the way. It is the efficient and elegant
  implementation from <verbatim|DSDP5.8> that guides <verbatim|HDSDP> to
  where it is.

  <\bibliography*|bib|tm-plain|sdpref|References>
    <\bib-list|29>
      <bibitem*|1><label|bib-aps2019mosek>Mosek ApS. <newblock>Mosek
      optimization toolbox for matlab. <newblock><with|font-shape|italic|User's
      Guide and Reference Manual, Version>, 4, 2019.<newblock>

      <bibitem*|2><label|bib-benson2008algorithm>Steven<nbsp>J
      Benson<localize| and >Yinyu Ye. <newblock>Algorithm 875:
      dsdp5\Vsoftware for semidefinite programming.
      <newblock><with|font-shape|italic|ACM Transactions on Mathematical
      Software (TOMS)>, 34(3):1\U20, 2008.<newblock>

      <bibitem*|3><label|bib-benson2000solving>Steven<nbsp>J Benson, Yinyu
      Ye<localize|, and >Xiong Zhang. <newblock>Solving large-scale sparse
      semidefinite programs for combinatorial optimization.
      <newblock><with|font-shape|italic|SIAM Journal on Optimization>,
      10(2):443\U461, 2000.<newblock>

      <bibitem*|4><label|bib-benson1999mixed>Steven<nbsp>J Benson, Yinyu
      Yeb<localize|, and >Xiong Zhang. <newblock>Mixed linear and
      semidefinite programming for combinatorial and quadratic optimization.
      <newblock><with|font-shape|italic|Optimization Methods and Software>,
      11(1-4):515\U544, 1999.<newblock>

      <bibitem*|5><label|bib-biswas2004semidefinite>Pratik Biswas<localize|
      and >Yinyu Ye. <newblock>Semidefinite programming for ad hoc wireless
      sensor network localization. <newblock><localize|In
      ><with|font-shape|italic|Proceedings of the 3rd international symposium
      on Information processing in sensor networks>, <localize|pages >46\U54.
      2004.<newblock>

      <bibitem*|6><label|bib-borchers1999sdplib>Brian Borchers.
      <newblock>Sdplib 1.2, a library of semidefinite programming test
      problems. <newblock><with|font-shape|italic|Optimization Methods and
      Software>, 11(1-4):683\U690, 1999.<newblock>

      <bibitem*|7><label|bib-borchers2006csdp>Brian Borchers. <newblock>Csdp
      user's guide. <newblock>2006.<newblock>

      <bibitem*|8><label|bib-copt>Cardinal optimizer.<newblock>

      <bibitem*|9><label|bib-fujisawa1997exploiting>Katsuki Fujisawa,
      Masakazu Kojima<localize|, and >Kazuhide Nakata. <newblock>Exploiting
      sparsity in primal-dual interior-point methods for semidefinite
      programming. <newblock><with|font-shape|italic|Mathematical
      Programming>, 79(1):235\U253, 1997.<newblock>

      <bibitem*|10><label|bib-goemans1995improved>Michel<nbsp>X
      Goemans<localize| and >David<nbsp>P Williamson. <newblock>Improved
      approximation algorithms for maximum cut and satisfiability problems
      using semidefinite programming. <newblock><with|font-shape|italic|Journal
      of the ACM (JACM)>, 42(6):1115\U1145, 1995.<newblock>

      <bibitem*|11><label|bib-hayashi2016quantum>Masahito Hayashi.
      <newblock><with|font-shape|italic|Quantum information theory>.
      <newblock>Springer, 2016.<newblock>

      <bibitem*|12><label|bib-kocvara2006pensdp>Michal Kocvara, Michael
      Stingl<localize|, and >PENOPT GbR. <newblock>Pensdp users guide
      (version 2.2). <newblock><with|font-shape|italic|PENOPT GbR>,
      1435:1436, 2006.<newblock>

      <bibitem*|13><label|bib-kwasniewiczimplementation>Tomasz
      Kwasniewicz<localize| and >François Glineur. <newblock>Implementation
      of a semidefinite optimization solver in the julia programming
      language.<newblock>

      <bibitem*|14><label|bib-laurent2009sums>Monique Laurent. <newblock>Sums
      of squares, moment matrices and optimization over polynomials.
      <newblock><localize|In ><with|font-shape|italic|Emerging applications
      of algebraic geometry>, <localize|pages >157\U270. Springer,
      2009.<newblock>

      <bibitem*|15><label|bib-laurent2005semidefinite>Monique
      Laurent<localize| and >Franz Rendl. <newblock>Semidefinite programming
      and integer programming. <newblock><with|font-shape|italic|Handbooks in
      Operations Research and Management Science>, 12:393\U514,
      2005.<newblock>

      <bibitem*|16><label|bib-majumdar2020recent>Anirudha Majumdar, Georgina
      Hall<localize|, and >Amir<nbsp>Ali Ahmadi. <newblock>Recent scalability
      improvements for semidefinite programming with applications in machine
      learning, control, and robotics. <newblock><with|font-shape|italic|Annual
      Review of Control, Robotics, and Autonomous Systems>, 3:331\U360,
      2020.<newblock>

      <bibitem*|17><label|bib-mittelmann2003independent>Hans<nbsp>D
      Mittelmann. <newblock>An independent benchmarking of sdp and socp
      solvers. <newblock><with|font-shape|italic|Mathematical Programming>,
      95(2):407\U430, 2003.<newblock>

      <bibitem*|18><label|bib-polik2007sedumi>Imre Polik, Tamas
      Terlaky<localize|, and >Yuriy Zinchenko. <newblock>Sedumi: a package
      for conic optimization. <newblock><localize|In
      ><with|font-shape|italic|IMA workshop on Optimization and Control,
      Univ. Minnesota, Minneapolis>. Citeseer, 2007.<newblock>

      <bibitem*|19><label|bib-potra1998superlinearly>Florian<nbsp>A
      Potra<localize| and >Rongqin Sheng. <newblock>A superlinearly
      convergent primal-dual infeasible-interior-point algorithm for
      semidefinite programming. <newblock><with|font-shape|italic|SIAM
      Journal on Optimization>, 8(4):1007\U1028, 1998.<newblock>

      <bibitem*|20><label|bib-potra1998homogeneous>Florian<nbsp>A
      Potra<localize| and >Rongqin Sheng. <newblock>On homogeneous
      interrior-point algorithms for semidefinite programming.
      <newblock><with|font-shape|italic|Optimization Methods and Software>,
      9(1-3):161\U184, 1998.<newblock>

      <bibitem*|21><label|bib-qu2020diagonal>Zhaonan Qu, Yinyu Ye<localize|,
      and >Zhengyuan Zhou. <newblock>Diagonal preconditioning: theory and
      algorithms. <newblock><with|font-shape|italic|ArXiv preprint
      arXiv:2003.07545>, 2020.<newblock>

      <bibitem*|22><label|bib-so2007theory>Anthony<nbsp>Man-Cho So<localize|
      and >Yinyu Ye. <newblock>Theory of semidefinite programming for sensor
      network localization. <newblock><with|font-shape|italic|Mathematical
      Programming>, 109(2):367\U384, 2007.<newblock>

      <bibitem*|23><label|bib-toh2002note>Kim-Chuan Toh. <newblock>A note on
      the calculation of step-lengths in interior-point methods for
      semidefinite programming. <newblock><with|font-shape|italic|Computational
      Optimization and Applications>, 21(3):301\U310, 2002.<newblock>

      <bibitem*|24><label|bib-toh2012implementation>Kim-Chuan Toh,
      Michael<nbsp>J Todd<localize|, and >Reha<nbsp>H Tütüncü. <newblock>On
      the implementation and usage of sdpt3\Ua matlab software package for
      semidefinite-quadratic-linear programming, version 4.0.
      <newblock><localize|In ><with|font-shape|italic|Handbook on
      semidefinite, conic and polynomial optimization>, <localize|pages
      >715\U754. Springer, 2012.<newblock>

      <bibitem*|25><label|bib-vandenberghe1996semidefinite>Lieven
      Vandenberghe<localize| and >Stephen Boyd. <newblock>Semidefinite
      programming. <newblock><with|font-shape|italic|SIAM review>,
      38(1):49\U95, 1996.<newblock>

      <bibitem*|26><label|bib-wolkowicz2005semidefinite>Henry Wolkowicz.
      <newblock>Semidefinite and cone programming bibliography/comments.
      <newblock><with|font-shape|italic|Http://orion. uwaterloo. ca/~
      hwolkowi/henry/book/fronthandbk. d/sdpbibliog. pdf>, 2005.<newblock>

      <bibitem*|27><label|bib-xu1996simplified>Xiaojie Xu, Pi-Fang
      Hung<localize|, and >Yinyu Ye. <newblock>A simplified homogeneous and
      self-dual linear programming algorithm and its implementation.
      <newblock><with|font-shape|italic|Annals of Operations Research>,
      62(1):151\U171, 1996.<newblock>

      <bibitem*|28><label|bib-yamashita2012latest>Makoto Yamashita, Katsuki
      Fujisawa, Mituhiro Fukuda, Kazuhiro Kobayashi, Kazuhide
      Nakata<localize|, and >Maho Nakata. <newblock>Latest developments in
      the sdpa family for solving large-scale sdps. <newblock><localize|In
      ><with|font-shape|italic|Handbook on semidefinite, conic and polynomial
      optimization>, <localize|pages >687\U713. Springer, 2012.<newblock>

      <bibitem*|29><label|bib-yang2015sdpnal>Liuqin Yang, Defeng
      Sun<localize|, and >Kim-Chuan Toh. <newblock>Sdpnal <math|+>+: a
      majorized semismooth newton-cg augmented lagrangian method for
      semidefinite programming with nonnegative constraints.
      <newblock><with|font-shape|italic|Mathematical Programming
      Computation>, 7(3):331\U366, 2015.<newblock>
    </bib-list>
  </bibliography*>

  <appendix|Mittlelmann's Benchmark Test>

  The following test of the benchmark dataset is run on an <verbatim|intel
  i11700K> with 128GB memory.

  <\tiny>
    <\with|par-columns|1>
      <\big-table>
        <\with|par-mode|center>
          <tabular*|<tformat|<cwith|1|-1|1|1|cell-lborder|1ln>|<cwith|1|-1|1|1|cell-halign|l>|<cwith|1|-1|1|1|cell-rborder|1ln>|<cwith|1|-1|2|2|cell-halign|r>|<cwith|1|-1|2|2|cell-rborder|1ln>|<cwith|1|-1|3|3|cell-halign|r>|<cwith|1|-1|3|3|cell-rborder|1ln>|<cwith|1|-1|4|4|cell-halign|r>|<cwith|1|-1|4|4|cell-rborder|1ln>|<cwith|1|-1|5|5|cell-halign|r>|<cwith|1|-1|5|5|cell-rborder|1ln>|<cwith|1|-1|6|6|cell-halign|r>|<cwith|1|-1|6|6|cell-rborder|1ln>|<cwith|1|-1|7|7|cell-halign|r>|<cwith|1|-1|7|7|cell-rborder|1ln>|<cwith|1|-1|8|8|cell-halign|r>|<cwith|1|-1|8|8|cell-rborder|1ln>|<cwith|1|-1|1|-1|cell-valign|c>|<cwith|1|1|1|-1|cell-tborder|1ln>|<cwith|1|1|1|-1|cell-bborder|1ln>|<cwith|76|76|1|-1|cell-bborder|1ln>|<table|<row|<cell|Instance>|<cell|Error
          1>|<cell|Error 2>|<cell|Error 3>|<cell|Error 4>|<cell|Error
          5>|<cell|Error 6>|<cell|Time>>|<row|<cell|1dc.1024>|<cell|3.440000e-08>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|5.400000e-06>|<cell|5.600000e-06>|<cell|2500.186>>|<row|<cell|1et.1024>|<cell|7.660000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|3.000000e-06>|<cell|2.900000e-06>|<cell|175.211>>|<row|<cell|1tc.1024>|<cell|6.790000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|2.600000e-06>|<cell|2.600000e-06>|<cell|142.422>>|<row|<cell|1zc.1024>|<cell|1.230000e-08>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.400000e-06>|<cell|1.300000e-06>|<cell|469.066>>|<row|<cell|AlH>|<cell|1.510000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.100000e-04>|<cell|1.100000e-04>|<cell|10162.681>>|<row|<cell|BH2>|<cell|2.850000e-11>|<cell|0.00000e+00>|<cell|3.800000e-09>|<cell|0.00000e+00>|<cell|1.300000e-07>|<cell|3.900000e-07>|<cell|315.241>>|<row|<cell|Bex2_1_5>|<cell|1.130000e-08>|<cell|0.00000e+00>|<cell|3.300000e-12>|<cell|0.00000e+00>|<cell|5.000000e-07>|<cell|1.200000e-07>|<cell|273.331>>|<row|<cell|Bst_jcbpaf2>|<cell|1.410000e-12>|<cell|0.00000e+00>|<cell|5.000000e-11>|<cell|0.00000e+00>|<cell|2.800000e-07>|<cell|1.700000e-07>|<cell|419.132>>|<row|<cell|CH2>|<cell|1.370000e-11>|<cell|0.00000e+00>|<cell|2.300000e-09>|<cell|0.00000e+00>|<cell|3.300000e-07>|<cell|4.000000e-07>|<cell|315.060>>|<row|<cell|G40_mb>|<cell|1.170000e-06>|<cell|1.200000e-17>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.800000e-05>|<cell|6.000000e-07>|<cell|7.025>>|<row|<cell|G48_mb>|<cell|1.310000e-04>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.100000e-12>|<cell|1.900000e-03>|<cell|6.800000e-07>|<cell|8.489>>|<row|<cell|G48mc>|<cell|7.080000e-11>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|2.400000e-06>|<cell|2.400000e-06>|<cell|2.681>>|<row|<cell|G55mc>|<cell|1.390000e-10>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|7.300000e-07>|<cell|7.300000e-07>|<cell|179.720>>|<row|<cell|G59mc>|<cell|1.220000e-10>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|4.900000e-07>|<cell|4.900000e-07>|<cell|264.597>>|<row|<cell|G60_mb>|<cell|1.030000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|2.000000e-12>|<cell|2.400000e-03>|<cell|2.400000e-03>|<cell|213.472>>|<row|<cell|G60mc>|<cell|1.030000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|2.000000e-12>|<cell|2.400000e-03>|<cell|2.400000e-03>|<cell|212.088>>|<row|<cell|H3O>|<cell|1.090000e-09>|<cell|0.00000e+00>|<cell|3.700000e-09>|<cell|0.00000e+00>|<cell|3.100000e-07>|<cell|3.600000e-07>|<cell|1344.764>>|<row|<cell|NH2>|<cell|1.820000e-09>|<cell|0.00000e+00>|<cell|2.300000e-09>|<cell|0.00000e+00>|<cell|3.800000e-07>|<cell|4.300000e-07>|<cell|290.156>>|<row|<cell|NH3>|<cell|1.360000e-09>|<cell|0.00000e+00>|<cell|2.700000e-09>|<cell|0.00000e+00>|<cell|3.200000e-07>|<cell|3.500000e-07>|<cell|1367.722>>|<row|<cell|NH4>|<cell|4.670000e-10>|<cell|0.00000e+00>|<cell|3.400000e-09>|<cell|0.00000e+00>|<cell|2.300000e-07>|<cell|3.900000e-07>|<cell|5202.509>>|<row|<cell|biggs>|<cell|3.620000e-09>|<cell|0.00000e+00>|<cell|1.900000e-12>|<cell|9.100000e-13>|<cell|3.500000e-08>|<cell|4.200000e-08>|<cell|14.316>>|<row|<cell|broyden25>|<cell|2.320000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|7.000000e-09>|<cell|7.300000e-09>|<cell|1774.592>>|<row|<cell|buck4>|<cell|1.390000e-11>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|4.400000e-07>|<cell|2.400000e-07>|<cell|21.185>>|<row|<cell|buck5>|<cell|4.920000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|5.400000e-04>|<cell|2.800000e-04>|<cell|194.738>>|<row|<cell|cancer_100>|<cell|5.250000e-13>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|3.400000e-15>|<cell|1.700000e-08>|<cell|3.400000e-08>|<cell|396.231>>|<row|<cell|checker_1.5>|<cell|2.830000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.300000e-06>|<cell|1.200000e-06>|<cell|41.693>>|<row|<cell|chs_5000>|<cell|1.210000e-17>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.500000e-07>|<cell|1.500000e-07>|<cell|36.825>>|<row|<cell|cnhil10>|<cell|3.980000e-08>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|8.500000e-09>|<cell|2.000000e-08>|<cell|63.071>>|<row|<cell|cphil12>|<cell|8.550000e-10>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|9.300000e-09>|<cell|9.300000e-09>|<cell|259.832>>|<row|<cell|diamond_patch>|<cell|5.000000e-01>|<cell|-0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|-9.400000e-01>|<cell|0.00000e+00>|<cell|Failed>>|<row|<cell|e_moment_quad>|<cell|2.370000e-10>|<cell|0.00000e+00>|<cell|5.400000e-09>|<cell|0.00000e+00>|<cell|-1.900000e-06>|<cell|1.100000e-08>|<cell|334.162>>|<row|<cell|e_moment_stable>|<cell|2.180000e-09>|<cell|7.700000e-19>|<cell|3.500000e-07>|<cell|8.500000e-10>|<cell|-1.100000e-05>|<cell|4.900000e-08>|<cell|256.368>>|<row|<cell|foot>|<cell|1.130000e-04>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|6.700000e-11>|<cell|-2.300000e-03>|<cell|1.900000e-04>|<cell|12.878>>|<row|<cell|hamming_8_3_4>|<cell|2.740000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.600000e-08>|<cell|1.600000e-08>|<cell|38.635>>|<row|<cell|hamming_9_5_6>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|Failed>>|<row|<cell|hand>|<cell|4.500000e-08>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|5.200000e-14>|<cell|1.100000e-06>|<cell|5.200000e-07>|<cell|2.431>>|<row|<cell|ice_2.0>|<cell|5.330000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.400000e-06>|<cell|1.200000e-06>|<cell|372.376>>|<row|<cell|inc_1200>|<cell|9.920000e-06>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|-2.300000e-04>|<cell|2.500000e-07>|<cell|128.205>>|<row|<cell|mater-5>|<cell|4.830000e-11>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|8.400000e-07>|<cell|8.400000e-07>|<cell|23.988>>|<row|<cell|mater-6>|<cell|1.670000e-10>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.100000e-06>|<cell|1.100000e-06>|<cell|61.504>>|<row|<cell|neosfbr25>|<cell|3.640000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|2.100000e-07>|<cell|2.100000e-07>|<cell|1084.074>>|<row|<cell|neosfbr30e8>|<cell|2.650000e-08>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|4.000000e-07>|<cell|4.200000e-07>|<cell|5924.588>>|<row|<cell|neu1>|<cell|3.650000e-06>|<cell|9.200000e-18>|<cell|4.800000e-07>|<cell|5.700000e-10>|<cell|-1.300000e-05>|<cell|5.000000e-06>|<cell|109.041>>|<row|<cell|neu1g>|<cell|7.850000e-10>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|4.000000e-09>|<cell|6.400000e-08>|<cell|89.921>>|<row|<cell|neu2>|<cell|2.480000e-08>|<cell|2.700000e-17>|<cell|3.900000e-07>|<cell|8.300000e-10>|<cell|-1.100000e-07>|<cell|4.600000e-08>|<cell|108.879>>|<row|<cell|neu2c>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|Failed>>|<row|<cell|neu2g>|<cell|4.280000e-10>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|9.900000e-09>|<cell|4.100000e-08>|<cell|90.972>>|<row|<cell|neu3>|<cell|8.110000e-09>|<cell|0.00000e+00>|<cell|2.600000e-07>|<cell|0.00000e+00>|<cell|-8.000000e-07>|<cell|3.000000e-08>|<cell|1731.641>>|<row|<cell|neu3g>|<cell|9.090000e-09>|<cell|1.600000e-17>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|6.500000e-09>|<cell|2.500000e-08>|<cell|1807.446>>|<row|<cell|p_auss2_3.0>|<cell|6.030000e-08>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|-2.200000e-06>|<cell|3.600000e-07>|<cell|451.407>>|<row|<cell|prob_2_4_0>|<cell|2.380000e-15>|<cell|0.00000e+00>|<cell|4.400000e-10>|<cell|0.00000e+00>|<cell|-7.100000e-07>|<cell|6.200000e-07>|<cell|221.371>>|<row|<cell|prob_2_4_1>|<cell|9.810000e-15>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|4.100000e-07>|<cell|4.100000e-07>|<cell|98.066>>|<row|<cell|rabmo>|<cell|9.160000e-10>|<cell|0.00000e+00>|<cell|1.600000e-07>|<cell|0.00000e+00>|<cell|-4.300000e-05>|<cell|1.500000e-08>|<cell|182.872>>|<row|<cell|reimer5>|<cell|1.740000e-09>|<cell|0.00000e+00>|<cell|1.100000e-07>|<cell|0.00000e+00>|<cell|-3.600000e-05>|<cell|4.600000e-09>|<cell|1862.190>>|<row|<cell|rendl>|<cell|1.530000e-10>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|4.100000e-07>|<cell|4.100000e-07>|<cell|7.351>>|<row|<cell|ros_2000>|<cell|5.510000e-16>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|1.400000e-07>|<cell|1.400000e-07>|<cell|3.760>>|<row|<cell|rose15>|<cell|8.550000e-08>|<cell|7.900000e-18>|<cell|1.300000e-07>|<cell|3.000000e-10>|<cell|-7.200000e-05>|<cell|1.500000e-08>|<cell|104.258>>|<row|<cell|sensor_1000b>|<cell|6.940000e-10>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|6.400000e-07>|<cell|1.400000e-07>|<cell|203.969>>|<row|<cell|shmup4>|<cell|1.770000e-08>|<cell|0.00000e+00>|<cell|1.400000e-13>|<cell|0.00000e+00>|<cell|9.700000e-07>|<cell|5.000000e-07>|<cell|63.314>>|<row|<cell|shmup5>|<cell|1.070000e-05>|<cell|0.00000e+00>|<cell|2.400000e-10>|<cell|0.00000e+00>|<cell|6.300000e-05>|<cell|5.100000e-07>|<cell|770.997>>|<row|<cell|spar060-020-1_LS>|<cell|7.770000e-08>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|9.100000e-07>|<cell|9.900000e-09>|<cell|779.882>>|<row|<cell|swissroll>|<cell|9.480000e-06>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|-7.300000e-03>|<cell|2.100000e-07>|<cell|37.456>>|<row|<cell|taha1a>|<cell|1.800000e-10>|<cell|0.00000e+00>|<cell|2.400000e-07>|<cell|0.00000e+00>|<cell|-1.900000e-07>|<cell|2.400000e-07>|<cell|196.912>>|<row|<cell|taha1b>|<cell|8.660000e-09>|<cell|0.00000e+00>|<cell|2.100000e-10>|<cell|0.00000e+00>|<cell|7.200000e-08>|<cell|9.100000e-08>|<cell|670.662>>|<row|<cell|taha1c>|<cell|9.900000e-09>|<cell|0.00000e+00>|<cell|4.100000e-07>|<cell|0.00000e+00>|<cell|1.800000e-05>|<cell|9.900000e-04>|<cell|2199.454>>|<row|<cell|theta12>|<cell|4.920000e-08>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|2.500000e-07>|<cell|2.900000e-07>|<cell|1095.150>>|<row|<cell|theta102>|<cell|7.370000e-11>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|6.600000e-07>|<cell|6.600000e-07>|<cell|5378.786>>|<row|<cell|theta123>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|1.00000e+00>|<cell|Failed>>|<row|<cell|tiger_texture>|<cell|9.230000e-04>|<cell|0.00000e+00>|<cell|1.200000e-12>|<cell|8.700000e-11>|<cell|1.600000e-03>|<cell|2.000000e-03>|<cell|53.752>>|<row|<cell|torusg3-15>|<cell|1.860000e-11>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|3.300000e-07>|<cell|3.300000e-07>|<cell|24.455>>|<row|<cell|trto4>|<cell|8.580000e-08>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|-1.400000e-06>|<cell|1.800000e-07>|<cell|7.633>>|<row|<cell|trto5>|<cell|1.430000e-06>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|-4.700000e-05>|<cell|2.400000e-07>|<cell|82.491>>|<row|<cell|vibra4>|<cell|6.030000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|2.000000e-03>|<cell|1.300000e-03>|<cell|33.061>>|<row|<cell|vibra5>|<cell|3.690000e-06>|<cell|0.00000e+00>|<cell|1.600000e-06>|<cell|0.00000e+00>|<cell|5.600000e-05>|<cell|1.100000e-04>|<cell|326.107>>|<row|<cell|yalsdp>|<cell|5.290000e-09>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|0.00000e+00>|<cell|2.400000e-06>|<cell|2.400000e-06>|<cell|125.136>>>>>
        </with>
      </big-table|Mittelmann's Benchmark Test>
    </with>
  </tiny>
</body>

<\initial>
  <\collection>
    <associate|font-base-size|10>
    <associate|info-flag|minimal>
    <associate|page-medium|paper>
    <associate|page-screen-margin|true>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|1|7>>
    <associate|auto-11|<tuple|5.3|7>>
    <associate|auto-12|<tuple|1|8>>
    <associate|auto-13|<tuple|6|9>>
    <associate|auto-14|<tuple|6.1|9>>
    <associate|auto-15|<tuple|2|9>>
    <associate|auto-16|<tuple|6.2|9>>
    <associate|auto-17|<tuple|3|10>>
    <associate|auto-18|<tuple|6.3|10>>
    <associate|auto-19|<tuple|4|10>>
    <associate|auto-2|<tuple|2|2>>
    <associate|auto-20|<tuple|6.4|10>>
    <associate|auto-21|<tuple|5|10>>
    <associate|auto-22|<tuple|6.5|11>>
    <associate|auto-23|<tuple|6.6|11>>
    <associate|auto-24|<tuple|7|11>>
    <associate|auto-25|<tuple|8|12>>
    <associate|auto-26|<tuple|8|12>>
    <associate|auto-27|<tuple|A|13>>
    <associate|auto-28|<tuple|6|13>>
    <associate|auto-3|<tuple|3|2>>
    <associate|auto-4|<tuple|3.1|2>>
    <associate|auto-5|<tuple|3.2|3>>
    <associate|auto-6|<tuple|4|5>>
    <associate|auto-7|<tuple|5|6>>
    <associate|auto-8|<tuple|5.1|6>>
    <associate|auto-9|<tuple|5.2|6>>
    <associate|bib-aps2019mosek|<tuple|1|12>>
    <associate|bib-benson1999mixed|<tuple|4|12>>
    <associate|bib-benson2000solving|<tuple|3|12>>
    <associate|bib-benson2008algorithm|<tuple|2|12>>
    <associate|bib-biswas2004semidefinite|<tuple|5|12>>
    <associate|bib-borchers1999sdplib|<tuple|6|12>>
    <associate|bib-borchers2006csdp|<tuple|7|12>>
    <associate|bib-copt|<tuple|8|12>>
    <associate|bib-fujisawa1997exploiting|<tuple|9|12>>
    <associate|bib-goemans1995improved|<tuple|10|12>>
    <associate|bib-hayashi2016quantum|<tuple|11|12>>
    <associate|bib-kocvara2006pensdp|<tuple|12|12>>
    <associate|bib-kwasniewiczimplementation|<tuple|13|12>>
    <associate|bib-laurent2005semidefinite|<tuple|15|12>>
    <associate|bib-laurent2009sums|<tuple|14|12>>
    <associate|bib-majumdar2020recent|<tuple|16|12>>
    <associate|bib-mittelmann2003independent|<tuple|17|12>>
    <associate|bib-polik2007sedumi|<tuple|18|12>>
    <associate|bib-potra1998homogeneous|<tuple|20|12>>
    <associate|bib-potra1998superlinearly|<tuple|19|12>>
    <associate|bib-qu2020diagonal|<tuple|21|12>>
    <associate|bib-so2007theory|<tuple|22|12>>
    <associate|bib-toh2002note|<tuple|23|12>>
    <associate|bib-toh2012implementation|<tuple|24|12>>
    <associate|bib-vandenberghe1996semidefinite|<tuple|25|12>>
    <associate|bib-wolkowicz2005semidefinite|<tuple|26|12>>
    <associate|bib-xu1996simplified|<tuple|27|12>>
    <associate|bib-yamashita2012latest|<tuple|28|13>>
    <associate|bib-yang2015sdpnal|<tuple|29|13>>
    <associate|sec2|<tuple|3|2>>
    <associate|sec3|<tuple|5|6>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      benson2008algorithm

      goemans1995improved

      laurent2005semidefinite

      vandenberghe1996semidefinite

      laurent2009sums

      hayashi2016quantum

      biswas2004semidefinite

      so2007theory

      wolkowicz2005semidefinite

      benson2008algorithm

      copt

      aps2019mosek

      polik2007sedumi

      toh2012implementation

      borchers2006csdp

      yamashita2012latest

      kocvara2006pensdp

      kwasniewiczimplementation

      yang2015sdpnal

      majumdar2020recent

      potra1998superlinearly

      potra1998homogeneous

      benson1999mixed

      benson2000solving

      benson2008algorithm

      benson1999mixed

      benson2000solving

      benson2008algorithm

      xu1996simplified

      potra1998homogeneous

      fujisawa1997exploiting

      toh2002note

      benson2008algorithm

      mittelmann2003independent

      borchers1999sdplib

      qu2020diagonal

      borchers1999sdplib

      mittelmann2003independent

      mittelmann2003independent

      benson2008algorithm
    </associate>
    <\associate|figure>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Pipeline of
      <with|font-family|<quote|tt>|language|<quote|verbatim>|HDSDP>>|<pageref|auto-10>>
    </associate>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Iteration monitor
      of Phase <with|font-family|<quote|tt>|language|<quote|verbatim>|A>>|<pageref|auto-12>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|2>||Max-cut
      problems>|<pageref|auto-15>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|3>||Graph partitioning
      problems>|<pageref|auto-17>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|4>||Optimal diagonal
      pre-conditioning problems>|<pageref|auto-19>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|5>||Feature of several
      benchmark problems>|<pageref|auto-21>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|6>||Mittelmann's
      Benchmark Test>|<pageref|auto-28>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Formulation
      and Notations> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Homogeneous
      Dual Scaling Algorithm > <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Dual-scaling Algorithm
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Homogeneous Dual-scaling
      Algorithm <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Initalization
      and Feasibility Certificate> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>HDSDP
      Software> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>

      <with|par-left|<quote|1tab>|5.1<space|2spc>Pre-solver
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|1tab>|5.2<space|2spc>Two-phase Optimization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|1tab>|5.3<space|2spc>Iteration Monitor
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Computational
      results> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.5fn>

      <with|par-left|<quote|1tab>|6.1<space|2spc>Maximum-cut
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|1tab>|6.2<space|2spc>Graph Partitioning
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|1tab>|6.3<space|2spc>Optimal Diagonal
      Pre-conditioning <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18>>

      <with|par-left|<quote|1tab>|6.4<space|2spc>Other Problems
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <with|par-left|<quote|1tab>|6.5<space|2spc>Mittelmann's Benchmark Test
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22>>

      <with|par-left|<quote|1tab>|6.6<space|2spc>When to use DSDP/HDSDP
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|7<space|2spc>Conclusions>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|8<space|2spc>Acknowledgement>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-25><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|References>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Appendix
      A<space|2spc>Mittlelmann's Benchmark Test>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-27><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>