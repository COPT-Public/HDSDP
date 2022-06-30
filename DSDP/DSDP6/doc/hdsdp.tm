<TeXmacs|2.1>

<style|<tuple|generic|centered-program>>

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
    \;

    \;

    <date|>
  </author-affiliation>>>>

  <abstract-data|<\abstract>
    HDSDP is a numerical software designed for solving the semidefinite
    programming problems (SDP). The main framework of HDSDP resembles
    DSDP<cite|benson2008algorithm> and several new features, especially a
    dual-scaling algorithm based on simplified homogeneous self-dual
    embedding, have been implemented. The self-dual embedding has several
    desirable features which enhance the numerical stability and several new
    heuristics and computational techniques are designed to accelerate its
    convergence. HDSDP aims to show how dual-scaling algorithms benefit from
    self-dual embedding and it is developed parallel to <verbatim|DSDP5.8>.
    Numerical experiments over several SDP benchmark datasets prove the
    robustness and efficiency of HDSDP.
  </abstract>>

  <section|Introduction>

  Semi-definite programming (SDP) is a mathematical programming problem
  defined by

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<X>>>|<cell|<around*|\<langle\>|<C>,<X>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<fa><X>=<b>>|<cell|>>|<row|<cell|>|<cell|<X>\<in\>\<bbb-S\><rsub|+><rsup|n>,>|<cell|>>>>
  </eqnarray*>

  where we do linear optimization subject to affine constraints over the cone
  of positive-semidefinite matrices. Due to its wide modeling ability, SDP
  has been employed as a powerful model in various communities including
  combinatorial optimization <cite|laurent2005semidefinite>, numerical
  analysis <cite|vandenberghe1999applications> and dynamic systems
  <cite|vandenberghe1996semidefinite>. Moreover, since the semial work of
  <cite|goemans1995improved>, SDP has also proved an efficient backend in
  modern discrete and nonconvex programming solvers to provide tight bounds
  using the technique of SDP relaxtion <cite|laurent2005semidefinite>.\ 

  While SDP is quite an useful tool in many applications, a fundamental issue
  is how to numerically solve them. Theoretically speaking, SDP is a convex
  conic problem which admits efficient polynomial-time algorithms and for
  general SDPs, the interior point method (IPM) is known as a most robust and
  efficient algorithm. Since 1990s, several high-performance SDPs softwares
  have been developed, including <verbatim|DSDP> <cite|benson2008algorithm>
  <verbatim|Mosek> <cite|aps2019mosek>, <verbatim|Sedumi>
  <cite|sturm1999using>, <verbatim|SDPT3> <cite|toh1999sdpt3>,
  <verbatim|CSDP> <cite|borchers1999csdp> and <verbatim|SDPA>
  <cite|yamashita2003implementation>. These numerical solvers are based on
  different variants of the IPM and enjoy nice convergence behavior both
  theoretically and practically. Among all the IPM-based SDP solvers, all of
  them implement the path-following primal-dual approach either using
  infeasible start <cite|potra1998superlinearly> or the homogeneous self-dual
  embedding <cite|potra1998homogeneous><cite|aps2019mosek> but with
  <verbatim|DSDP> being an exception. <verbatim|DSDP>, as its name suggests,
  implements a purely dual IPM algorithm and is based on the potential
  reduction framework proposed in <cite|benson2001dsdp3>. For semi-definite
  programs, dual method enjoys sparsity and low-rank structure of data
  <cite|benson2001dsdp3> at the cost of weaker ability to detect
  infeasibility. One natural idea is:

  <center|<em|Can we let dual algorithm detect infeasibility better as in
  primal-dual algorithms?>>

  In this work, we follow the track of <verbatim|DSDP> and make several
  practical improvements to improve the original solver both in infeasibility
  detection and solution efficiency and more detailedly,

  <\itemize>
    <item>The intuition of homogeneous self-dual model is incorporated into
    the dual-scaling solver

    <item>Several new computationally tricks are added to improve its
    convergence
  </itemize>

  and a new piece of SDP software, named <verbatim|HDSDP> (Homogeneous
  Dual-scaling SDP) is presented. In parallel to <verbatim|DSDP>, the
  software is written in ANSI C standard and can be called as a sub-routine
  library.

  The rest of the paper is organized as follows. In <strong|Section
  <reference|sec2>>, we introduce the algorithm implemented in
  <verbatim|HDSDP>. In <strong|Section <reference|sec3>>, we give some
  details of the newly-added computational tricks, some of which may be
  presented beyond the scope just for SDPs. Last we give some preliminary
  numerical results of the <verbatim|HDSDP> solver.

  <\remark>
    We need to clarify that this paper is currently formatted more like a
    technical document for the solver. When the solver is finally complete,
    there would be more detailed descriptions of the tricks and algorithms.
  </remark>

  <section|Homogeneous Dual Scaling Algorithm ><label|sec2>

  <subsection|<strong|Notations>>

  Without loss of generality, we use capitalized letter <math|<A>> to denote
  matrix and <math|<math-bf|><a>> to denote vector; <math|<a><rsup|j>> and
  <math|<a><rsub|i>> denote corresponding the column and row of the matrix
  respectively. Specially we define <math|<fa><X>:\<bbb-S\><rsup|n>\<rightarrow\>\<bbb-R\><rsup|m>=<matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<A>,<X><rsub|1>|\<rangle\>>>>|<row|<cell|\<vdots\>>>|<row|<cell|<around*|\<langle\>|<A><rsub|m>,<X><rsub|m>|\<rangle\>>>>>>>>
  and <math|<fa><rsup|\<ast\>><y>\<assign\><big|sum><rsub|i=1><rsup|m><A><rsub|i>y<rsub|i>>
  denotes its adjoint. We denote <math|<around*|\<langle\>|<A>,<math-bf|B>|\<rangle\>>\<assign\><big|sum><rsub|i
  j>a<rsub|i j>b<rsub|i j>> to be the matrix inner product and
  <math|<around*|\<\|\|\>|<A>|\<\|\|\>><rsub|F>\<assign\><sqrt|<big|sum><rsub|i
  j>a<rsub|i j><rsup|2>>> denotes matrix Frobenius norm.
  <math|\<bbb-S\><rsub|+><rsup|n>> denotes the cone of positive semidefinite
  matrices.

  <subsection|Standard Form>

  HDSDP solves the following standard semi-definite programming problem

  <strong|Dual>

  <\equation*>
    <tabular|<tformat|<table|<row|<cell|<around*|(|P|)>>|<cell|min<rsub|<X>>>|<cell|<around*|\<langle\>|<C>,<X>|\<rangle\>>>|<cell|>>|<row|<cell|>|<cell|<text|subject
    to>>|<cell|<fa><X>=<b>>|<cell|>>|<row|<cell|>|<cell|>|<cell|<X>\<in\>\<bbb-S\><rsub|+><rsup|n>>|<cell|>>>>><space|2em><tabular|<tformat|<table|<row|<cell|<around*|(|D|)>>|<cell|max<rsub|<y>,<bs>>>|<cell|<b><rsup|\<top\>><y>>>|<row|<cell|>|<cell|<text|subject
    to>>|<cell|<fa><rsup|\<ast\>><y>+<bs>=<C>>>|<row|<cell|>|<cell|>|<cell|<bs>\<in\>\<bbb-S\><rsub|+><rsup|n>>>>>>
  </equation*>

  And we assume that the primal-dual pair is well-defined. i.e., if both
  primal and dual problems are feasible, then a zero duality gap is assumed.

  <subsection|HSD and Dual-scaling Algorithm>

  Given standard for SDP, HSD (homogeneous self-dual embedding) is a skew
  symmetric system \ whose non-trivial interior point solution certificates
  primal-dual feasibility. Here we consider the (simplified) HSD
  <cite|xu1996simplified> for standard form SDP

  <\eqnarray*>
    <tformat|<table|<row|<cell|<fa><X>-<b>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|-<fa><rsup|\<ast\>><y>+<C>\<tau\>-<bs>>|<cell|=>|<cell|<0>>>|<row|<cell|<b><rsup|\<top\>><y>-<around*|\<langle\>|<C>,<X>|\<rangle\>>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|<X><bs>=\<mu\><I>,>|<cell|>|<cell|\<kappa\>\<tau\>=\<mu\>>>|<row|<cell|<X>,<bs>\<in\>\<bbb-S\><rsub|+><rsup|n>.>|<cell|>|<cell|\<kappa\>,\<tau\>\<geq\>0,>>>>
  </eqnarray*>

  where <math|\<kappa\>,\<tau\>\<geq\>0> are two homogenizing variables. The
  central path <math|\<cal-C\><around*|(|\<mu\>|)>> of HSD is characterizes
  by <math|<around*|{|<around*|(|<X><around*|(|\<mu\>|)>,<bs><around*|(|\<mu\>|)>,<y><around*|(|\<mu\>|)>,\<kappa\><around*|(|\<mu\>|)>,\<tau\><around*|(|\<mu\>|)>|)>|}>>
  such that

  <\eqnarray*>
    <tformat|<table|<row|<cell|<fa><X>-<b>\<tau\>>|<cell|=>|<cell|<0>>>|<row|<cell|-<fa><rsup|\<ast\>><y>+<C>\<tau\>-<bs>>|<cell|=>|<cell|<0>>>|<row|<cell|<b><rsup|\<top\>><y>-<around*|\<langle\>|<C>,<X>|\<rangle\>>-\<kappa\>>|<cell|=>|<cell|0>>|<row|<cell|<X>=\<mu\><bs><rsup|-1>,>|<cell|>|<cell|\<kappa\>=\<mu\>\<tau\><rsup|-1>>>|<row|<cell|<X>,<bs>\<in\>\<bbb-S\><rsub|+><rsup|n>,>|<cell|>|<cell|\<kappa\>,\<tau\>\<geq\>0>>>>
  </eqnarray*>

  HDSDP basically follows the simplified homogeneous self-dual embedding
  incorporating infeasible-start interior point method. Given
  <math|<around*|(|<X>\<in\>\<bbb-S\><rsub|+><rsup|n>,<bs>\<in\>\<bbb-S\><rsub|+><rsup|n>,<y>,\<kappa\>\<geq\>0,\<tau\>\<geq\>0|)>>,
  and defining HSD residuals by,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<r><rsub|<p>>>|<cell|\<assign\>>|<cell|\<cal-A\><X>-<b>\<tau\>>>|<row|<cell|<R><rsub|<y>>>|<cell|\<assign\>>|<cell|-<fa><rsup|\<ast\>><y>+<C>\<tau\>-<bs>>>|<row|<cell|r<rsub|\<kappa\>>>|<cell|\<assign\>>|<cell|<b><rsup|\<top\>><y>-<around*|\<langle\>|<C>,<X>|\<rangle\>>-\<kappa\>>>|<row|<cell|<R><rsub|\<mu\>>>|<cell|\<assign\>>|<cell|<X>-\<mu\><bs><rsup|-1>>>|<row|<cell|r<rsub|\<mu\>>>|<cell|\<assign\>>|<cell|\<kappa\>-\<mu\>\<tau\><rsup|-1>,>>>>
  </eqnarray*>

  <verbatim|HDSDP> takes Newton's step towards feasibility of the whole HSD
  and solves

  <strong|HSD-Newton>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\>\<Delta\><X>-<b>\<Delta\>\<tau\>>|<cell|=>|<cell|-<r><rsub|<p>>>>|<row|<cell|-\<cal-A\><rsup|\<ast\>>\<Delta\><y>+\<Delta\>\<tau\><C>-\<Delta\><bs>>|<cell|=>|<cell|-\<gamma\><R><rsub|<y>>>>|<row|<cell|<b><rsup|\<top\>>\<Delta\><y>-<around*|\<langle\>|<C>,\<Delta\><X>|\<rangle\>>-\<Delta\>\<kappa\>>|<cell|=>|<cell|-r<rsub|\<kappa\>>>>|<row|<cell|\<Delta\><X>+\<mu\><bs><rsup|-1>\<Delta\><bs><bs><rsup|-1>>|<cell|=>|<cell|-<R><rsub|\<mu\><rsub|1>>>>|<row|<cell|\<kappa\>+\<mu\>\<tau\><rsup|-2>\<Delta\>\<tau\>>|<cell|=>|<cell|-r<rsub|\<mu\><rsub|2>>.>>>>
  </eqnarray*>

  to obtain the Newton direction. In practice the system is computed via
  Schur complement matrix <math|<M>>

  <\equation*>
    <M>\<assign\><matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<A><rsub|1>,<bs><rsup|-1><A><rsub|1><bs><rsup|-1>|\<rangle\>>>|<cell|\<cdots\>>|<cell|<around*|\<langle\>|<A><rsub|1>,<bs><rsup|-1><A><rsub|m><bs><rsup|-1>|\<rangle\>>>>|<row|<cell|\<vdots\>>|<cell|\<ddots\>>|<cell|\<vdots\>>>|<row|<cell|<around*|\<langle\>|<A><rsub|m>,<bs><rsup|-1><A><rsub|1><bs><rsup|-1>|\<rangle\>>>|<cell|\<cdots\>>|<cell|<around*|\<langle\>|<A><rsub|m>,<bs><rsup|-1><A><rsub|m><bs><rsup|-1>|\<rangle\>>>>>>>\<in\>\<bbb-R\><rsup|m\<times\>m>
  </equation*>

  and we record it for later use. After obtaining a Newton's direction,
  <verbatim|HDSDP> chooses a proper stepsize
  <math|\<alpha\>\<in\><around*|(|0,1|]>> such that
  <math|<bs>+\<alpha\>\<Delta\><bs>\<in\>\<bbb-S\><rsub|++><rsup|n>,\<tau\>+\<alpha\>\<Delta\>\<tau\>\<geq\>0>
  and updates

  <\eqnarray*>
    <tformat|<table|<row|<cell|<y>>|<cell|\<leftarrow\>>|<cell|<y>+\<alpha\>\<Delta\><y>>>|<row|<cell|\<tau\>>|<cell|\<leftarrow\>>|<cell|\<tau\>+\<alpha\>\<Delta\>\<tau\>>>|<row|<cell|<bs>>|<cell|\<leftarrow\>>|<cell|<bs>+\<alpha\>\<Delta\><bs>.>>>>
  </eqnarray*>

  Here we note that <math|<X>> and <math|\<kappa\>> do not appear in the
  updates since they are regarded as the primal variables. After taking the
  Newton step, the dual infeasibility is decreased by <math|1-\<alpha\>>. The
  above procedure is repeated till the dual infeasibility is eliminated or
  the HSD certificates dual infeasibility.

  <subsection|Barrier and Centrality Charactrization>

  One main difference between the dual algorithm from primal-dual methods is
  that dual algorithm has no reference of the the primal objective value. At
  a cost, we need to instead evaluate how far the iteration is deviating from
  the central path to determine the strategy of iteration. For a given
  barrier parameter <math|\<mu\>\<gtr\>0>, the proximity measure
  <math|\<delta\><around*|(|<bs>\<comma\>\<tau\>,\<mu\>|)>> is obtained from
  the following orthogonal projection

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<X>>>|<cell|<around*|\<\|\|\>|<frac|<bs><rsup|1/2><X><bs><rsup|1/2>|\<mu\>>-<I>|\<\|\|\>><rsub|F><rsup|2>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<fa><X>=<b>\<tau\>.>|<cell|>>>>
  </eqnarray*>

  Intuitively, the problem projects the dual variable onto the central path
  characterized by <math|\<mu\>> buts ignores the conic constraint. The
  solution to the projection, by optimality condition, is given by

  <\equation*>
    <X><around*|(|\<mu\>|)>\<assign\>\<mu\><around*|(|<bs><rsup|-1>+<bs><rsup|-1><around*|(|<fa><rsup|\<ast\>>\<Delta\><rsub|\<mu\>><y>|)><bs><rsup|-1>|)>=\<mu\><bs><rsup|-1><around*|(|<C>\<tau\>-<fa><rsup|\<ast\>><around*|(|<y>-\<Delta\><rsub|\<mu\>><y>|)>-<R><rsub|<y>>|)><bs><rsup|-1>
  </equation*>

  where <math|\<Delta\><rsub|\<mu\>><y>=<frac|\<tau\>|\<mu\>><M><rsup|-1><b>-<M><rsup|-1><fa><bs><rsup|-1>>
  and one desirable feature of the dual algorithm is that if
  <math|<C>\<tau\>-<fa><rsup|\<ast\>><around*|(|<y>-\<Delta\><rsub|\<mu\>><y>|)>-<R><rsub|<y>>>
  is positive definite, we can implicitly derive the corresponding primal
  objective value without actual reference to the primal variable.

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|\<langle\>|<C>,<X><around*|(|\<mu\>|)>|\<rangle\>>>|<cell|=>|<cell|<frac|1|\<tau\>><around*|\<langle\>|<R><rsub|<y>>+<fa><rsup|\<ast\>><y>+<bs>,<X><around*|(|\<mu\>|)>|\<rangle\>>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<tau\>><around*|\<langle\>|<R><rsub|<y>>+<bs>,<X><around*|(|\<mu\>|)>|\<rangle\>>+<b><rsup|\<top\>><y>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<tau\>><around*|\<langle\>|<R><rsub|<y>>+<bs>,\<mu\><around*|(|<bs><rsup|-1>+<bs><rsup|-1><around*|(|<fa><rsup|\<ast\>>\<Delta\><rsub|\<mu\>><y>|)><bs><rsup|-1>|)>|\<rangle\>>+<b><rsup|\<top\>><y>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<tau\>><around*|\<langle\>|<R><rsub|<y>>,\<mu\><around*|(|<bs><rsup|-1>+<bs><rsup|-1><around*|(|<fa><rsup|\<ast\>>\<Delta\><rsub|\<mu\>><y>|)><bs><rsup|-1>|)>|\<rangle\>>+<around*|\<langle\>|<bs>,\<mu\><around*|(|<bs><rsup|-1>+<bs><rsup|-1><around*|(|<fa><rsup|\<ast\>>\<Delta\><rsub|\<mu\>><y>|)><bs><rsup|-1>|)>|\<rangle\>>+<b><rsup|\<top\>><y>>>|<row|<cell|>|<cell|=>|<cell|<frac|\<mu\>|\<tau\>><around*|{|<around*|\<langle\>|<R><rsub|<y>>,<bs><rsup|-1>|\<rangle\>>+<around*|(|<fa><bs><rsup|-1><R><rsub|<y>><bs><rsup|-1>+<fa><bs><rsup|-1>|)><rsup|\<top\>>\<Delta\><rsub|\<mu\>><y>+n|}>+<b><rsup|\<top\>><y>.>>>>
  </eqnarray*>

  After getting a new objective value, we update the barrier parameter
  according to the duality gap as in primal-dual methods. Here we summarize
  the homogeneous dual-scaling algorithm as follows.

  <\named-algorithm|1. HDSDP Dual Infeasibility Elimination<space|11em>>
    <strong|input> <math|<y>,<fa>,<C>,<b>,<bs>,\<kappa\>,\<tau\>>

    <space|2em><strong|for> <math|i=1,\<ldots\>,K>

    <space|4em><strong|update> <math|\<mu\>>

    <space|4em><strong|compute> <math|\<Delta\><y>,\<Delta\>\<tau\>,\<Delta\><bs>>

    <space|4em><strong|choose> <math|\<alpha\>>

    <\eqnarray*>
      <tformat|<table|<row|<cell|<bs>>|<cell|\<leftarrow\>>|<cell|<bs>+\<alpha\>\<Delta\><bs>>>|<row|<cell|\<tau\>>|<cell|\<leftarrow\>>|<cell|\<tau\>+\<alpha\>\<Delta\>\<tau\>>>|<row|<cell|<y>>|<cell|\<leftarrow\>>|<cell|<y>+\<alpha\>\<Delta\><y>>>>>
    </eqnarray*>

    <space|4em><strong|check dual feasibility>

    <space|2em><strong|end for>

    <strong|end>
  </named-algorithm>

  Although the above algorithm succinctly describes the implemented
  algorithm, the true implementation is generally much more complicated than
  described and we present some of the interesting tricks in the following
  section.

  <section|Implementation Details><label|sec3>

  In this section we describe some implementation details of the
  <verbatim|HDSDP> solver. In <verbatim|HDSDP> We divided the solution into
  three consecutive phases, namely

  <\itemize>
    <item>presolving

    scales the coefficient of the problem data and performs decomposition

    <item>dual infeasibility elimination

    use HSD to obtain a dual feasible solution

    <item>primal feasibility check

    use dual potential reduction to pursue optimality.
  </itemize>

  <verbatim|HDSDP> combines re-implements most of the tricks from
  <verbatim|DSDP> with some newly added features. In this section the new
  features will be marked <verbatim|*> and the last phase, as largely the
  same as <verbatim|DSDP>, is not presented.

  <big-figure|<image|procedure.png|800px|||>|The inner procedure of HDSDP>

  <subsection|Presolve>

  In the current version DSDP performs rank-one detection and coefficient
  scaling operations.

  <subsubsection|<strong|Rank-one detection>*>

  This operation mainly detects the rank-one properties of the coefficient
  matrices in each block. Given a matrix <math|<A>> as DSDP data, we check
  whether <math|<A>=\<pm\><a><a><rsup|\<top\>>> for <math|<a>> extracted from
  <math|<A>> assuming it is rank-one. This is done by computing
  <math|<around*|\<\|\|\>|<A>-<a><a><rsup|\<top\>>|\<\|\|\>><rsub|F>> and
  takes <math|\<cal-O\><around*|(|n<rsup|2>|)>/\<cal-O\><around*|(|nnz|)>>
  time for dense and sparse structures respectively.

  <subsubsection|Sparse eigen-decomposition>

  <verbatim|HDSDP> implements a state-of-the-art eigen-decomposition solver
  from <verbatim|DSDP>. Given a sparse matrix <math|<A>>, <verbatim|DSDP>
  first finds a permutation that collects all its non-zero values to a small
  dense block, calls <verbatim|Lapack> routines and finally uses the inverse
  permutation to restore the true eigen-decomposition

  <\equation*>
    <matrix|<tformat|<table|<row|<cell|\<ast\>>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|<space|1em>>|<cell|<space|1em>>|<cell|\<ddots\>>|<cell|>|<cell|<space|1em>>>|<row|<cell|>|<cell|\<ddots\>>|<cell|>|<cell|\<ddots\>>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<ddots\>>|<cell|<space|1em>>|<cell|\<ast\>>>|<row|<cell|>|<cell|>|<cell|<space|1em>>|<cell|\<ast\>>|<cell|>>>>>\<Rightarrow\><matrix|<tformat|<cwith|2|2|1|3|cell-bborder|0ln>|<cwith|3|3|1|3|cell-tborder|0ln>|<cwith|1|1|1|3|cell-tborder|1ln>|<cwith|3|3|1|3|cell-bborder|1ln>|<cwith|4|4|1|3|cell-tborder|1ln>|<cwith|1|3|1|1|cell-lborder|1ln>|<cwith|1|3|3|3|cell-rborder|1ln>|<cwith|1|3|4|4|cell-lborder|1ln>|<table|<row|<cell|\<ast\>>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|<space|1em>>|<cell|<space|1em>>|<cell|\<ast\>>|<cell|>|<cell|>>|<row|<cell|>|<cell|\<ast\>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|\<ddots\>>|<cell|>>|<row|<cell|>|<cell|>|<cell|<space|1em>>|<cell|>|<cell|\<ddots\>>>>>>.
  </equation*>

  As a special case, if the number of non-zeros in the matrix falls below
  some threshold, <verbatim|HDSDP> will automatically apply Given's rotation
  to compute the eigen-decomposition.

  <subsubsection|Coefficient scaling>

  The second step of the presolving procedure scales the data to improve the
  numerical stability.

  In the coefficient scaling procedure, we scale matrices from different
  blocks. i.e., we scale <math|<A><rsub|1k>,<A><rsub|2 k>> simultaneously by

  <\equation*>
    <A><rsup|\<pi\>><rsub|1k>=<frac|<A><rsub|1k>|<sqrt|max<around*|{|<around*|\<\|\|\>|<A><rsub|1k>|\<\|\|\>><rsub|F>,<around*|\<\|\|\>|<A><rsub|2k>|\<\|\|\>><rsub|F>,<around*|\||b<rsub|k>|\|>|}>\<cdummy\>min<around*|{|<around*|\<\|\|\>|<A><rsub|1k>|\<\|\|\>><rsub|F>,<around*|\<\|\|\>|<A><rsub|2k>|\<\|\|\>><rsub|F>,<around*|\||b<rsub|k>|\|>|}>>>=\<pi\><rsub|k><rsup|-1><A><rsub|1k>,\<forall\>k\<in\><around*|[|m|]>
  </equation*>

  <\equation*>
    b<rsup|\<pi\>><rsub|k>=\<pi\><rsub|k><rsup|-1>b<rsub|k>,\<forall\>k\<in\><around*|[|m|]>
  </equation*>

  <subsection|Dual infeasibility elimination>

  Dual infeasibility elimination in DSDP aims to get a dual feasible solution
  as a starting point of the dual potential reduction.\ 

  <subsubsection|Initialization>

  <verbatim|HDSDP> initializes <math|<R><rsub|<y>>=-\<alpha\><I>> and
  <math|<bs><rsup|0>=<C>\<tau\><rsup|0>-<R><rsub|<y>>=<C>\<tau\><rsup|0>+\<alpha\><I>>,
  where <math|\<alpha\>> is chosen such that
  <math|<C>\<tau\><rsup|0>+\<alpha\><I>\<succ\><0>>.

  <subsubsection|Damped Dual Infeasibility Elimination*>

  Recall that in <verbatim|HDSDP> we compute a Newton's step towards dual
  feasibility and centrality

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\><around*|(|<X>+\<Delta\><X>|)>>|<cell|=>|<cell|<0>>>|<row|<cell|-\<cal-A\><rsup|\<ast\>>\<Delta\><y>-\<Delta\><bs>>|<cell|=>|<cell|-<R><rsub|<d>>>>|<row|<cell|\<Delta\><X>+\<mu\><bs><rsup|-1>\<Delta\><bs><bs><rsup|-1>>|<cell|=>|<cell|\<mu\><bs><rsup|-1>-<X>>>>>
  </eqnarray*>

  and for brevity we temporarily fix <math|\<tau\>=\<kappa\>=1>. Since the
  third line of the system linearizes the complementarity condition, it is
  inevitable that the first order approximation results in an error and, if
  we always pursues eliminating infeasibility completely, it is quite likely
  that we hit the barrier of the cone and never progress. To overcome this
  drawback, we adopt an adaptive damping stratgey that adjusts how aggressive
  the algorithm should be. To illustrate how the method works, we introduce a
  damping factor <math|\<gamma\>\<in\><around*|[|0,1|]>> into the Newton step
  and get

  <strong|Damped Dual Infeasibility>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-A\><around*|(|<X>+\<Delta\><X>|)>>|<cell|=>|<cell|<0>>>|<row|<cell|-\<cal-A\><rsup|\<ast\>>\<Delta\><y>-\<Delta\><bs>>|<cell|=>|<cell|-\<gamma\><R><rsub|<d>>>>|<row|<cell|\<Delta\><X>+\<mu\><bs><rsup|-1>\<Delta\><bs><bs><rsup|-1>>|<cell|=>|<cell|\<mu\><bs><rsup|-1>-<X>.>>>>
  </eqnarray*>

  Then we can re-write the damped update as

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<Delta\><y>>|<cell|=>|<cell|\<gamma\><M><rsup|-1>\<cal-A\><bs><rsup|-1><R><rsub|<d>><bs><rsup|-1>-<M><rsup|-1>\<cal-A\><bs><rsup|-1>>>|<row|<cell|>|<cell|\<backassign\>>|<cell|\<gamma\>\<Delta\><y><rsup|r>-\<Delta\><y><rsup|c>>>|<row|<cell|<y>+\<alpha\>\<Delta\><y>>|<cell|=>|<cell|<y>-\<alpha\>\<Delta\><y><rsup|c>+\<alpha\>\<gamma\>\<Delta\><y><rsup|r>>>|<row|<cell|<bs>+\<alpha\>\<Delta\><bs>>|<cell|=>|<cell|<bs>+\<alpha\>\<gamma\><R><rsub|<d>>-\<alpha\>\<cal-A\><rsup|\<ast\>>\<Delta\><y>>>|<row|<cell|>|<cell|=>|<cell|<bs>+\<alpha\>\<cal-A\><rsup|\<ast\>>\<Delta\><y><rsup|c>+\<alpha\>\<gamma\><around*|(|<R><rsub|<d>>-\<cal-A\><rsup|\<ast\>>\<Delta\><y><rsup|r>|)>>>>>
  </eqnarray*>

  To decide the damping factor, we first line-search of the barrier function
  till <math|-log det <around*|(|<bs>+\<alpha\><rsub|c>\<Delta\><bs>|)>\<leq\>-log
  det <bs>>, which guarantees the centrality. Then another line-search
  determines <math|<bs>+\<alpha\><rsub|c>\<cal-A\><rsup|\<ast\>>\<Delta\><y><rsup|c>+\<alpha\><rsub|c>\<gamma\><around*|(|<R><rsub|<d>>-\<cal-A\><rsup|\<ast\>>\<Delta\><y><rsup|r>|)>\<in\>\<bbb-S\><rsub|++><rsup|n>>
  so that the dual infeasibility is eliminated without compromising too much
  centrality.

  <subsubsection|Schur complement tricks*>

  This section is specially devoted to the Schur complement in
  <verbatim|HDSDP>, whose cost is often dominating when the number of
  constraints is large. The Schur complement in <verbatim|HDSDP> is given by

  \;

  <\equation*>
    <M>\<assign\><matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<A><rsub|1>,<bs><rsup|-1><A><rsub|1><bs><rsup|-1>|\<rangle\>>>|<cell|\<cdots\>>|<cell|<around*|\<langle\>|<A><rsub|1>,<bs><rsup|-1><A><rsub|m><bs><rsup|-1>|\<rangle\>>>>|<row|<cell|\<vdots\>>|<cell|\<ddots\>>|<cell|\<vdots\>>>|<row|<cell|<around*|\<langle\>|<A><rsub|m>,<bs><rsup|-1><A><rsub|1><bs><rsup|-1>|\<rangle\>>>|<cell|\<cdots\>>|<cell|<around*|\<langle\>|<A><rsub|m>,<bs><rsup|-1><A><rsub|m><bs><rsup|-1>|\<rangle\>>>>>>>.
  </equation*>

  In real practice, since <math|<M>> is symmetric, we only need the upper
  (low) triangular part of <math|<M>>

  <\equation*>
    <matrix|<tformat|<table|<row|<cell|<around*|\<langle\>|<A><rsub|1>,<bs><rsup|-1><A><rsub|1><bs><rsup|-1>|\<rangle\>>>|<cell|\<cdots\>>|<cell|<around*|\<langle\>|<A><rsub|1>,<bs><rsup|-1><A><rsub|m><bs><rsup|-1>|\<rangle\>>>>|<row|<cell|\<longrightarrow\>>|<cell|\<longrightarrow\>>|<cell|\<longrightarrow\>>>|<row|<cell|>|<cell|\<ddots\>>|<cell|\<vdots\>>>|<row|<cell|>|<cell|>|<cell|<around*|\<langle\>|<A><rsub|m>,<bs><rsup|-1><A><rsub|m><bs><rsup|-1>|\<rangle\>>>>>>>
  </equation*>

  and it is set up row by row. To accomodate the additional computation
  resulting from HSD, we adopt five strategies, two exploiting low-rank
  structure (from <verbatim|DSDP>) and three exploiting sparsity (from
  <verbatim|SDPA>) to compute <math|<M>>. Technique <strong|M1> and
  <strong|M2> require eigen-decomposition of
  <math|<A><rsub|i>=<big|sum><rsub|p>\<lambda\><rsub|p><a><rsub|i,p><a><rsub|i,p><rsup|\<top\>>>;
  <strong|M3>, <strong|M4> and <strong|M5> requires the inverse of
  <math|<bs><rsup|-1>>. Given the number of nonzero elements in
  <math|<around*|{|<A><rsub|i>|}>> (say <math|f<rsub|1>,\<ldots\>,f<rsub|m>>),
  we can compute a permutation <math|<around*|{|\<sigma\><rsub|<around*|(|1|)>>,\<ldots\>,\<sigma\><rsub|<around*|(|m|)>>|}>>
  of <math|<around*|{|1,\<ldots\>,m|}>> using <cite|fujisawa1997exploiting>
  and then set up the Schur matrix using one of the following five techniques
  that minimizes the number of float point operations.

  <strong|Technique M1>

  <\enumerate>
    <item><strong|Setup> <math|><math|<math-bf|B><rsub|\<sigma\><around*|(|i|)>>=<big|sum><rsub|p=1><rsup|r<rsub|\<sigma\><around*|(|i|)>>>\<lambda\><rsub|p><around*|(|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,p>|)><around*|(|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,p>|)><rsup|\<top\>>>

    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>
    \<sigma\><around*|(|j|)>>=<around*|\<langle\>|<math-bf|B><rsub|\<sigma\><around*|(|i|)>>,<A><rsub|\<sigma\><around*|(|j|)>>|\<rangle\>>,\<forall\>j\<geq\>i>
  </enumerate>

  The first step in <strong|M1> technique takes
  <math|r<rsub|\<sigma\><around*|(|i|)>><around*|(|n<rsup|2>+2\<kappa\>n<rsup|2>|)>>
  and the second takes <math|\<kappa\><big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>.
  The overall cost is then <math|r<rsub|\<sigma\><around*|(|i|)>><around*|(|n<rsup|2>+2\<kappa\>n<rsup|2>|)>+\<kappa\><big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>.

  <strong|Technique M2>

  <\enumerate>
    <item><strong|Setup> <math|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,p>,p=1,\<ldots\>,r<rsub|\<sigma\><around*|(|i|)>>>

    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>\<sigma\><around*|(|j|)>>=<big|sum><rsub|p=1><rsup|r<rsub|\<sigma\><around*|(|i|)>>><around*|(|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,p>|)><rsup|\<top\>><A><rsub|\<sigma\><around*|(|j|)>><around*|(|<bs><rsup|-1><a><rsub|\<sigma\><around*|(|i|)>,p>|)>>
  </enumerate>

  The first step in <strong|M2> technique takes
  <math|r<rsub|\<sigma\><around*|(|i|)>>n<rsup|2>> and the second takes
  <math|\<kappa\>r<rsub|\<sigma\><around*|(|i|)>><big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>.
  Hence the total complexity is <math|r<rsub|\<sigma\><around*|(|i|)>>n<rsup|2>+\<kappa\>r<rsub|\<sigma\><around*|(|i|)>><big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>.

  <strong|Technique M3>

  <\enumerate>
    <item><strong|Setup> <math|<math-bf|B><rsub|\<sigma\><around*|(|i|)>>=<bs><rsup|-1><A><rsub|\<sigma\><around*|(|i|)>><bs><rsup|-1>>

    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>\<sigma\><around*|(|j|)>>=<around*|\<langle\>|<math-bf|B><rsub|\<sigma\><around*|(|i|)>>,<A><rsub|\<sigma\><around*|(|j|)>>|\<rangle\>>,\<forall\>j\<geq\>i>
  </enumerate>

  The first step in <strong|M3> technique takes <math|n
  f<rsub|\<sigma\><around*|(|i|)>>+n<rsup|3>> and the second takes
  <math|<big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>. The total
  complexity is <math|n \<kappa\>f<rsub|\<sigma\><around*|(|i|)>>+n<rsup|3>+<big|sum><rsub|j\<geq\>i>\<kappa\>f<rsub|\<sigma\><around*|(|j|)>>>.

  <strong|Technique M4>

  <\enumerate>
    <item><strong|Setup> <math|<math-bf|G><rsub|\<sigma\><around*|(|i|)>>=<bs><rsup|-1><A><rsub|\<sigma\><around*|(|i|)>>>

    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>\<sigma\><around*|(|j|)>>=<around*|\<langle\>|<math-bf|G><rsub|\<sigma\><around*|(|i|)>><bs><rsup|-1>,<A><rsub|\<sigma\><around*|(|j|)>>|\<rangle\>>,\<forall\>j\<leq\>i>
  </enumerate>

  The first step in <strong|M4> technique takes <math|n
  f<rsub|\<sigma\><around*|(|i|)>>> and the second takes
  <math|<big|sum><rsub|j\<geq\>i>\<kappa\><around*|(|n+1|)>f<rsub|\<sigma\><around*|(|j|)>>>,
  Hence the total cost is <math|><math|n \<kappa\>f<rsub|\<sigma\><around*|(|i|)>>+<big|sum><rsub|j\<geq\>i>\<kappa\><around*|(|n+1|)>f<rsub|\<sigma\><around*|(|j|)>>>.

  <strong|Technique M5>

  <\enumerate>
    <item><strong|Compute> <math|M<rsub|\<sigma\><around*|(|i|)>\<sigma\><around*|(|j|)>>=<around*|\<langle\>|<bs><rsup|-1><A><rsub|\<sigma\><around*|(|i|)>><bs><rsup|-1>,<A><rsub|\<sigma\><around*|(|j|)>>|\<rangle\>>>
    directly
  </enumerate>

  The above operation takes <math|2\<kappa\><around*|(|f<rsub|\<sigma\><around*|(|i|)>>+1|)>
  <big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>.

  <big-table|<block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<cwith|1|-1|1|-1|cell-valign|c>|<table|<row|<cell|Method>|<cell|Schur>|<cell|<math|<bs><rsup|-1>>>>|<row|<cell|<strong|M1>>|<cell|<math|r<rsub|\<sigma\><around*|(|i|)>><around*|(|n<rsup|2>+2n<rsup|2>|)>+\<kappa\><big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>>|<cell|0>>|<row|<cell|<strong|M2>>|<cell|<math|r<rsub|\<sigma\><around*|(|i|)>><around*|(|n<rsup|2>+\<kappa\><big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>|)>>>|<cell|0>>|<row|<cell|<strong|M3>>|<cell|<math|n\<kappa\>f<rsub|\<sigma\><around*|(|i|)>>+n<rsup|3>+<big|sum><rsub|j\<geq\>i>\<kappa\>f<rsub|\<sigma\><around*|(|j|)>>>>|<cell|<math|n<rsup|3>>>>|<row|<cell|<strong|M4>>|<cell|<math|n\<kappa\>
  f<rsub|\<sigma\><around*|(|i|)>>+<big|sum><rsub|j\<geq\>i>\<kappa\><around*|(|n+1|)>f<rsub|\<sigma\><around*|(|j|)>>>>|<cell|<math|n<rsup|3>>>>|<row|<cell|<strong|M5>>|<cell|<math|\<kappa\><around*|(|2\<kappa\>f<rsub|\<sigma\><around*|(|i|)>>+1|)>
  <big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>>|<cell|<math|n<rsup|3>>>>>>>|Summary
  of the techniques>

  Based on the above discussion, we now formally propose the re-ordering
  heuristic used in DSDP that exploits both sparsity and the low-rank
  structure.

  <\named-algorithm|2. DSDP Re-ordered Schur Matrix Setup<space|14em>>
    <strong|input> <math|<around*|{|<A><rsub|i>|}>,<bs>,\<kappa\>\<geq\>1>

    \;

    <strong|compute> Number of nonzeros <math|f<rsub|1>,\<ldots\>,f<rsub|m>>

    <strong|compute> Permutation <math|<around*|{|\<sigma\><rsub|<around*|(|i|)>>|}>>
    by by <math|<around*|{|f<rsub|1>,\<ldots\>,f<rsub|m>|}>> in descending
    order.

    <strong|compute>

    <space|12em><math|<tabular|<tformat|<table|<row|<cell|d<rsub|i
    1>>|<cell|=>|<cell|r<rsub|\<sigma\><around*|(|i|)>><around*|(|n<rsup|2>+2\<kappa\>n<rsup|2>|)>+\<kappa\><big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>>|<row|<cell|d<rsub|i2>>|<cell|=>|<cell|r<rsub|\<sigma\><around*|(|i|)>><around*|(|n<rsup|2>+\<kappa\><big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>|)>>>|<row|<cell|d<rsub|i3>>|<cell|=>|<cell|n\<kappa\>f<rsub|\<sigma\><around*|(|i|)>>+n<rsup|3>+<big|sum><rsub|j\<geq\>i>\<kappa\>f<rsub|\<sigma\><around*|(|j|)>>>>|<row|<cell|d<rsub|i4>>|<cell|=>|<cell|n\<kappa\>
    f<rsub|\<sigma\><around*|(|i|)>>+<big|sum><rsub|j\<geq\>i>\<kappa\><around*|(|n+1|)>f<rsub|\<sigma\><around*|(|j|)>>>>|<row|<cell|d<rsub|i5>>|<cell|=>|<cell|\<kappa\><around*|(|2\<kappa\>f<rsub|\<sigma\><around*|(|i|)>>+1|)>
    <big|sum><rsub|j\<geq\>i>f<rsub|\<sigma\><around*|(|j|)>>>>>>>>

    <strong|compute>

    \ \ \ \ <space|11em><math|<tabular|<tformat|<table|<row|<cell|\<xi\><rsub|1>>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|m>min<around*|{|d<rsub|i
    1>,d<rsub|i2>|}>>>|<row|<cell|\<xi\><rsub|2>>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|m>min<around*|{|d<rsub|i
    1>,d<rsub|i2>,d<rsub|i3>,d<rsub|i4>,d<rsub|i5>|}>+n<rsup|3>>>>>>>

    <strong|if> <math|\<xi\><rsub|1>\<less\>\<xi\><rsub|2>>

    <space|2em><strong|use M1 & M2>

    <strong|else>

    <space|2em><strong|use M1 & M2 & M3 & M4 & M5>

    <strong|end>

    \;

    <strong|for> <math|i=1,\<ldots\>,m>

    <space|2em><strong|choose> <strong|MX> that minimizes
    <math|d<rsub|\<sigma\><around*|(|i|)>k>>

    <space|2em><strong|for> <math|j=i,\<ldots\>,m>

    <space|4em><strong|compute> <math|M<rsub|\<sigma\><around*|(|i|)>\<sigma\><around*|(|j|)>>>
    by <strong|MX>

    <space|2em><strong|end>

    <strong|end>

    \;

    <strong|output> <math|<M>>
  </named-algorithm>

  <section|Numerical Experiments>

  In this section, we present the preliminary numerical test on two classical
  SDP benchmarks

  <\itemize>
    <item><verbatim|SDPLIB> <cite|borchers1999sdplib>

    A collection of 90 small and medium-sized SDP problems

    <item><verbatim|Mittelmann's Sparse SDP Benchmark>

    A collection of 75 harder SDP instances

    <item><verbatim|Mittelmann's Infeasible SDP Benchmark>

    A collection of 400 (weakly) infeasible SDP problems
  </itemize>

  We present the number of <verbatim|solved> instances by
  <verbatim|Mittelmann> criterion.

  <subsection|SDPLIB>

  The detailed solution log over SDPLIB is left in the supplementart
  materials.

  <\big-table>
    <block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<table|<row|<cell|Instances>|<cell|Number>>|<row|<cell|All>|<cell|90>>|<row|<cell|Solved>|<cell|87/90>>|<row|<cell|Failed>|<cell|3/90>>>>>
  </big-table|SDPLIB>

  <subsection|Mittelmann's Sparse SDP Benchmark>

  <\frame>
    <\small>
      <\code>
        \ \ \ \ \ \ \ \ \ \ \ Instance \ \ \ \ \ \ \ \ CSDP \ \ \ \ \ \ HDSDP
        \ \ \ \ \ \ MOSEK \ \ \ \ \ \ \ SDPA \ \ \ \ \ \ SDPT3
        \ \ \ \ \ SEDUMI\ 

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 2.42
        \ \ \ \ \ \ \ 5.11 \ \ \ \ \ \ \ 1.00 \ \ \ \ \ \ \ 1.71
        \ \ \ \ \ \ \ 1.07 \ \ \ \ \ \ \ 5.29\ 

        ----------------------------------------------------------------------------------------------------

        \ \ \ \ \ \ \ \ \ \ \ \ 1dc.1024 \ \ \ \ \ \ \ 1033
        \ \ \ \ \ \ \ 1955 \ \ \ \ \ \ \ \ 283 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ \ 403 \ \ \ \ \ \ 22087

        \ \ \ \ \ \ \ \ \ \ \ \ 1et.1024 \ \ \ \ \ \ \ \ \ 88
        \ \ \ \ \ \ \ \ 237 \ \ \ \ \ \ \ \ \ 45 \ \ \ \ \ \ \ \ \ 38
        \ \ \ \ \ \ \ \ \ 45 \ \ \ \ \ \ \ 1447

        \ \ \ \ \ \ \ \ \ \ \ \ 1tc.1024 \ \ \ \ \ \ \ \ \ 66
        \ \ \ \ \ \ \ \ 142 \ \ \ \ \ \ \ \ \ 34 \ \ \ \ \ \ \ \ \ 25
        \ \ \ \ \ \ \ \ \ 32 \ \ \ \ \ \ \ \ 937

        \ \ \ \ \ \ \ \ \ \ \ \ 1zc.1024 \ \ \ \ \ \ \ \ 282
        \ \ \ \ \ \ \ \ 396 \ \ \ \ \ \ \ \ \ 95 \ \ \ \ \ \ \ \ 105
        \ \ \ \ \ \ \ \ 116 \ \ \ \ \ \ \ 6110

        \ \ \ \ \ \ \ \ \ \ \ Alh_1.r20 \ \ \ \ \ \ \ 5909 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ 3380 \ \ \ \ \ \ \ 7126 \ \ \ \ \ \ 18001
        \ \ \ \ \ \ \ 2765

        \ \ \ \ \ \ \ \ \ \ \ \ \ BH2.r14 \ \ \ \ \ \ \ \ 171
        \ \ \ \ \ \ \ \ 260 \ \ \ \ \ \ \ \ \ 76 \ \ \ \ \ \ \ \ \ 94
        \ \ \ \ \ \ \ \ 153 \ \ \ \ \ \ \ \ \ 90

        \ \ \ \ \ \ \ \ \ \ \ \ Bex2_1_5 \ \ \ \ \ \ \ \ 301
        \ \ \ \ \ \ \ \ 534 \ \ \ \ \ \ \ \ \ 17 \ \ \ \ \ \ \ \ \ 37
        \ \ \ \ \ \ \ \ 287 \ \ \ \ \ \ \ \ 120

        \ \ \ \ \ \ \ \ \ Bst_jcbpaf2 \ \ \ \ \ \ \ \ 324 \ \ \ \ \ \ \ \ 668
        \ \ \ \ \ \ \ \ \ 27 \ \ \ \ \ \ \ \ \ 42 \ \ \ \ \ \ \ \ 328
        \ \ \ \ \ \ \ \ 107

        \ \ \ \ \ \ \ \ \ \ \ CH2_1.r14 \ \ \ \ \ \ \ \ 147
        \ \ \ \ \ \ \ \ 284 \ \ \ \ \ \ \ \ \ 72 \ \ \ \ \ \ \ \ \ 90
        \ \ \ \ \ \ \ \ 144 \ \ \ \ \ \ \ \ \ 79

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ G40_mb \ \ \ \ \ \ \ \ \ 84
        \ \ \ \ \ \ \ \ \ 77 \ \ \ \ \ \ \ \ 174 \ \ \ \ \ \ \ \ \ 17
        \ \ \ \ \ \ \ \ \ 25 \ \ \ \ \ \ \ \ 487

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ G48_mb \ \ \ \ \ \ \ \ 802
        \ \ \ \ \ \ \ \ 800 \ \ \ \ \ \ \ \ 191 \ \ \ \ \ \ \ \ \ 55
        \ \ \ \ \ \ \ \ \ 49 \ \ \ \ \ \ \ 1539

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ G48mc \ \ \ \ \ \ \ \ \ 35
        \ \ \ \ \ \ \ \ \ 13 \ \ \ \ \ \ \ \ \ 71 \ \ \ \ \ \ \ \ \ 37
        \ \ \ \ \ \ \ \ \ 24 \ \ \ \ \ \ \ \ 641

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ G55mc \ \ \ \ \ \ \ \ 155
        \ \ \ \ \ \ \ \ 363 \ \ \ \ \ \ \ \ 679 \ \ \ \ \ \ \ \ 199
        \ \ \ \ \ \ \ \ 191 \ \ \ \ \ \ \ 4876

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ G59mc \ \ \ \ \ \ \ \ 173
        \ \ \ \ \ \ \ \ 583 \ \ \ \ \ \ \ \ 646 \ \ \ \ \ \ \ \ 261
        \ \ \ \ \ \ \ \ 256 \ \ \ \ \ \ \ 5533

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ G60_mb \ \ \ \ \ \ \ 3779
        \ \ \ \ \ \ \ 2849 \ \ \ \ \ \ \ 7979 \ \ \ \ \ \ \ \ 580
        \ \ \ \ \ \ \ \ 592 \ \ \ \ \ \ 31410

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ G60mc \ \ \ \ \ \ \ 4918
        \ \ \ \ \ \ 40000 \ \ \ \ \ \ \ 8005 \ \ \ \ \ \ \ \ 575
        \ \ \ \ \ \ \ \ 590 \ \ \ \ \ \ 31530

        \ \ \ \ \ \ \ \ \ \ \ \ H30_.r16 \ \ \ \ \ \ \ \ 366
        \ \ \ \ \ \ \ 1175 \ \ \ \ \ \ \ \ 310 \ \ \ \ \ \ \ \ 701
        \ \ \ \ \ \ \ \ 834 \ \ \ \ \ \ \ \ 247

        \ \ \ \ \ \ \ \ \ \ \ \ NH2-.r14 \ \ \ \ \ \ \ \ 148
        \ \ \ \ \ \ \ \ 250 \ \ \ \ \ \ \ \ \ 73 \ \ \ \ \ \ \ \ \ 43
        \ \ \ \ \ \ \ \ 188 \ \ \ \ \ \ \ \ \ 71

        \ \ \ \ \ \ \ \ \ \ \ \ \ NH3.r16 \ \ \ \ \ \ \ \ 370
        \ \ \ \ \ \ \ 1170 \ \ \ \ \ \ \ \ 297 \ \ \ \ \ \ \ \ 718
        \ \ \ \ \ \ \ \ 833 \ \ \ \ \ \ \ \ 247

        \ \ \ \ \ \ \ \ \ \ \ \ NH4*.r18 \ \ \ \ \ \ \ 2223
        \ \ \ \ \ \ \ 3600 \ \ \ \ \ \ \ 1465 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ 5057 \ \ \ \ \ \ \ \ 930

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ biggs \ \ \ \ \ \ \ \ 190
        \ \ \ \ \ \ \ \ \ 14 \ \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ \ \ \ 8 \ \ \ \ \ \ \ \ \ 17

        \ \ \ \ \ \ \ \ \ \ \ broyden25 \ \ \ \ \ \ \ \ 372
        \ \ \ \ \ \ \ 1687 \ \ \ \ \ \ \ \ \ 32 \ \ \ \ \ \ \ \ \ 43
        \ \ \ \ \ \ \ \ 561 \ \ \ \ \ \ \ \ 211

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ buck4 \ \ \ \ \ \ \ \ 107
        \ \ \ \ \ \ \ \ \ 34 \ \ \ \ \ \ \ \ \ 13 \ \ \ \ \ \ \ \ \ \ 4
        \ \ \ \ \ \ \ \ \ 11 \ \ \ \ \ \ \ \ 100

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ buck5 \ \ \ \ \ \ \ \ 222
        \ \ \ \ \ \ \ \ 775 \ \ \ \ \ \ \ \ 400 \ \ \ \ \ \ \ \ \ 37
        \ \ \ \ \ \ \ \ \ 93 \ \ \ \ \ \ \ 2320

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ cancer \ \ \ \ \ \ \ \ 111
        \ \ \ \ \ \ \ \ 245 \ \ \ \ \ \ \ \ \ 75 \ \ \ \ \ \ \ \ \ 19
        \ \ \ \ \ \ \ \ \ 45 \ \ \ \ \ \ \ 1672

        \ \ \ \ \ \ \ \ \ \ \ \ \ checker \ \ \ \ \ \ \ \ 144
        \ \ \ \ \ \ \ \ 115 \ \ \ \ \ \ \ \ \ 72 \ \ \ \ \ \ \ \ 109
        \ \ \ \ \ \ \ \ \ 71 \ \ \ \ \ \ \ 3397

        \ \ \ \ \ \ \ \ \ \ \ \ chs_5000 \ \ \ \ \ \ 40000 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ \ \ \ 4 \ \ \ \ \ \ \ \ \ 21 \ \ \ \ \ \ \ \ \ 24
        \ \ \ \ \ \ \ \ \ 31

        \ \ \ \ \ \ \ \ \ \ \ \ \ cnhil10 \ \ \ \ \ \ \ \ \ 38
        \ \ \ \ \ \ \ \ \ 96 \ \ \ \ \ \ \ \ \ \ 9 \ \ \ \ \ \ \ \ \ 17
        \ \ \ \ \ \ \ \ \ \ 8 \ \ \ \ \ \ \ \ 161

        \ \ \ \ \ \ \ \ \ \ \ \ \ cphil12 \ \ \ \ \ \ \ \ \ 23
        \ \ \ \ \ \ \ \ 418 \ \ \ \ \ \ \ \ \ 22 \ \ \ \ \ \ \ \ 484
        \ \ \ \ \ \ \ \ \ 31 \ \ \ \ \ \ \ 1145

        \ \ \ \ \ \ \ \ \ \ \ dia_patch \ \ \ \ \ \ 40000 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ 3381 \ \ \ \ \ \ 40000 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ 38947

        \ \ \ \ \ \ \ \ \ \ \ \ e*quad*2 \ \ \ \ \ \ \ \ 125
        \ \ \ \ \ \ \ 1386 \ \ \ \ \ \ \ \ \ 14 \ \ \ \ \ \ \ \ \ 10
        \ \ \ \ \ \ \ \ \ 49 \ \ \ \ \ \ \ \ 347

        \ \ \ \ \ \ \ \ \ \ e*stable*2 \ \ \ \ \ \ \ \ \ 84
        \ \ \ \ \ \ \ \ 233 \ \ \ \ \ \ \ \ \ 11 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ \ \ 40 \ \ \ \ \ \ \ \ 293

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ foot \ \ \ \ \ \ \ \ 186
        \ \ \ \ \ \ \ \ 274 \ \ \ \ \ \ \ \ 533 \ \ \ \ \ \ \ \ \ 37
        \ \ \ \ \ \ \ \ \ 32 \ \ \ \ \ \ 40000

        \ \ \ \ \ \ \ \ \ \ \ hamming_8 \ \ \ \ \ \ \ \ 177
        \ \ \ \ \ \ \ \ \ 65 \ \ \ \ \ \ \ \ \ 41 \ \ \ \ \ \ \ \ \ 62
        \ \ \ \ \ \ \ \ \ 65 \ \ \ \ \ \ \ 3811

        \ \ \ \ \ \ \ \ \ \ \ hamming_9 \ \ \ \ \ \ \ 6508 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ 40000 \ \ \ \ \ \ \ 1863 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ 40000

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ hand \ \ \ \ \ \ \ \ \ 46
        \ \ \ \ \ \ \ \ \ 67 \ \ \ \ \ \ \ \ \ 76 \ \ \ \ \ \ \ \ \ \ 6
        \ \ \ \ \ \ \ \ \ \ 8 \ \ \ \ \ \ \ \ 159

        \ \ \ \ \ \ \ \ \ \ \ \ \ ice_2.0 \ \ \ \ \ \ \ 1406
        \ \ \ \ \ \ \ \ 549 \ \ \ \ \ \ \ 4584 \ \ \ \ \ \ \ \ 955
        \ \ \ \ \ \ \ \ 484 \ \ \ \ \ \ 46335

        \ \ \ \ \ \ \ \ \ \ \ \ inc_1200 \ \ \ \ \ \ \ \ 177
        \ \ \ \ \ \ 40000 \ \ \ \ \ \ \ \ \ 98 \ \ \ \ \ \ \ \ \ 12
        \ \ \ \ \ \ \ \ \ 34 \ \ \ \ \ \ \ \ 759

        \ \ \ \ \ \ \ \ \ \ \ \ \ mater-5 \ \ \ \ \ \ \ \ 120
        \ \ \ \ \ \ \ \ 500 \ \ \ \ \ \ \ \ \ \ 3 \ \ \ \ \ \ \ \ \ \ 7
        \ \ \ \ \ \ \ \ \ \ 7 \ \ \ \ \ \ \ \ \ 59

        \ \ \ \ \ \ \ \ \ \ \ \ \ mater-6 \ \ \ \ \ \ \ \ 826
        \ \ \ \ \ \ 40000 \ \ \ \ \ \ \ \ \ \ 8 \ \ \ \ \ \ \ \ \ 17
        \ \ \ \ \ \ \ \ \ 18 \ \ \ \ \ \ \ \ 115

        \ \ \ \ \ \ \ \ \ \ \ neosfbr25 \ \ \ \ \ \ \ \ 601
        \ \ \ \ \ \ \ \ 800 \ \ \ \ \ \ \ \ 113 \ \ \ \ \ \ \ \ 158
        \ \ \ \ \ \ \ \ 455 \ \ \ \ \ \ \ 5685

        \ \ \ \ \ \ \ \ \ neosfbr30e8 \ \ \ \ \ \ \ 2659 \ \ \ \ \ \ \ 2000
        \ \ \ \ \ \ \ \ 607 \ \ \ \ \ \ \ 1057 \ \ \ \ \ \ \ 1698
        \ \ \ \ \ \ \ 8009

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ neu1 \ \ \ \ \ \ \ \ 260
        \ \ \ \ \ \ \ \ 152 \ \ \ \ \ \ \ \ \ \ 6 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ \ \ 73 \ \ \ \ \ \ \ \ \ 55

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ neu1g \ \ \ \ \ \ \ \ 226
        \ \ \ \ \ \ \ \ 110 \ \ \ \ \ \ \ \ \ \ 7 \ \ \ \ \ \ \ \ \ \ 8
        \ \ \ \ \ \ \ \ \ 50 \ \ \ \ \ \ \ \ \ 44

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ neu2 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ \ 157 \ \ \ \ \ \ \ \ \ 10 \ \ \ \ \ \ \ \ \ \ 1
        \ \ \ \ \ \ \ \ \ 73 \ \ \ \ \ \ \ \ \ 58

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ neu2c \ \ \ \ \ \ \ \ 746
        \ \ \ \ \ \ \ \ 395 \ \ \ \ \ \ \ \ \ \ 8 \ \ \ \ \ \ \ \ \ 17
        \ \ \ \ \ \ \ \ 184 \ \ \ \ \ \ \ \ 105

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ neu2g \ \ \ \ \ \ \ \ 106
        \ \ \ \ \ \ \ \ 132 \ \ \ \ \ \ \ \ \ 31 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ \ \ 59 \ \ \ \ \ \ \ \ \ 58

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ neu3 \ \ \ \ \ \ \ \ 387
        \ \ \ \ \ \ \ 2073 \ \ \ \ \ \ \ \ \ 44 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ \ 105 \ \ \ \ \ \ \ \ 571

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ neu3g \ \ \ \ \ \ \ \ 381
        \ \ \ \ \ \ \ 3609 \ \ \ \ \ \ \ \ \ 49 \ \ \ \ \ \ \ \ 254
        \ \ \ \ \ \ \ \ \ 94 \ \ \ \ \ \ \ \ 901

        \ \ \ \ \ \ \ \ \ \ \ \ \ p_auss2 \ \ \ \ \ \ \ 1885
        \ \ \ \ \ \ \ \ 519 \ \ \ \ \ \ \ 5948 \ \ \ \ \ \ \ 1276
        \ \ \ \ \ \ \ \ 640 \ \ \ \ \ \ \ \ 895

        \ \ \ \ \ \ \ \ \ \ prob_2_4_0 \ \ \ \ \ \ \ \ \ 82
        \ \ \ \ \ \ \ 1215 \ \ \ \ \ \ 40000 \ \ \ \ \ \ \ \ \ 32
        \ \ \ \ \ \ \ \ \ 14 \ \ \ \ \ \ \ \ 156

        \ \ \ \ \ \ \ \ \ \ prob_2_4_1 \ \ \ \ \ \ \ \ \ 97
        \ \ \ \ \ \ \ \ 111 \ \ \ \ \ \ 40000 \ \ \ \ \ \ \ \ \ 82
        \ \ \ \ \ \ \ \ \ 11 \ \ \ \ \ \ \ \ 152

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ rabmo \ \ \ \ \ \ \ \ \ 52
        \ \ \ \ \ \ 40000 \ \ \ \ \ \ \ \ \ \ 8 \ \ \ \ \ \ \ \ \ 16
        \ \ \ \ \ \ \ \ \ 65 \ \ \ \ \ \ \ \ 169

        \ \ \ \ \ \ \ \ \ \ \ \ \ reimer5 \ \ \ \ \ \ \ \ 292
        \ \ \ \ \ \ 40000 \ \ \ \ \ \ \ \ \ 30 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ 40000 \ \ \ \ \ \ \ \ 436

        \ \ \ \ \ \ \ \ \ \ \ \ \ r1_2000 \ \ \ \ \ \ \ \ \ 49
        \ \ \ \ \ \ \ \ 101 \ \ \ \ \ \ \ \ 333 \ \ \ \ \ \ \ \ \ 20
        \ \ \ \ \ \ \ \ \ 20 \ \ \ \ \ \ \ \ 419

        \ \ \ \ \ \ \ \ \ \ \ \ ros_2000 \ \ \ \ \ \ \ \ 855
        \ \ \ \ \ \ \ 1556 \ \ \ \ \ \ \ \ \ \ 1 \ \ \ \ \ \ \ \ \ \ 2
        \ \ \ \ \ \ \ \ \ \ 2 \ \ \ \ \ \ \ \ \ 12

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ rose15 \ \ \ \ \ \ \ \ 236
        \ \ \ \ \ \ \ \ 141 \ \ \ \ \ \ \ \ \ \ 5 \ \ \ \ \ \ 40000
        \ \ \ \ \ \ \ \ \ 30 \ \ \ \ \ \ \ \ \ 87

        \ \ \ \ \ \ \ \ \ sensor_1000 \ \ \ \ \ \ \ \ 264 \ \ \ \ \ \ \ \ 290
        \ \ \ \ \ \ \ \ \ 24 \ \ \ \ \ \ \ \ \ 18 \ \ \ \ \ \ \ \ \ 32
        \ \ \ \ \ \ \ \ 414

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ shmup4 \ \ \ \ \ \ \ \ 217
        \ \ \ \ \ \ \ \ 273 \ \ \ \ \ \ \ \ 222 \ \ \ \ \ \ \ \ \ 30
        \ \ \ \ \ \ \ \ 105 \ \ \ \ \ \ \ 1580

        \ \ \ \ \ \ \ \ \ \ \ \ \ \ shmup5 \ \ \ \ \ \ \ \ 855
        \ \ \ \ \ \ \ 2590 \ \ \ \ \ \ \ 3247 \ \ \ \ \ \ \ \ 268
        \ \ \ \ \ \ \ \ 789 \ \ \ \ \ \ 21601

        \ \ \ \ \ \ \ \ \ \ \ \ \ spar060 \ \ \ \ \ \ \ \ \ 58
        \ \ \ \ \ \ \ 1092 \ \ \ \ \ \ \ \ \ 18 \ \ \ \ \ \ \ \ \ 17
        \ \ \ \ \ \ \ \ \ 39 \ \ \ \ \ \ \ 1162

        ---------------------------------------------------------------------------------------------------
      </code>
    </small>
  </frame>

  <\big-table>
    <block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<table|<row|<cell|Instances>|<cell|Number>>|<row|<cell|All>|<cell|75>>|<row|<cell|Solved>|<cell|63/75>>|<row|<cell|Failed>|<cell|12/75>>>>>
  </big-table|Mittelmann's Sparse SDP Benchmark>

  <subsection|Mittelmann's Infeasible SDP Benchmark>

  <big-table|<block|<tformat|<cwith|1|5|1|6|cell-halign|c>|<cwith|1|5|1|1|cell-halign|c>|<cwith|1|5|1|1|cell-halign|c>|<table|<row|<cell|<verbatim|HDSDP>>|<cell|<code|MOSEK>>|<cell|<code|SDPA>>|<cell|<code|SDPT3>>|<cell|<code|SeDuMi>>|<cell|<code|PENSDP>>>|<row|<cell|83>|<cell|15>|<cell|100>|<cell|10>|<cell|87>|<cell|0>>|<row|<cell|100>|<cell|100>|<cell|100>|<cell|100>|<cell|100>|<cell|0>>|<row|<cell|8>|<cell|4>|<cell|30>|<cell|4>|<cell|25>|<cell|0>>|<row|<cell|97>|<cell|84>|<cell|100>|<cell|96>|<cell|100>|<cell|0>>>>>|Mittelmann's
  Infeasible SDP Benchmark>

  From the above experiments, <verbatim|HDSDP>, in terms of robustness and
  solution, is comparable to the primal-dual solvers. Currently
  <verbatim|HDSDP> is still under active development and better efficiency is
  expected to be achieved in the future.

  <section|Conclusions>

  In this paper we present a dual-scaling algorithm based on homogeneous
  self-dual model. The resultant solver, <verbatim|HDSDP>, is entailed with a
  more robust infeasibility detection module. Moreover, several newly added
  computational techniques further enhance the efficiency of the solver. The
  first formal version of <verbatim|HDSDP> is expected to be released in June
  and please contact the author for more details.

  <section|Acknowledgement>

  I sincerely appreciate the instructions from\ 

  <\itemize>
    <item>Yinyu Ye

    <item>Dongdong Ge

    <item>Qi Huangfu
  </itemize>

  who provided instructions and helpful suggestions during the development of
  the solver.

  <\bibliography*|bib|tm-alpha|sdpref|References>
    <\bib-list|16>
      <bibitem*|ApS19><label|bib-aps2019mosek>Mosek ApS. <newblock>Mosek
      optimization toolbox for matlab. <newblock><with|font-shape|italic|User's
      Guide and Reference Manual, Version>, 4, 2019.<newblock>

      <bibitem*|Bor99a><label|bib-borchers1999csdp>Brian Borchers.
      <newblock>Csdp, ac library for semidefinite programming.
      <newblock><with|font-shape|italic|Optimization methods and Software>,
      11(1-4):613\U623, 1999.<newblock>

      <bibitem*|Bor99b><label|bib-borchers1999sdplib>Brian Borchers.
      <newblock>Sdplib 1.2, a library of semidefinite programming test
      problems. <newblock><with|font-shape|italic|Optimization Methods and
      Software>, 11(1-4):683\U690, 1999.<newblock>

      <bibitem*|BY01><label|bib-benson2001dsdp3>Steven<nbsp>J
      Benson<localize| and >Yinyu Ye. <newblock>Dsdp3: dual scaling algorithm
      for general positive semidefinite programming.
      <newblock><localize|Technical Report>, Argonne National Lab., IL (US),
      2001.<newblock>

      <bibitem*|BY08><label|bib-benson2008algorithm>Steven<nbsp>J
      Benson<localize| and >Yinyu Ye. <newblock>Algorithm 875:
      dsdp5\Vsoftware for semidefinite programming.
      <newblock><with|font-shape|italic|ACM Transactions on Mathematical
      Software (TOMS)>, 34(3):1\U20, 2008.<newblock>

      <bibitem*|FKN97><label|bib-fujisawa1997exploiting>Katsuki Fujisawa,
      Masakazu Kojima<localize|, and >Kazuhide Nakata. <newblock>Exploiting
      sparsity in primal-dual interior-point methods for semidefinite
      programming. <newblock><with|font-shape|italic|Mathematical
      Programming>, 79(1):235\U253, 1997.<newblock>

      <bibitem*|GW95><label|bib-goemans1995improved>Michel<nbsp>X
      Goemans<localize| and >David<nbsp>P Williamson. <newblock>Improved
      approximation algorithms for maximum cut and satisfiability problems
      using semidefinite programming. <newblock><with|font-shape|italic|Journal
      of the ACM (JACM)>, 42(6):1115\U1145, 1995.<newblock>

      <bibitem*|LR05><label|bib-laurent2005semidefinite>Monique
      Laurent<localize| and >Franz Rendl. <newblock>Semidefinite programming
      and integer programming. <newblock><with|font-shape|italic|Handbooks in
      Operations Research and Management Science>, 12:393\U514,
      2005.<newblock>

      <bibitem*|PS98a><label|bib-potra1998superlinearly>Florian<nbsp>A
      Potra<localize| and >Rongqin Sheng. <newblock>A superlinearly
      convergent primal-dual infeasible-interior-point algorithm for
      semidefinite programming. <newblock><with|font-shape|italic|SIAM
      Journal on Optimization>, 8(4):1007\U1028, 1998.<newblock>

      <bibitem*|PS98b><label|bib-potra1998homogeneous>Florian<nbsp>A
      Potra<localize| and >Rongqin Sheng. <newblock>On homogeneous
      interrior-point algorithms for semidefinite programming.
      <newblock><with|font-shape|italic|Optimization Methods and Software>,
      9(1-3):161\U184, 1998.<newblock>

      <bibitem*|Stu99><label|bib-sturm1999using>Jos<nbsp>F Sturm.
      <newblock>Using sedumi 1.02, a matlab toolbox for optimization over
      symmetric cones. <newblock><with|font-shape|italic|Optimization methods
      and software>, 11(1-4):625\U653, 1999.<newblock>

      <bibitem*|TTT99><label|bib-toh1999sdpt3>Kim-Chuan Toh, Michael<nbsp>J
      Todd<localize|, and >Reha<nbsp>H Tütüncü. <newblock>Sdpt3\Va matlab
      software package for semidefinite programming, version 1.3.
      <newblock><with|font-shape|italic|Optimization methods and software>,
      11(1-4):545\U581, 1999.<newblock>

      <bibitem*|VB96><label|bib-vandenberghe1996semidefinite>Lieven
      Vandenberghe<localize| and >Stephen Boyd. <newblock>Semidefinite
      programming. <newblock><with|font-shape|italic|SIAM review>,
      38(1):49\U95, 1996.<newblock>

      <bibitem*|VB99><label|bib-vandenberghe1999applications>Lieven
      Vandenberghe<localize| and >Stephen Boyd. <newblock>Applications of
      semidefinite programming. <newblock><with|font-shape|italic|Applied
      Numerical Mathematics>, 29(3):283\U299, 1999.<newblock>

      <bibitem*|XHY96><label|bib-xu1996simplified>Xiaojie Xu, Pi-Fang
      Hung<localize|, and >Yinyu Ye. <newblock>A simplified homogeneous and
      self-dual linear programming algorithm and its implementation.
      <newblock><with|font-shape|italic|Annals of Operations Research>,
      62(1):151\U171, 1996.<newblock>

      <bibitem*|YFK03><label|bib-yamashita2003implementation>Makoto
      Yamashita, Katsuki Fujisawa<localize|, and >Masakazu Kojima.
      <newblock>Implementation and evaluation of sdpa 6.0 (semidefinite
      programming algorithm 6.0). <newblock><with|font-shape|italic|Optimization
      Methods and Software>, 18(4):491\U505, 2003.<newblock>
    </bib-list>
  </bibliography*>
</body>

<\initial>
  <\collection>
    <associate|info-flag|minimal>
    <associate|page-medium|paper>
    <associate|page-screen-margin|true>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|3.1.1|6>>
    <associate|auto-11|<tuple|3.1.2|6>>
    <associate|auto-12|<tuple|3.1.3|6>>
    <associate|auto-13|<tuple|3.2|6>>
    <associate|auto-14|<tuple|3.2.1|6>>
    <associate|auto-15|<tuple|3.2.2|6>>
    <associate|auto-16|<tuple|3.2.3|7>>
    <associate|auto-17|<tuple|1|8>>
    <associate|auto-18|<tuple|4|9>>
    <associate|auto-19|<tuple|4.1|9>>
    <associate|auto-2|<tuple|2|2>>
    <associate|auto-20|<tuple|2|9>>
    <associate|auto-21|<tuple|4.2|10>>
    <associate|auto-22|<tuple|3|11>>
    <associate|auto-23|<tuple|4.3|11>>
    <associate|auto-24|<tuple|4|11>>
    <associate|auto-25|<tuple|5|11>>
    <associate|auto-26|<tuple|6|11>>
    <associate|auto-27|<tuple|<with|mode|<quote|math>|\<bullet\>>|11>>
    <associate|auto-3|<tuple|2.1|2>>
    <associate|auto-4|<tuple|2.2|2>>
    <associate|auto-5|<tuple|2.3|2>>
    <associate|auto-6|<tuple|2.4|4>>
    <associate|auto-7|<tuple|3|5>>
    <associate|auto-8|<tuple|1|5>>
    <associate|auto-9|<tuple|3.1|5>>
    <associate|bib-aps2019mosek|<tuple|ApS19|11>>
    <associate|bib-benson2001dsdp3|<tuple|BY01|11>>
    <associate|bib-benson2008algorithm|<tuple|BY08|11>>
    <associate|bib-borchers1999csdp|<tuple|Bor99a|11>>
    <associate|bib-borchers1999sdplib|<tuple|Bor99b|11>>
    <associate|bib-fujisawa1997exploiting|<tuple|FKN97|11>>
    <associate|bib-goemans1995improved|<tuple|GW95|12>>
    <associate|bib-laurent2005semidefinite|<tuple|LR05|12>>
    <associate|bib-potra1998homogeneous|<tuple|PS98b|12>>
    <associate|bib-potra1998superlinearly|<tuple|PS98a|12>>
    <associate|bib-sturm1999using|<tuple|Stu99|12>>
    <associate|bib-toh1999sdpt3|<tuple|TTT99|12>>
    <associate|bib-vandenberghe1996semidefinite|<tuple|VB96|12>>
    <associate|bib-vandenberghe1999applications|<tuple|VB99|12>>
    <associate|bib-xu1996simplified|<tuple|XHY96|12>>
    <associate|bib-yamashita2003implementation|<tuple|YFK03|12>>
    <associate|sec2|<tuple|2|2>>
    <associate|sec3|<tuple|3|5>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      laurent2005semidefinite

      vandenberghe1999applications

      vandenberghe1996semidefinite

      goemans1995improved

      laurent2005semidefinite

      benson2008algorithm

      aps2019mosek

      sturm1999using

      toh1999sdpt3

      borchers1999csdp

      yamashita2003implementation

      potra1998superlinearly

      potra1998homogeneous

      aps2019mosek

      benson2001dsdp3

      benson2001dsdp3

      xu1996simplified

      fujisawa1997exploiting

      borchers1999sdplib
    </associate>
    <\associate|figure>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||The inner procedure
      of HDSDP>|<pageref|auto-8>>
    </associate>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Summary of the
      techniques>|<pageref|auto-17>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|2>||SDPLIB>|<pageref|auto-20>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|3>||Mittelmann's Sparse
      SDP Benchmark>|<pageref|auto-22>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|4>||Mittelmann's
      Infeasible SDP Benchmark>|<pageref|auto-24>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Homogeneous
      Dual Scaling Algorithm > <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Notations>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>Standard Form
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1tab>|2.3<space|2spc>HSD and Dual-scaling
      Algorithm <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|1tab>|2.4<space|2spc>Barrier and Centrality
      Charactrization <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Implementation
      Details> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Presolve
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|2tab>|3.1.1<space|2spc><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Rank-one
      detection>* <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|2tab>|3.1.2<space|2spc>Sparse eigen-decomposition
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|2tab>|3.1.3<space|2spc>Coefficient scaling
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Dual infeasibility
      elimination <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|2tab>|3.2.1<space|2spc>Initialization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|2tab>|3.2.2<space|2spc>Damped Dual Infeasibility
      Elimination* <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|2tab>|3.2.3<space|2spc>Schur complement tricks*
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Numerical
      Experiments> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.5fn>

      <with|par-left|<quote|1tab>|4.1<space|2spc>SDPLIB
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19>>

      <with|par-left|<quote|1tab>|4.2<space|2spc>Mittelmann's Sparse SDP
      Benchmark <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>>

      <with|par-left|<quote|1tab>|4.3<space|2spc>Mittelmann's Infeasible SDP
      Benchmark <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Conclusions>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-25><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Acknowledge>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|References>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-27><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>