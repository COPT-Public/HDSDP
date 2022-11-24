<TeXmacs|2.1>

<style|<tuple|generic|triangle-list>>

<\body>
  <doc-data|<doc-title|On the design of a dual potential reduction
  solver>|<doc-author|<author-data|<\author-affiliation>
    \;

    HDSDP Developers' Group

    \;

    \;

    <date|>
  </author-affiliation>>>>

  In this note we describe the implementation of a dual potential reduction
  solver that exploits either the embedding or the big-<math|M> potential
  reduction method to solve\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|x>>|<cell|<around*|\<langle\>|c,x|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|\<cal-A\> x=b>|<cell|>>|<row|<cell|>|<cell|x\<succeq\><rsub|\<cal-K\>>0>|<cell|>>>>
  </eqnarray*>

  via its dual

  <\eqnarray*>
    <tformat|<table|<row|<cell|max<rsub|y,s>>|<cell|b<rsup|\<top\>>y>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|\<cal-A\><rsup|\<ast\>>y+s=c>|<cell|>>|<row|<cell|>|<cell|s\<succeq\><rsub|\<cal-K\><rsup|\<ast\>>>0.>|<cell|>>>>
  </eqnarray*>

  The main solver contains the following components

  <\itemize>
    <item>Solver/data

    Interface: set data, set parameter, set solution, optimize, get solution

    Sparse, dense, rank-one

    <item>Algorithm

    HSD/infeasible start/dual potential reduction

    Presolve, phase A, phase B

    Schur complement setup*, potential line-search, barrier line-search

    <item>Linear algebra

    sparse, dense, low rank; eigen, trace, decomposition (Cholesky)

    Lanczos, conjugate gradient, block buffer computation*

    <item>Other utilities

    Parameter tuner, IO and things like that
  </itemize>

  <section|SDP Data structures>

  <subsection|Factorized data>

  The following structure stores the eigen-decomposition of a data matrix
  <math|A=<big|sum><rsub|i=1><rsup|r>\<lambda\><rsub|i>u<rsub|i>u<rsub|i><rsup|\<top\>>>.
  The structure should support the following operations.

  <\itemize>
    <item><math|<around*|\<langle\>|A,B|\<rangle\>>=<big|sum><rsub|i=1><rsup|r>\<lambda\><rsub|i>u<rsub|i><rsup|\<top\>>B
    u<rsub|i>>

    <item><math|B=S<rsup|-1>A S<rsup|-1>=<big|sum><rsub|i=1><rsup|r>\<lambda\><rsub|i><around*|(|S<rsup|-1>u<rsub|i>|)><around*|(|S<rsup|-1>u<rsub|i>|)><rsup|\<top\>>>

    In this case the LHS serves as a buffer
  </itemize>

  <\listing>
    <\verbatim>
      typedef struct {
    </verbatim>

    \;

    <\verbatim>
      \ \ \ \ int \ \ \ \ rank;
    </verbatim>

    <space|2em><verbatim|double *evals;>

    <space|2em><verbatim|double *evecs;>

    \;

    <verbatim|} eigFactor;>
  </listing>

  <subsection|SDP coefficient matrix>

  <with|color|red|NOTE: Only lower triangular is stored.>\ 

  We use the following structures to store <math|A> and <math|C> matrices
  from SDP coefficients. They should support the following functionalities

  <\itemize>
    <item><math|B\<leftarrow\>\<alpha\>A+B>

    <item><math|<around*|\<langle\>|A<rsub|i>,A<rsub|j>|\<rangle\>>> (TODO)

    <item><math|<around*|\<\|\|\>|A|\<\|\|\>><rsub|F>>

    <item><math|<big|sum><rsub|i j><around*|\||a<rsub|i j>|\|>>

    <item><math|A\<leftarrow\>\<alpha\>A>

    <item><verbatim|[V, e] = eig(A)> (TODO)

    <item><verbatim|full(A)>
  </itemize>

  <\listing>
    <\verbatim>
      typedef struct {
    </verbatim>

    \;

    <space|2em><verbatim|int \ \ \ \ \ \ nCol;>

    <space|2em><verbatim|void \ \ \ \ \ *dataMat;>

    <space|2em><verbatim|eigFactor *eig;>

    <space|2em><verbatim|void \ \ \ \ \ (*dataMataApB) \ (void *, double,
    void *);>

    <space|2em><verbatim|void \ \ \ \ \ (*dataMatScal) \ (void *, double);>

    <space|2em><verbatim|double \ \ \ (*dataMatNorm) \ (void *, int);>

    <space|2em><verbatim|int \ \ \ \ \ \ (*dataMatEig) \ \ (void *, void
    **);>

    <space|2em><verbatim|int \ \ \ \ \ \ (*dataMatGetNnz)(void *);>

    <space|2em><verbatim|void \ \ \ \ \ (*dataMatDump) \ (void *, double *);>

    \;

    <verbatim|} sdpCoeffMat;>
  </listing>

  <subsubsection|Sparse matrix>

  <\listing>
    <\verbatim>
      typedef struct {
    </verbatim>

    \;

    <space|1em><verbatim|int \ \ \ \ nSDPCol;>

    <space|1em><verbatim|int \ \ \ \ nTriMatElem;>

    <space|1em><verbatim|int \ \ \ *triMatCol;>

    <space|1em><verbatim|int \ \ \ *triMatRow;>

    <space|1em><verbatim|double *triMatElem;>

    \;

    <verbatim|} sdpSparseData;>
  </listing>

  <subsubsection|Dense matrix>

  <\listing>
    <\verbatim>
      typedef struct {
    </verbatim>

    \;

    <space|1em><verbatim|int \ \ \ \ nSDPCol;>

    <space|1em><verbatim|double *dsMatElem;>

    \;

    <verbatim|} sdpDenseData;>
  </listing>

  <subsubsection|Rank-one sparse matrix>

  <\listing>
    <\verbatim>
      typedef struct {
    </verbatim>

    \;

    <space|1em><verbatim|int \ \ \ \ \ nSDPCol;>

    <space|1em><verbatim|int \ \ \ \ \ nSpR1FactorElem;>

    <space|1em><verbatim|int \ \ \ \ *spR1MatIdx;>

    <space|1em><verbatim|double \ *spR1MatElem;>

    \;

    <verbatim|} sdpRankOneSpData;>
  </listing>

  <subsubsection|Rank-one dense matrix>

  <\listing>
    <\verbatim>
      typedef struct {
    </verbatim>

    \;

    <space|1em><verbatim|int \ \ \ \ nSDPCol;>

    <space|1em><verbatim|double *r1MatFactor;>

    \;

    <verbatim|} sdpRankOneSpData;>
  </listing>

  <\subsection>
    SDP variable and step
  </subsection>

  We use the following structures to store <math|S>.

  <\itemize>
    <item><math|S<rsup|-1>>

    <item><verbatim|L = chol(S)>

    <item><verbatim|L \\ z, L' \\ z>

    <item><math|y\<leftarrow\>\<alpha\>S x+y>
  </itemize>

  coming<text-dots>

  <subsection|Schur complement matrix>

  <section|Contribution and formats>

  <\itemize>
    <item>Indentation, bracket

    Default as in Xcode, following the samples below

    <item>Doxygen string and comments

    Using <verbatim|@file, @brief, /* */>

    <item>Function with <verbatim|void> return value should
    <verbatim|return;>

    <item>Name style

    Bottom-level routine: <verbatim|extern void csp_Axpby>

    Medium-level routine: <verbatim|hdsdpSpMatTrace>

    <item>Use <verbatim|assert> whenever necessary

    <item>Static before extern

    <item><text-dots>
  </itemize>

  <\frame>
    <\cpp-code>
      static int pdsCreate( void **pldl, int n ) {

      \ \ \ \ 

      \ \ \ \ int retcode = RETCODE_OK;

      \ \ \ \ pds_linsys *pds = NULL;

      \ \ \ \ 

      \ \ \ \ POTLP_INIT(pds, pds_linsys, 1);

      \ \ \ \ 

      \ \ \ \ if ( !pds ) {

      \ \ \ \ \ \ \ \ retcode = RETCODE_FAILED;

      \ \ \ \ \ \ \ \ goto exit_cleanup;

      \ \ \ \ }

      \ \ \ \ 

      \ \ \ \ \ pds-\<gtr\>n = n;

      \ \ \ \ *pldl = pds;

      \ \ \ \ 

      \ \ \ \ /* Initialize pardiso */

      \ \ \ \ POTLP_ZERO(pds-\<gtr\>pt, void *, 64);

      \ \ \ \ POTLP_ZERO(pds-\<gtr\>iparm, int, 64);

      \ \ \ \ 

      \ \ \ \ int mtype = PARDISO_SYM_INDEFINITE;

      \ \ \ \ pardisoinit(pds-\<gtr\>pt, &mtype, pds-\<gtr\>iparm);

      \ \ \ \ 

      \ \ \ \ set_pardiso_param(pds-\<gtr\>iparm, PARDISO_PARAM_NONDEFAULT,
      1);

      \ \ \ \ set_pardiso_param(pds-\<gtr\>iparm, PARDISO_PARAM_SYMBOLIC,
      PARDISO_PARAM_SYMBOLIC_MMD);

      \ \ \ \ set_pardiso_param(pds-\<gtr\>iparm, PARDISO_PARAM_PERTURBATION,
      3);

      \ \ \ \ set_pardiso_param(pds-\<gtr\>iparm, PARDISO_PARAM_INPLACE, 1);

      \ \ \ \ set_pardiso_param(pds-\<gtr\>iparm, PARDISO_PARAM_INDEX,
      PARDISO_PARAM_INDEX_C);

      \ \ \ \ 

      exit_cleanup:

      \ \ \ \ return retcode;

      }
    </cpp-code>
  </frame>

  \;

  <\frame>
    <\cpp-code>
      extern void potVecScal( pot_vec *pVexX, double sVal ) {

      \ \ \ \ 

      \ \ \ \ scal(&pVexX-\<gtr\>n, &sVal, pVexX-\<gtr\>x,
      &potIntConstantOne);

      \ \ \ \ 

      \ \ \ \ if ( pVexX -\<gtr\>nrm != -1.0 ) {

      \ \ \ \ \ \ \ \ pVexX-\<gtr\>nrm = pVexX-\<gtr\>nrm * fabs(sVal);

      \ \ \ \ }

      \ \ \ \ 

      \ \ \ \ return;

      }
    </cpp-code>
  </frame>
</body>

<\initial>
  <\collection>
    <associate|font-base-size|9>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|2|3>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|1.2|2>>
    <associate|auto-4|<tuple|1.2.1|2>>
    <associate|auto-5|<tuple|1.2.2|3>>
    <associate|auto-6|<tuple|1.2.3|3>>
    <associate|auto-7|<tuple|1.2.4|3>>
    <associate|auto-8|<tuple|1.3|3>>
    <associate|auto-9|<tuple|1.4|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>SDP
      Data structures> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Factorized data
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>SDP coefficient matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|2tab>|1.2.1<space|2spc>Sparse matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|2tab>|1.2.2<space|2spc>Dense matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|1.2.3<space|2spc>Rank-one sparse matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|2tab>|1.2.4<space|2spc>Rank-one dense matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|1tab>|1.3<space|2spc>SDP variable and step
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|1tab>|1.4<space|2spc>Schur complement matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Contribution
      and formats> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>