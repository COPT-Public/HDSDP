<TeXmacs|2.1.1>

<style|<tuple|manual|padded-paragraphs|compact-list>>

<\body>
  <\table-of-contents|toc>
    1<space|2spc>Introduction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-1>

    <with|par-left|1tab|1.1<space|2spc>Release note
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-2>>

    2<space|2spc>Installation and requirement
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-3>

    <with|par-left|1tab|2.1<space|2spc>Package dependence
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-4>>

    <with|par-left|4tab|Linking with LAPACK and BLAS.
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-5><vspace|0.15fn>>

    <with|par-left|1tab|2.2<space|2spc>Building the HDSDP library and
    SDPA/MPS solver <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-6>>

    <with|par-left|1tab|2.3<space|2spc>Using the executable binary
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-7>>

    3<space|2spc>Data and solver interface
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-8>

    <with|par-left|1tab|3.1<space|2spc>Notations
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-9>>

    <with|par-left|1tab|3.2<space|2spc>Return code
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-10>>

    <with|par-left|1tab|3.3<space|2spc>Data format
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-11>>

    <with|par-left|2tab|3.3.1<space|2spc>SDP cone data
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-12>>

    <with|par-left|2tab|3.3.2<space|2spc>LP cone data
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-13>>

    <with|par-left|4tab|Scalar bound cone.
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-14><vspace|0.15fn>>

    <with|par-left|1tab|3.4<space|2spc>Data interface
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-15>>

    <with|par-left|2tab|3.4.1<space|2spc>Create user data structure
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-16>>

    <with|par-left|2tab|3.4.2<space|2spc>Set cone data
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-17>>

    <with|par-left|2tab|3.4.3<space|2spc>Determine cone type
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-18>>

    <with|par-left|2tab|3.4.4<space|2spc>Clear user data
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-19>>

    <with|par-left|2tab|3.4.5<space|2spc>Destroy user data structure
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-20>>

    <with|par-left|1tab|3.5<space|2spc>Solver interface
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-21>>

    <with|par-left|2tab|3.5.1<space|2spc>Create the solver structure
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-22>>

    <with|par-left|2tab|3.5.2<space|2spc>Solver initialization
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-23>>

    <with|par-left|2tab|3.5.3<space|2spc>Set conic data
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-24>>

    <with|par-left|2tab|3.5.4<space|2spc>Set dual objective
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-25>>

    <with|par-left|2tab|3.5.5<space|2spc>Set dual starting point
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-26>>

    <with|par-left|2tab|3.5.6<space|2spc>Set integer parameter
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-27>>

    <with|par-left|2tab|3.5.7<space|2spc>Set double parameter
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-28>>

    <with|par-left|2tab|3.5.8<space|2spc>Optimization
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-29>>

    <with|par-left|2tab|3.5.9<space|2spc>Get dual solution
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-30>>

    <with|par-left|2tab|3.5.10<space|2spc>Get cone solution
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-31>>

    <with|par-left|2tab|3.5.11<space|2spc>Check solution and compute solution
    metrics <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-32>>

    <with|par-left|2tab|3.5.12<space|2spc>Release solver memory
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-33>>

    <with|par-left|2tab|3.5.13<space|2spc>Destroy solver structure
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-34>>

    4<space|2spc>Parameters <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-35>

    <with|par-left|1tab|4.1<space|2spc>Integer parameters
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-36>>

    <with|par-left|2tab|4.1.1<space|2spc>INT_PARAM_CORRECTORA
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-37>>

    <with|par-left|2tab|4.1.2<space|2spc>INT_PARAM_CORRECTORB
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-38>>

    <with|par-left|2tab|4.1.3<space|2spc>INT_PARAM_MAXITER
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-39>>

    <with|par-left|2tab|4.1.4<space|2spc>INT_PARAM_PSDP
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-40>>

    <with|par-left|2tab|4.1.5<space|2spc>INT_PARAM_PRELEVEL
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-41>>

    <with|par-left|2tab|4.1.6<space|2spc>INT_PARAM_THREADS
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-42>>

    <with|par-left|1tab|4.2<space|2spc>Double parameters
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-43>>

    <with|par-left|2tab|4.2.1<space|2spc>DBL_PARAM_RELFEASTOL
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-44>>

    <with|par-left|2tab|4.2.2<space|2spc>DBL_PARAM_RELOPTTOL
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-45>>

    <with|par-left|2tab|4.2.3<space|2spc>DBL_PARAM_ABSFEASTOL
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-46>>

    <with|par-left|2tab|4.2.4<space|2spc>DBL_PARAM_ABSOPTTOL
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-47>>

    <with|par-left|2tab|4.2.5<space|2spc>DBL_PARAM_TIMELIMIT
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-48>>

    <with|par-left|2tab|4.2.6<space|2spc>DBL_PARAM_POTRHOVAL
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-49>>

    <with|par-left|2tab|4.2.7<space|2spc>DBL_PARAM_HSDGAMMA
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-50>>

    <with|par-left|2tab|4.2.8<space|2spc>DBL_PARAM_DUALBOX_LOW
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-51>>

    <with|par-left|2tab|4.2.9<space|2spc>DBL_PARAM_DUALBOX_UP
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-52>>

    <with|par-left|2tab|4.2.10<space|2spc>DBL_PARAM_BARMUSTART
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-53>>

    <with|par-left|2tab|4.2.11<space|2spc>DBL_PARAM_DUALSTART
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-54>>

    <with|par-left|2tab|4.2.12<space|2spc>DBL_PARAM_POBJSTART
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-55>>

    <with|par-left|2tab|4.2.13<space|2spc>DBL_PARAM_TRXESTIMATE
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-56>>

    <with|par-left|2tab|4.2.14<space|2spc>DBL_PARAM_PRECORDACC
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-57>>

    5<space|2spc>Logging and examples <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-58>

    <with|par-left|1tab|5.1<space|2spc>Presolving
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-59>>

    <with|par-left|1tab|5.2<space|2spc>Solution logging
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-60>>

    <with|par-left|1tab|5.3<space|2spc>Solution quality check
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-61>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|References>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-62><vspace|0.5fn>
  </table-of-contents>

  <new-page*><section|Introduction>

  HDSDP is a numerical software for general-purpose semidefinite programming
  problems (SDP).

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|X>>|<cell|<around*|\<langle\>|C,X|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|\<cal-A\>X=b>|<cell|>>|<row|<cell|>|<cell|X\<succeq\>0>|<cell|>>>>
  </eqnarray*>

  The solver implements the dual interior point method
  <cite|benson2008algorithm> and incorporates several important extensions
  that enhance the performance of the dual method. HDSDP has been tested over
  massive datasets and proves a robust and e\^cient SDP solver
  <cite|gao2022hdsdp>.

  <subsection|Release note>

  <\itemize>
    <item>Sep 2022

    Initial release of <samp|HDSDP> 1.0

    <item>June 2023

    The whole repo rewritten based on new design structure

    <item>March 2024

    Added a dedicated LP solver

    <item>July 2024

    Added a primal SDP solver and hybrid primal-primal-dual LP solver

    <item>Oct 2024

    Added new linear system support
  </itemize>

  <new-page*><section|Installation and requirement>

  HDSDP is written in ANSI C and is freely available. Since version 2.0,
  HDSDP is maintained as personal repository

  <center|<samp|<center|><hlink|https://github.com/COPT-Public/HDSDP|>>>

  and all the source files for HDSDP are directly available.\ 

  <subsection|Package dependence>

  HDSDP is a standalone subroutine library, and the only external dependence
  is <samp|LAPACK> and <samp|BLAS>, which are easily accessible for most
  platforms. Alternatively, HDSDP can be linked with <samp|Intel MKL> and
  leverage its multi-thread computing features.

  <paragraph|Linking with LAPACK and BLAS.>Different distributions of
  <samp|BLAS> and <samp|LAPACK> sometimes have different ways of naming. For
  example, <with|font-family|tt|ddot> can appear as
  <with|font-family|tt|ddot>, <with|font-family|tt|ddot_>,
  <with|font-family|tt|DDOT>, or <with|font-family|tt|DDOT_>. Users can
  enable the following macros in <with|font-family|tt|external/lapack_names>
  to ensure that the symbols are correctly detected:

  <\itemize>
    <item><with|font-family|tt|#define UNDERBLAS>

    For <with|font-family|tt|ddot_>

    <item><with|font-family|tt|#define CAPBLAS>

    For <with|font-family|tt|DDOT>

    <item><with|font-family|tt|#define UNDERCAPBLAS>

    For <with|font-family|tt|DDOT_>
  </itemize>

  <subsection|Building the HDSDP library and SDPA/MPS solver>

  HDSDP can be used in two ways: as a subroutine library or as an SDPA
  solver. The subroutine library allows users to call HDSDP within their
  customized applications, while the SDPA solver mainly serves for
  benchmarking purpose. Both the library and the SDPA solver can be built
  using the <samp|CMAKE> build system.\ 

  The <samp|CMakeLists.txt> file located in the repository should be modified
  based on the user's local environment.

  <\itemize>
    <item>Line 2. Target architecture

    If a MacOS user is using HDSDP linked with MKL, uncomment the following
    line

    <\code>
      # set(CMAKE_OSX_ARCHITECTURES x86_64)
    </code>

    <item>Line 30. MKL path

    If a user wants to use HDSDP linked with MKL, set the path to the MKL
    libraries according to the instructions in <samp|CMakeLists.txt>
  </itemize>

  By default HDSDP is compiled without MKL, and the only requirements are
  <samp|LAPACK> and <samp|BLAS> routines. Then the library can be built with

  <\listing>
    <\verbatim>
      mkdir build
    </verbatim>

    <\verbatim>
      cd build
    </verbatim>

    <\verbatim>
      cmake .. \ # Add "-DLINSYS_PARDISO=ON" if MKL is linked
    </verbatim>

    <\verbatim>
      make
    </verbatim>
  </listing>

  and the build system will generate three files

  <\itemize>
    <item>Static library <samp|<with|font-family|tt|hdsdp>>

    <item>Dynamic library <samp|<with|font-family|tt|hdsdp>>

    <item>Executable binary <samp|<with|font-family|tt|sdpasolve>>
  </itemize>

  <subsection|Using the executable binary>

  The <with|font-family|tt|sdpasolve> executable is used for benchmarking
  purpose and accepts two types of input

  <\listing>
    <\verbatim>
      ./sdpasolve data.dat-s # For SDP problems
    </verbatim>

    <\verbatim>
      ./sdpasolve data.mps \ \ # For standard-form LP problems
    </verbatim>
  </listing>

  <new-page*><section|Data and solver interface>

  HDSDP is designed as a stand-alone optimization solver and provides a
  self-contained interface. After the data is input, the solver will initiate
  the subsequent solution phases automatically until the solution procedure
  ends or fails. With the compiled library, the user can call the solver by
  including the header <with|font-family|tt|interface/hdsdp.h> and the data
  interface <with|font-family|tt|interface/hdsdp_user_data.h>. See
  <with|font-family|tt|test/test_file_io.c> for more examples of calling
  <with|font-family|tt|HDSDP>.

  The header <with|font-family|tt|hdsdp.h> defines several utility macros,
  including return code, status code, parameters, constants and the solver
  interface.

  The header <with|font-family|tt|hdsdp_user_data.h> defines the user input
  data interface of <with|font-family|tt|HDSDP>.

  <subsection|Notations>

  <with|font-family|tt|HDSDP> internally solves the following mixed SDP-LP
  problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|X,x>>|<cell|<big|sum><rsub|j=1><rsup|p><around*|\<langle\>|C<rsub|j>,X<rsub|j>|\<rangle\>>+<around*|\<langle\>|c,x|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<big|sum><rsub|j=1><rsup|p><around*|\<langle\>|A<rsub|i
    j>,X<rsub|j>|\<rangle\>>+<around*|\<langle\>|a<rsub|i>,x|\<rangle\>>=b<rsub|i>,>|<cell|i=1,\<ldots\>,m>>|<row|<cell|>|<cell|X<rsub|j>\<in\>\<bbb-S\><rsub|+><rsup|n<rsub|j>>,>|<cell|j=1,\<ldots\>,p>>|<row|<cell|>|<cell|x\<in\>\<bbb-R\><rsup|n<rsub|l>><rsub|+>.>|<cell|>>>>
  </eqnarray*>

  We define the following notations

  <\itemize>
    <item><math|m> is the number of constraints

    <item><math|p> is the number of SDP cones

    <item><math|n<rsub|j>> is the cone dimension of the <math|j>-th SDP cone

    <item><math|n<rsub|s>=<big|sum><rsub|j=1><rsup|p>n<rsub|j>> is the total
    SDP cone dimension

    <item><math|n<rsub|l>> is the LP cone dimension.
  </itemize>

  <subsection|Return code>

  HDSDP subroutines use return code to handle error and exceptions.
  Subroutines with return code <with|font-family|tt|hdsdp_retcode> returns
  one of the codes below

  <\itemize>
    <item><with|font-family|tt|HDSDP_RETCODE_OK>

    The subroutine successfully exits

    <item><with|font-family|tt|HDSDP_RETCODE_FAILED>

    The subroutine exits due to exception or error

    <item><with|font-family|tt|HDSDP_RETCODE_MEMORY>

    The subroutine exits due to memory allocation issues
  </itemize>

  Users can use the <with|font-family|tt|HDSDP_CALL(func)> macro from
  <with|font-family|tt|interface/hdsdp_utils.h> when calling HDSDP routines,
  and the macro will automatically jump to <with|font-family|tt|exit_cleanup>
  label when a non-OK status code is detected: for example

  <\cpp-code>
    HDSDP_CALL(HDSDPOptimize(hdsdp, 1));

    exit_cleanup:

    \ \ \ \ HDSDPDestroy(&hdsdp);
  </cpp-code>

  <subsection|Data format>

  The data interface of HDSDP is defined in
  <with|font-family|tt|interface/hdsdp_user_data.h> and allows users to input
  <math|<around*|{|C<rsub|j>|}>,<around*|{|A<rsub|i
  j>|}>,<around*|{|a<rsub|i>|}>,c> using sparse matrix input.

  <subsubsection|SDP cone data>

  Each SDP cone in HDSDP is determined by <math|m> constraint matrices and an
  objective matrix\ 

  <\equation*>
    <around*|{|A<rsub|1>,\<ldots\>,A<rsub|m>,C|}>,
  </equation*>

  where <math|<around*|{|A<rsub|i>|}>> and <math|C> are all symmetric
  matrices. Given a symmetric matrix <math|A\<in\>\<bbb-S\><rsup|n>>:\ 

  <\equation*>
    A=<matrix|<tformat|<table|<row|<cell|\<downarrow\><with|color|red|a<rsub|1
    1>>>|<cell|>|<cell|>|<cell|>>|<row|<cell|\<downarrow\>a<rsub|2
    1>>|<cell|\<downarrow\>a<rsub|2 2>>|<cell|>|<cell|>>|<row|<cell|\<vdots\>>|<cell|\<vdots\>>|<cell|\<ddots\>>|<cell|>>|<row|<cell|<with|color|red|a<rsub|n
    1>>>|<cell|a<rsub|n 2>>|<cell|\<cdots\>>|<cell|a<rsub|n n>>>>>>,
  </equation*>

  we define its vectorization <math|vec<around*|(|A|)>:\<bbb-S\><rsup|n>\<rightarrow\>\<bbb-R\><rsup|<tfrac|n\<times\><around*|(|n+1|)>|2>>>
  as follows

  <\equation*>
    vec<around*|(|A|)>\<assign\><around*|(|<with|color|red|a<rsub|11>>,a<rsub|21>,\<cdots\><with|color|red|a<rsub|n1>>,\<ldots\>,a<rsub|22>,a<rsub|23>\<ldots\>,a<rsub|n
    n>|)><rsup|\<top\>>.
  </equation*>

  The vectorized version of <math|A> only contains the lower triangular part
  of <math|A>. Then, for each SDP block we can concatenate the coefficient
  matrices by

  <\equation*>
    V<rsub|SDP>\<assign\><matrix|<tformat|<table|<row|<cell|\|>|<cell|\|>|<cell|\|>|<cell|>|<cell|\|>>|<row|<cell|vec<around*|(|C|)>>|<cell|vec<around*|(|A<rsub|1>|)>>|<cell|vec<around*|(|A<rsub|2>|)>>|<cell|\<cdots\>>|<cell|vec<around*|(|A<rsub|m>|)>>>|<row|<cell|\|>|<cell|\|>|<cell|\|>|<cell|>|<cell|\|>>>>>\<in\>\<bbb-R\><rsup|<tfrac|n\<times\><around*|(|n+1|)>|2>\<times\><around*|(|m+1|)>><rsup|>.
  </equation*>

  Then the matrix <math|V<rsub|SDP>> is compressed using standard compressed
  sparse column (<with|font-family|tt|CSC>) format.

  <\example>
    Consider the following 2D SDP

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|X>>|<cell|<around*|\<langle\>|<with|color|red|<matrix|<tformat|<table|<row|<cell|6>|<cell|2>>|<row|<cell|2>|<cell|3>>>>>>,X|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<around*|\<langle\>|<with|color|blue|<matrix|<tformat|<table|<row|<cell|1>|<cell|0>>|<row|<cell|0>|<cell|1>>>>>>,X|\<rangle\>>=1>|<cell|>>|<row|<cell|>|<cell|<around*|\<langle\>|<with|color|orange|<matrix|<tformat|<table|<row|<cell|0>|<cell|1>>|<row|<cell|1>|<cell|0>>>>>>,X|\<rangle\>>=2>|<cell|>>|<row|<cell|>|<cell|X\<succeq\>0.>|<cell|>>>>
    </eqnarray*>

    The problem has 2 constraints and an SDP cone. Then the SDP cone can be
    represented by

    <\equation*>
      V<rsub|SDP>=<matrix|<tformat|<cwith|1|-1|2|2|color|blue>|<cwith|1|-1|3|3|color|orange>|<cwith|1|3|1|1|color|red>|<table|<row|<cell|6>|<cell|1>|<cell|0>>|<row|<cell|2>|<cell|0>|<cell|2>>|<row|<cell|3>|<cell|1>|<cell|0>>>>>
    </equation*>

    and the cone can be represented by the following sparse matrix format.\ 

    <\listing>
      <\verbatim>
        Vp[4] = {0, 3, 5, 6};
      </verbatim>

      <\verbatim>
        Vi[6] = {0, 1, 2, 0, 2, 1};
      </verbatim>

      <\verbatim>
        Vx[6] = {6, 2, 3, 1, 1, 2};
      </verbatim>
    </listing>
  </example>

  <subsubsection|LP cone data>

  LP cone data is similar to SDP, where we concatenate
  <math|<around*|{|a<rsub|i>|}>> by

  <\equation*>
    V<rsub|LP>=<matrix|<tformat|<table|<row|<cell|c>|<cell|a<rsub|1>>|<cell|a<rsub|2>>|<cell|\<cdots\>>|<cell|a<rsub|m>>>>>>\<in\>\<bbb-R\><rsup|n\<times\><around*|(|m+1|)>>
  </equation*>

  and <math|V<rsub|LP>> is also compressed by <with|font-family|tt|CSC>
  format.

  <paragraph|Scalar bound cone.>HDSDP also provides a scalar bound cone for
  constraints on the dual variables <math|l\<cdummy\>e\<leq\>y\<leq\>u\<cdummy\>e>.
  HDSDP will detect in the presolve phase whether an LP cone is a scalar
  bound cone, and will switch to bound cone when it is detected.

  <subsection|Data interface>

  <subsubsection|Create user data structure>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern hdsdp_retcode HUserDataCreate( user_data **pHdata );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Create an HDSDP user data pointer

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|pHdata>

    Type: <with|font-family|tt|user_data **>

    Address of the pointer to which the user data structure is to be
    allocated
  </itemize>

  <with|font-series|bold|Return>

  Status of execution

  <subsubsection|Set cone data>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HUserDataSetConeData( user_data *Hdata,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ cone_type
    cone,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ int
    nRow,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ int
    nCol,

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ int
    *coneMatBeg,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ int
    *coneMatIdx,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ double
    *coneMatElem );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Set cone data in an existing user data structure

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|Hdata>

    Type: <with|font-family|tt|user_data *>

    Pointer of the user data structure

    <item><with|font-family|tt|cone>

    Type: <with|font-family|tt|cone_type *>

    Type of the cone specified by the user. Possible choices are\ 

    <\itemize>
      <item><with|font-family|tt|HDSDP_CONETYPE_LP>

      LP cone

      <item><with|font-family|tt|HDSDP_CONETYPE_DENSE_SDP>

      SDP cone

      <item><with|font-family|tt|HDSDP_CONETYPE_SPARSE_SDP>

      SDP cone where most of the <math|<around*|{|A<rsub|i>|}>> are all-zero
      matrices
    </itemize>

    <item><with|font-family|tt|nRow>

    Type: <with|font-family|tt|int>

    Number of constraints of the cone. Should be the same for all the cones

    <item><with|font-family|tt|nCol>

    Type: <with|font-family|tt|int>

    Dimension of the cone, <math|n<rsub|l>> or <math|n<rsub|j>>

    <item><with|font-family|tt|coneMatBeg>

    Type: <with|font-family|tt|int *>

    Column pointer for the <with|font-family|tt|CSC> sparse matrix

    <item><with|font-family|tt|coneMatIdx>

    Type: <with|font-family|tt|int *>

    Row indices for the <with|font-family|tt|CSC> sparse matrix

    <item><with|font-family|tt|coneMatElem>

    Type: <with|font-family|tt|double *>

    Elements for the <with|font-family|tt|CSC> sparse matrix
  </itemize>

  <with|font-series|bold|Return>

  No return

  <subsubsection|Determine cone type>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern cone_type HUserDataChooseCone( user_data *Hdata );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Automatically detect and adjust the cone type based on problem data. The
  routine can determine <with|font-series|bold|1)> whether a dense SDP cone
  should instead be specified as a sparse SDP cone <with|font-series|bold|2)>
  whether an LP cone is a bound cone. This routine will be automatically
  invoked within the HDSDP solver routine

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|Hdata>

    Type: <with|font-family|tt|user_data *>

    Pointer of the user data structure with cone data available
  </itemize>

  <with|font-series|bold|Return>

  Type of the cone detected by the solver

  <subsubsection|Clear user data>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HUserDataClear( user_data *Hdata );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Clear the user data

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|Hdata>

    Type: <with|font-family|tt|user_data *>

    Pointer of the user data
  </itemize>

  <with|font-series|bold|Return>

  No return

  <subsubsection|Destroy user data structure>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HUserDataDestroy( user_data **pHdata );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Clear and release all the memory of user data

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|pHdata>

    Type: <with|font-family|tt|user_data **>

    Address of the pointer to which the user data structure is allocated
  </itemize>

  <with|font-series|bold|Return>

  No return

  <subsection|Solver interface>

  <subsubsection|Create the solver structure>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern hdsdp_retcode HDSDPCreate( hdsdp **pHSolver );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Create an HDSDP solver data structure

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|pHSolver>

    Type: <with|font-family|tt|hdsdp **>

    Address of the pointer to which the solver is to be allocated
  </itemize>

  <with|font-series|bold|Return>

  Status of execution

  <subsubsection|Solver initialization>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern hdsdp_retcode HDSDPInit( hdsdp *HSolver, int nRows, int nCones );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Initialize solver and allocate solver memory

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|nRows>

    Type: <with|font-family|tt|int>

    Number of constraints of the problem

    <item><with|font-family|tt|nCones>

    Type: <with|font-family|tt|int nCones>

    Number of cones of the problem
  </itemize>

  <with|font-series|bold|Return>

  Status of execution

  <subsubsection|Set conic data>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern hdsdp_retcode HDSDPSetCone( hdsdp *HSolver, int iCone, void
    *userCone );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Set cone data in the solver

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|iCone>

    Type: <with|font-family|tt|int>

    Index of the cone. Range between 0 and <with|font-family|tt|nCones>
    specified when calling <with|font-family|tt|HDSDPInit>

    <item><with|font-family|tt|userCone>

    Type: <with|font-family|tt|void *>

    Pointer to user cone data
  </itemize>

  <with|font-series|bold|Return>

  Status of execution

  <subsubsection|Set dual objective>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HDSDPSetDualObjective( hdsdp *HSolver, double *dObj );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Set the primal constraint right-hand-side and the dual objective <math|b>

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|dObj>

    Type: <with|font-family|tt|double *>

    Dual objective. A vector of length <math|m>
  </itemize>

  <with|font-series|bold|Return>

  No return

  <subsubsection|Set dual starting point>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HDSDPSetDualStart( hdsdp *HSolver, double *dStart );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Set the starting point for the dual problem <math|y<rsub|0>>.

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|dStart>

    Type: <with|font-family|tt|double *>

    Dual variable starting point. A vector of length <math|m>
  </itemize>

  <with|font-series|bold|Return>

  No return

  <subsubsection|Set integer parameter>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HDSDPSetIntParam( hdsdp *HSolver,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ int intParam,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ int
    intParamVal );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Set an integer parameter (see the next section for the available
  parameters)

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|intParam>

    Type: <with|font-family|tt|int>

    Integer parameter name

    <item><with|font-family|tt|intParamVal>

    Type: <with|font-family|tt|int>

    Value of the integer parameter
  </itemize>

  <with|font-series|bold|Return>

  No return

  <subsubsection|Set double parameter>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HDSDPSetDblParam( hdsdp *HSolver,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ int dblParam,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ double
    dblParamVal );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Set a double parameter (see the next section for the available parameters)

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|dblParam>

    Type: <with|font-family|tt|int>

    Double parameter name

    <item><with|font-family|tt|dblParamVal>

    Type: <with|font-family|tt|double>

    Value of the double parameter
  </itemize>

  <with|font-series|bold|Return>

  No return

  <subsubsection|Optimization>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern hdsdp_retcode HDSDPOptimize( hdsdp *HSolver, int dOptOnly );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Invoke HDSDP optimization routine

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|dOptOnly>

    Type: <with|font-family|tt|int>

    Whether HDSDP should only optimize the dual problem. If
    <with|font-family|tt|dOptOnly> is nonzero, then the algorithm will
    perform Phase 1 (infeasible-start embedding phase) and will not go into
    primal solution extraction.
  </itemize>

  <with|font-series|bold|Return>

  Status of execution

  <subsubsection|Get dual solution>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern hdsdp_retcode HDSDPGetRowDual( hdsdp *HSolver,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ double
    *pObjVal,

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ double
    *dObjVal,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ double
    *dualVal );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Extract dual optimal value after optimization finishes

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|pObjVal>

    Type: <with|font-family|tt|double *>

    Optimal primal objective <em|value>. Can be <with|font-family|tt|NULL> if
    objective value is not needed

    <item><with|font-family|tt|dObjVal>

    Type: <with|font-family|tt|double *>

    Optimal dual objective <em|value>. Can be <with|font-family|tt|NULL> if
    objective value is not needed

    <item><with|font-family|tt|dualVal>

    Type: <with|font-family|tt|double *>

    Optimal dual objective solution <math|y<rsup|\<star\>>>. Can be
    <with|font-family|tt|NULL> if solution is not needed
  </itemize>

  <with|font-series|bold|Return>

  Status of execution

  <subsubsection|Get cone solution>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HDSDPGetConeValues( hdsdp *HSolver,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ int
    iCone,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ double
    *conePrimal,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ double
    *coneDual,\ 

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ double
    *coneAuxi );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Extract the cone primal optimal solution. Auxiliary memory is required for
  extraction.

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|iCone>

    Type: <with|font-family|tt|int>

    Index of the cone

    <item><with|font-family|tt|conePrimal>

    Type: <with|font-family|tt|double *>

    Optimal primal solution <math|X<rsup|\<star\>>>

    <item><with|font-family|tt|coneDual>

    Type: <with|font-family|tt|double *>

    Optimal dual slack <math|S<rsup|\<star\>>>

    <item><with|font-family|tt|coneAuxi>

    Type: <with|font-family|tt|double *>

    Auxiliary memory for extraction of size <math|n<rsub|j>> for SDP) or
    <math|n<rsub|l>> (for LP).
  </itemize>

  <with|font-series|bold|Return>

  No return

  <subsubsection|Check solution and compute solution metrics>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern hdsdp_retcode HDSDPCheckSolution( hdsdp *HSolver, double
    diErrors[6] );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Check the solution quality and compute the DIMACS errors

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver

    <item><with|font-family|tt|diErrors>

    Type: <with|font-family|tt|double *>

    The six DIMACS errors indicating the quality of the solution
  </itemize>

  <with|font-series|bold|Return>

  Status of execution

  <subsubsection|Release solver memory>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HDSDPClear( hdsdp *HSolver );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Release all the internal memory allocated by HDSDP

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|HSolver>

    Type: <with|font-family|tt|hdsdp *>

    Pointer of the solver
  </itemize>

  <with|font-series|bold|Return>

  No return

  <subsubsection|Destroy solver structure>

  <with|font-series|bold|Function>

  <\cpp-code>
    extern void HDSDPDestroy( hdsdp **pHSolver );
  </cpp-code>

  <with|font-series|bold|Explanation>

  Clear and release all the memory of HDSDP.

  <with|font-series|bold|Argument>

  <\itemize>
    <item><with|font-family|tt|pHSolver>

    Type: <with|font-family|tt|hdsdp **>

    Pointer of the solver
  </itemize>

  <with|font-series|bold|Return>

  No return

  <new-page*><section|Parameters>

  <subsection|Integer parameters>

  <subsubsection|INT_PARAM_CORRECTORA>

  <with|font-series|bold|Explanation>: Number of corrector steps in HDSDP
  infeasible embedding phase. More corrector steps reduce the number of
  iterations but can increase time spent per iteration

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: Adjusted by problem feature

  <subsubsection|INT_PARAM_CORRECTORB>

  <with|font-series|bold|Explanation>: Number of corrector steps in HDSDP
  feasible dual potential reduction phase. More corrector steps reduce the
  number of iterations but can increase time spent per iteration

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: Adjusted by problem feature

  <subsubsection|INT_PARAM_MAXITER>

  <with|font-series|bold|Explanation>: Maximum number of potential reduction
  iterations

  <with|font-series|bold|Range>: <math|\<geq\>1>

  <with|font-series|bold|Default>: 500

  <subsubsection|INT_PARAM_PSDP>

  <with|font-series|bold|Explanation>: <em|Experiment feature>. Whether
  primal IPM solver is turned on when dual SDP is close to convergence

  <with|font-series|bold|Range>: <math|<around*|{|0,1|}>>

  <with|font-series|bold|Default>: 0

  <subsubsection|INT_PARAM_PRELEVEL>

  <with|font-series|bold|Explanation>: Presolving level for problem data and
  KKT solver. In general no benefit using low level presolving

  <with|font-series|bold|Range>: <math|<around*|{|0,1,2|}>>

  <with|font-series|bold|Default>: 2

  <subsubsection|INT_PARAM_THREADS>

  <with|font-series|bold|Explanation>: Number of threads for linear algebra
  subroutines. Only effective when Intel MKL is linked to HDSDP.

  <with|font-series|bold|Range>: <math|\<geq\>1>

  <with|font-series|bold|Default>: 12

  <subsection|Double parameters>

  <subsubsection|DBL_PARAM_RELFEASTOL>

  <with|font-series|bold|Explanation>: Relative feasibility tolerance for
  <math|max<around*|{|<frac|<around*|\<\|\|\>|\<cal-A\>X-b|\<\|\|\>>|1+<around*|\<\|\|\>|b|\<\|\|\>><rsub|1>>,<frac|<around*|\<\|\|\>|\<cal-A\><rsup|\<ast\>>y+S-C|\<\|\|\>><rsub|F>|1+<around*|\<\|\|\>|C|\<\|\|\>><rsub|sum>>|}>>

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: <math|10<rsup|-8>>, also adjusted by
  problem feature

  <subsubsection|DBL_PARAM_RELOPTTOL>

  <with|font-series|bold|Explanation>: Relative optimality tolerance for
  <math|<frac|<around*|\||<around*|\<langle\>|C,X|\<rangle\>>-<around*|\<langle\>|b,y|\<rangle\>>|\|>|<around*|\||<around*|\<langle\>|C,X|\<rangle\>>|\|>+<around*|\||<around*|\<langle\>|b,y|\<rangle\>>|\|>+1>>

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: <math|10<rsup|-8>>, also adjusted by
  problem feature

  <subsubsection|DBL_PARAM_ABSFEASTOL>

  <with|font-series|bold|Explanation>: Absolute feasibility tolerance for
  <math|max<around*|{|<around*|\<\|\|\>|\<cal-A\>X-b|\<\|\|\>>,<around*|\<\|\|\>|\<cal-A\><rsup|\<ast\>>y+S-C|\<\|\|\>><rsub|F>|}>>

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: <math|10<rsup|-8>>, also adjusted by
  problem feature

  <subsubsection|DBL_PARAM_ABSOPTTOL>

  <with|font-series|bold|Explanation>: Absolute optimality tolerance for
  <math|<around*|\||<around*|\<langle\>|C,X|\<rangle\>>-<around*|\<langle\>|b,y|\<rangle\>>|\|>>

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: <math|10<rsup|-8>>, also adjusted by
  problem feature

  <subsubsection|DBL_PARAM_TIMELIMIT>

  <with|font-series|bold|Explanation>: Solver running time limit

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: 3600.0

  <subsubsection|DBL_PARAM_POTRHOVAL>

  <with|font-series|bold|Explanation>: Potential function parameter\ 

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: Adjusted by problem feature

  <subsubsection|DBL_PARAM_HSDGAMMA>

  <with|font-series|bold|Explanation>: Infeasibility elimination
  aggressiveness in embedding phase

  <with|font-series|bold|Range>: <math|<around*|[|0,1|]>>

  <with|font-series|bold|Default>: Adjusted by problem feature

  <subsubsection|DBL_PARAM_DUALBOX_LOW>

  <with|font-series|bold|Explanation>: Scalar lower bound on the dual
  variables <math|y\<geq\>l\<cdummy\>e>

  <with|font-series|bold|Range>: <math|\<bbb-R\>>

  <with|font-series|bold|Default>: Adjusted by problem feature

  <subsubsection|DBL_PARAM_DUALBOX_UP>

  <with|font-series|bold|Explanation>: Scalar upper bound on the dual
  variables <math|y\<leq\>u\<cdummy\>e>

  <with|font-series|bold|Range>: <math|\<bbb-R\>>

  <with|font-series|bold|Default>: Adjusted by problem feature

  <subsubsection|DBL_PARAM_BARMUSTART>

  <with|font-series|bold|Explanation>: Starting barrier parameter

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: <math|10<rsup|5>>

  <subsubsection|DBL_PARAM_DUALSTART>

  <with|font-series|bold|Explanation>: Starting dual slack size:
  <math|\<theta\>> in <math|S<rsub|0>=C-\<cal-A\><rsup|\<ast\>>y<rsub|0>+\<theta\>
  I>

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: <math|10<rsup|5>>

  <subsubsection|DBL_PARAM_POBJSTART>

  <with|font-series|bold|Explanation>: Initial primal objective guess

  <with|font-series|bold|Range>: <math|\<bbb-R\>>

  <with|font-series|bold|Default>: <math|10<rsup|10>>

  <subsubsection|DBL_PARAM_TRXESTIMATE>

  <with|font-series|bold|Explanation>: Primal optimal solution trace estimate

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: <math|10<rsup|8>>

  <subsubsection|DBL_PARAM_PRECORDACC>

  <with|font-series|bold|Explanation>: Accuracy for preparing for primal
  solution recovery. When the dual solution is too accurate, extraction of
  primal solution can fail due to significant numerical errors. Therefore,
  this parameter should not be set too small.

  <with|font-series|bold|Range>: <math|\<geq\>0>

  <with|font-series|bold|Default>: <math|10<rsup|-8>>

  <section|Logging and examples>

  HDSDP has a logging system that allows users to keep track of the
  optimization progress.

  <subsection|SDP logging>

  <\center>
    <\frame>
      <\small>
        <\shell-code>
          Filename: ../examples/mcp100.dat-s

          Reading SDPA file in 0.000535 seconds<nbsp>

          \;

          HDSDP: software for semi-definite programming<nbsp>

          \;

          Wenzhi Gao, Dongdong Ge, Yinyu Ye, 2024

          ---------------------------------------------

          Pre-solver starts<nbsp>

          <nbsp> Processing the cones<nbsp>

          <nbsp> <nbsp> M1: 0<nbsp> <nbsp> M2: 100<nbsp> M3: 0<nbsp> <nbsp>
          M4: 0<nbsp> <nbsp> M5: 0 <nbsp> <nbsp>

          <nbsp> Starting KKT analysis<nbsp>

          <nbsp> <nbsp> Using dense Schur complement

          <nbsp> Collecting statistics<nbsp>

          <nbsp> Making adjustments

          <nbsp> <nbsp> Scale cone objective by 1.0e+00<nbsp>

          <nbsp> <nbsp> Scale rhs by 1.0e+00<nbsp>

          <nbsp> <nbsp> Hardware has 1 thread(s)

          Pre-solver ends. Elapsed time: 0.0 seconds<nbsp>

          \;

          Statistics<nbsp>

          <nbsp> Number of rows: 100

          <nbsp> Number of cones: 1

          <nbsp> Number of sparse SDP cones: 0

          ...

          <nbsp> Norm of objective: 2.69e+02

          <nbsp> Norm of SDP data: 1.00e+02

          <nbsp> Norm of RHS: 1.00e+02

          \;

          Parameters

          <nbsp> Maximum iteration <nbsp> : 500

          <nbsp> Infeasible corrector: 4

          <nbsp> Feasible corrector<nbsp> : 0

          \ \ ...

          <nbsp> Dual box low<nbsp> <nbsp> <nbsp> <nbsp> : -1.0e+06

          <nbsp> Dual box up <nbsp> <nbsp> <nbsp> <nbsp> : 1.0e+06

          <nbsp> Starting primal <nbsp> <nbsp> : 1.0e+08

          \;

          This is a trace-implied SDP problem

          Optimizing over 1 thread(s)<nbsp>

          Initialize with dual residual 1.0e+05

          HDSDP starts. Using infeasible dual method<nbsp>

          \;

          <nbsp> <nbsp> nIter <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> pObj
          <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> dObj<nbsp> <nbsp> <nbsp>
          dInf<nbsp> <nbsp> <nbsp> <nbsp> Mu <nbsp> Step <nbsp> <nbsp> \|P\|
          <nbsp> T [D]<nbsp> \ \ \ \ \ \ \ 

          <nbsp> <nbsp> <nbsp> <nbsp> 1<nbsp> +7.33331918e+07<nbsp>
          -7.35913354e+07<nbsp> 0.00e+00<nbsp> 3.56e+04 <nbsp> 0.76<nbsp>
          7.2e+00 <nbsp> 0.0<nbsp>

          \;

          Infeasible method finds a dual feasible solution<nbsp>

          HDSDP re-starts. Using feasible dual method<nbsp>

          \;

          <nbsp> <nbsp> nIter <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> pObj
          <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> dObj<nbsp> <nbsp> <nbsp>
          pInf<nbsp> <nbsp> <nbsp> <nbsp> Mu <nbsp> Step <nbsp> <nbsp> \|P\|
          <nbsp> T [P]<nbsp>

          <nbsp> <nbsp> <nbsp> <nbsp> 1<nbsp> +7.33331918e+07<nbsp>
          -4.81595948e+07<nbsp> 5.23e-02<nbsp> 4.26e+05 <nbsp> 1.00<nbsp>
          1.0e+01 <nbsp> 0.0<nbsp>

          ...

          <nbsp><nbsp> <nbsp> <nbsp> 34<nbsp> -2.26157351e+02<nbsp>
          -2.26157351e+02<nbsp> 0.00e+00<nbsp> 5.28e-12 <nbsp> 0.12<nbsp>
          1.2e+01 <nbsp> 0.1<nbsp>

          \;

          Optimization time: 0.1 seconds

          DIMACS error metric:

          <nbsp> <nbsp> 1.10e-13 0.00e+00 0.00e+00 0.00e+00 6.50e-09
          6.50e-09<nbsp>

          \;

          SDP Status: Primal dual optimal<nbsp>

          <nbsp> pObj -2.2615734854e+02

          <nbsp> dObj -2.2615735148e+02

          PD Gap +2.9464339377e-06

          <nbsp> Time 0.1 seconds
        </shell-code>
      </small>
    </frame>
  </center>

  <subsection|LP logging>

  <\center>
    <frame|<\small>
      <\shell-code>
        Reading specialized standard form mps ../examples/afiro.mps<nbsp>

        Optimizing an LP of 23 variables and 15 constraints

        Data statistics: \|A\| = 1.06e+02 \|b\| = 9.82e+03 \|c\| = 2.08e+00
        Nnz = 49<nbsp>

        Using Hybrid Primal-Primal-Dual solver<nbsp>

        \;

        <nbsp> <nbsp> nIter <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> pObj
        <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> <nbsp> dObj<nbsp> <nbsp> <nbsp>
        pInf<nbsp> <nbsp> <nbsp> dInf<nbsp> <nbsp> <nbsp> <nbsp> Mu<nbsp>
        <nbsp> P/D Step <nbsp> T [PD]<nbsp>

        <nbsp> <nbsp> <nbsp> <nbsp> 0<nbsp> -2.76027927e+02<nbsp>
        -1.04931952e+02<nbsp> 8.59e-01<nbsp> 1.91e+00 <nbsp> 0.00<nbsp> 0.00
        <nbsp> 0.0<nbsp>

        <nbsp> <nbsp> <nbsp> <nbsp> 1<nbsp> -2.54976427e+02<nbsp>
        -4.55135868e+02<nbsp> 1.83e-02<nbsp> 3.93e-01 <nbsp> 0.98<nbsp> 0.79
        <nbsp> 0.0<nbsp>

        <nbsp> <nbsp> <nbsp> <nbsp> 2<nbsp> -4.33824047e+02<nbsp>
        -5.31380377e+02<nbsp> 4.30e-03<nbsp> 4.39e-03 <nbsp> 0.76<nbsp> 0.99
        <nbsp> 0.0<nbsp>

        <nbsp> <nbsp> <nbsp> <nbsp> 3<nbsp> -4.63764901e+02<nbsp>
        -4.83113579e+02<nbsp> 5.05e-04<nbsp> 9.48e-04 <nbsp> 0.88<nbsp> 0.78
        <nbsp> 0.0<nbsp>

        <nbsp> <nbsp> <nbsp> <nbsp> 4<nbsp> -4.64643942e+02<nbsp>
        -4.65030078e+02<nbsp> 3.98e-15<nbsp> 1.41e-05 <nbsp> 1.00<nbsp> 0.99
        <nbsp> 0.0<nbsp>

        <nbsp> <nbsp> <nbsp> <nbsp> 5<nbsp> -4.64752596e+02<nbsp>
        -4.64754528e+02<nbsp> 4.93e-16<nbsp> 7.06e-08 <nbsp> 0.99<nbsp> 0.99
        <nbsp> 0.0<nbsp>

        <nbsp> <nbsp> <nbsp> <nbsp> 6<nbsp> -4.64753140e+02<nbsp>
        -4.64753150e+02<nbsp> 2.70e-15<nbsp> 3.53e-10 <nbsp> 1.00<nbsp> 1.00
        <nbsp> 0.0<nbsp>

        <nbsp> <nbsp> <nbsp> <nbsp> 7<nbsp> -4.64753143e+02<nbsp>
        -4.64753143e+02<nbsp> 3.22e-15<nbsp> 1.77e-12 <nbsp> 1.00<nbsp> 1.00
        <nbsp> 0.0<nbsp>

        \;

        LP Status: Primal dual optimal<nbsp>

        \;

        LP Solution statistic<nbsp>

        pObj: -4.648e+02 <nbsp> dObj: -4.648e+02<nbsp>

        Abs. pInf:<nbsp> 1.152e-11<nbsp> <nbsp> Rel. pInf:<nbsp>
        3.221e-15<nbsp>

        Abs. dInf:<nbsp> 4.266e-12<nbsp> <nbsp> Rel. dInf:<nbsp>
        1.765e-12<nbsp>

        Abs. Gap :<nbsp> 4.830e-08<nbsp> <nbsp> Rel. gap :<nbsp> 5.191e-11

        \;

        Linear solver statistic<nbsp>

        Num. Factor: <nbsp> 7 (+1)<nbsp> <nbsp> Factor Time:<nbsp> 0.00<nbsp>

        Num. Solves:<nbsp> 15 (+1)<nbsp> <nbsp> Solve<nbsp> Time:<nbsp>
        0.00<nbsp>

        \;

        Elapsed Time:<nbsp> 0.000 seconds<nbsp>
      </shell-code>
    </small>>
  </center>

  <\bibliography*|bib|tm-plain|../../gwz/hdsdp/hdsdp/doc/sdpref|References>
    <\bib-list|2>
      <bibitem*|1><label|bib-benson2008algorithm>Steven<nbsp>J
      Benson<localize| and >Yinyu Ye. <newblock>Algorithm 875:
      dsdp5\Vsoftware for semidefinite programming.
      <newblock><with|font-shape|italic|ACM Transactions on Mathematical
      Software (TOMS)>, 34(3):1\U20, 2008.<newblock>

      <bibitem*|2><label|bib-gao2022hdsdp>Wenzhi Gao, Dongdong Ge<localize|,
      and >Yinyu Ye. <newblock>Hdsdp: software for semidefinite programming.
      <newblock><with|font-shape|italic|ArXiv preprint arXiv:2207.13862>,
      2022.<newblock>
    </bib-list>
  </bibliography*>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|3>>
    <associate|auto-10|<tuple|3.2|6>>
    <associate|auto-11|<tuple|3.3|7>>
    <associate|auto-12|<tuple|3.3.1|7>>
    <associate|auto-13|<tuple|3.3.2|8>>
    <associate|auto-14|<tuple|3.3.2.1|8>>
    <associate|auto-15|<tuple|3.4|8>>
    <associate|auto-16|<tuple|3.4.1|8>>
    <associate|auto-17|<tuple|3.4.2|8>>
    <associate|auto-18|<tuple|3.4.3|9>>
    <associate|auto-19|<tuple|3.4.4|10>>
    <associate|auto-2|<tuple|1.1|3>>
    <associate|auto-20|<tuple|3.4.5|10>>
    <associate|auto-21|<tuple|3.5|10>>
    <associate|auto-22|<tuple|3.5.1|10>>
    <associate|auto-23|<tuple|3.5.2|11>>
    <associate|auto-24|<tuple|3.5.3|11>>
    <associate|auto-25|<tuple|3.5.4|12>>
    <associate|auto-26|<tuple|3.5.5|12>>
    <associate|auto-27|<tuple|3.5.6|13>>
    <associate|auto-28|<tuple|3.5.7|13>>
    <associate|auto-29|<tuple|3.5.8|14>>
    <associate|auto-3|<tuple|2|4>>
    <associate|auto-30|<tuple|3.5.9|14>>
    <associate|auto-31|<tuple|3.5.10|15>>
    <associate|auto-32|<tuple|3.5.11|16>>
    <associate|auto-33|<tuple|3.5.12|16>>
    <associate|auto-34|<tuple|3.5.13|16>>
    <associate|auto-35|<tuple|4|17>>
    <associate|auto-36|<tuple|4.1|17>>
    <associate|auto-37|<tuple|4.1.1|17>>
    <associate|auto-38|<tuple|4.1.2|17>>
    <associate|auto-39|<tuple|4.1.3|17>>
    <associate|auto-4|<tuple|2.1|4>>
    <associate|auto-40|<tuple|4.1.4|17>>
    <associate|auto-41|<tuple|4.1.5|17>>
    <associate|auto-42|<tuple|4.1.6|17>>
    <associate|auto-43|<tuple|4.2|18>>
    <associate|auto-44|<tuple|4.2.1|18>>
    <associate|auto-45|<tuple|4.2.2|18>>
    <associate|auto-46|<tuple|4.2.3|18>>
    <associate|auto-47|<tuple|4.2.4|18>>
    <associate|auto-48|<tuple|4.2.5|18>>
    <associate|auto-49|<tuple|4.2.6|18>>
    <associate|auto-5|<tuple|1|4>>
    <associate|auto-50|<tuple|4.2.7|18>>
    <associate|auto-51|<tuple|4.2.8|18>>
    <associate|auto-52|<tuple|4.2.9|19>>
    <associate|auto-53|<tuple|4.2.10|19>>
    <associate|auto-54|<tuple|4.2.11|19>>
    <associate|auto-55|<tuple|4.2.12|19>>
    <associate|auto-56|<tuple|4.2.13|19>>
    <associate|auto-57|<tuple|4.2.14|19>>
    <associate|auto-58|<tuple|5|19>>
    <associate|auto-59|<tuple|5.1|20>>
    <associate|auto-6|<tuple|2.2|4>>
    <associate|auto-60|<tuple|5.2|21>>
    <associate|auto-61|<tuple|5.2|23>>
    <associate|auto-62|<tuple|5.3|21>>
    <associate|auto-7|<tuple|2.3|5>>
    <associate|auto-8|<tuple|3|6>>
    <associate|auto-9|<tuple|3.1|6>>
    <associate|bib-benson2008algorithm|<tuple|1|23>>
    <associate|bib-gao2022hdsdp|<tuple|2|23>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      benson2008algorithm

      gao2022hdsdp
    </associate>
    <\associate|toc>
      1<space|2spc>Introduction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Release note
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      2<space|2spc>Installation and requirement
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Package dependence
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|4tab>|Linking with LAPACK and BLAS.
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>Building the HDSDP library
      and SDPA/MPS solver <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|2.3<space|2spc>Using the executable binary
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      3<space|2spc>Data and solver interface
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Notations
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Return code
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|1tab>|3.3<space|2spc>Data format
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|2tab>|3.3.1<space|2spc>SDP cone data
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|2tab>|3.3.2<space|2spc>LP cone data
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|4tab>|Scalar bound cone.
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|3.4<space|2spc>Data interface
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|2tab>|3.4.1<space|2spc>Create user data structure
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|2tab>|3.4.2<space|2spc>Set cone data
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <with|par-left|<quote|2tab>|3.4.3<space|2spc>Determine cone type
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18>>

      <with|par-left|<quote|2tab>|3.4.4<space|2spc>Clear user data
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19>>

      <with|par-left|<quote|2tab>|3.4.5<space|2spc>Destroy user data
      structure <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <with|par-left|<quote|1tab>|3.5<space|2spc>Solver interface
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>>

      <with|par-left|<quote|2tab>|3.5.1<space|2spc>Create the solver
      structure <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22>>

      <with|par-left|<quote|2tab>|3.5.2<space|2spc>Solver initialization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>

      <with|par-left|<quote|2tab>|3.5.3<space|2spc>Set conic data
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>>

      <with|par-left|<quote|2tab>|3.5.4<space|2spc>Set dual objective
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-25>>

      <with|par-left|<quote|2tab>|3.5.5<space|2spc>Set dual starting point
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26>>

      <with|par-left|<quote|2tab>|3.5.6<space|2spc>Set integer parameter
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-27>>

      <with|par-left|<quote|2tab>|3.5.7<space|2spc>Set double parameter
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-28>>

      <with|par-left|<quote|2tab>|3.5.8<space|2spc>Optimization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-29>>

      <with|par-left|<quote|2tab>|3.5.9<space|2spc>Get dual solution
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-30>>

      <with|par-left|<quote|2tab>|3.5.10<space|2spc>Get cone solution
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-31>>

      <with|par-left|<quote|2tab>|3.5.11<space|2spc>Check solution and
      compute solution metrics <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-32>>

      <with|par-left|<quote|2tab>|3.5.12<space|2spc>Release solver memory
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-33>>

      <with|par-left|<quote|2tab>|3.5.13<space|2spc>Destroy solver structure
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-34>>

      4<space|2spc>Parameters <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-35>

      <with|par-left|<quote|1tab>|4.1<space|2spc>Integer parameters
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-36>>

      <with|par-left|<quote|2tab>|4.1.1<space|2spc>INT_PARAM_CORRECTORA
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-37>>

      <with|par-left|<quote|2tab>|4.1.2<space|2spc>INT_PARAM_CORRECTORB
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-38>>

      <with|par-left|<quote|2tab>|4.1.3<space|2spc>INT_PARAM_MAXITER
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-39>>

      <with|par-left|<quote|2tab>|4.1.4<space|2spc>INT_PARAM_PSDP
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-40>>

      <with|par-left|<quote|2tab>|4.1.5<space|2spc>INT_PARAM_PRELEVEL
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-41>>

      <with|par-left|<quote|2tab>|4.1.6<space|2spc>INT_PARAM_THREADS
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-42>>

      <with|par-left|<quote|1tab>|4.2<space|2spc>Double parameters
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-43>>

      <with|par-left|<quote|2tab>|4.2.1<space|2spc>DBL_PARAM_RELFEASTOL
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-44>>

      <with|par-left|<quote|2tab>|4.2.2<space|2spc>DBL_PARAM_RELOPTTOL
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-45>>

      <with|par-left|<quote|2tab>|4.2.3<space|2spc>DBL_PARAM_ABSFEASTOL
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-46>>

      <with|par-left|<quote|2tab>|4.2.4<space|2spc>DBL_PARAM_ABSOPTTOL
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-47>>

      <with|par-left|<quote|2tab>|4.2.5<space|2spc>DBL_PARAM_TIMELIMIT
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-48>>

      <with|par-left|<quote|2tab>|4.2.6<space|2spc>DBL_PARAM_POTRHOVAL
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-49>>

      <with|par-left|<quote|2tab>|4.2.7<space|2spc>DBL_PARAM_HSDGAMMA
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-50>>

      <with|par-left|<quote|2tab>|4.2.8<space|2spc>DBL_PARAM_DUALBOX_LOW
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-51>>

      <with|par-left|<quote|2tab>|4.2.9<space|2spc>DBL_PARAM_DUALBOX_UP
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-52>>

      <with|par-left|<quote|2tab>|4.2.10<space|2spc>DBL_PARAM_BARMUSTART
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-53>>

      <with|par-left|<quote|2tab>|4.2.11<space|2spc>DBL_PARAM_DUALSTART
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-54>>

      <with|par-left|<quote|2tab>|4.2.12<space|2spc>DBL_PARAM_POBJSTART
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-55>>

      <with|par-left|<quote|2tab>|4.2.13<space|2spc>DBL_PARAM_TRXESTIMATE
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-56>>

      <with|par-left|<quote|2tab>|4.2.14<space|2spc>DBL_PARAM_PRECORDACC
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-57>>

      5<space|2spc>Logging and examples <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-58>

      <with|par-left|<quote|1tab>|5.1<space|2spc>SDP logging
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-59>>

      <with|par-left|<quote|1tab>|5.2<space|2spc>LP logging
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-60>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|References>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-61><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>