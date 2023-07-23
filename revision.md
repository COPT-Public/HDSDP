Wenzhi:

Comments from the TOMS Associate Editor and referees have now been prepared for the following manuscript submitted to the ACM Transactions on Mathematical Software:

Identifier: TOMS-2022-0052
Title: Algorithm xxxx: HDSDP: Software for Semi-definite Programming
Author(s): Gao, Wenzhi; Ge, Dongdong; Ye, Yinyu
Type:  2 Algorithm Papers

Comments (included at the foot of this letter, or as PDF attachments) indicate that a major revision to your manuscript will be necessary before it can be considered for publication.

The revised manuscript will be once again subjected to a careful review by referees to determine its suitability for publication.  Your revision should be submitted via ACM Manuscript Central at https://mc.manuscriptcentral.com/toms.  We will expect to receive it within six months time.

Best wishes,

Tim
=======================================================================
Algorithms Editor's Comments to Author:
Algorithms Editor: Hopkins, Tim
Comments to the Author:

In the comments from Reviewer 2 below they refer to a lack of detail about the software. This type of material should NOT appear in the article as it is already supplied in the user manual  (I think). If the reviewer has pointed out something that is missing in the user manual, please update the manual and not the article. 

=======================================================================
Reviewer Comments to Author:
Referee: 1

Recommendation: Needs Major Revision

Comments:
The authors describe HDSDP, a semidefinite programming solver that
uses a two-phase approach.  In the first phase, the SDP is
reformulated using the homogeneous self-dual embedding.  The embedded
problem is then solved to resolve primal and dual feasibility.  In
cases of infeasibility, the self dual embedding provides a certificate
of infeasibility.  If a dual feasible solution is obtained, then the
problem is solved to optimality using a dual interior point method.

Some details of the algorithm are vague and specific parameters are not provided.  For example, in section 4, the authors write "By default S is initialized to a multiple of the identity matrix." but that multiple is not specified.

The paper does not give a detailed analysis of the convergence of the algorithm but since this has already been done for DSDP, that seems
unnecessary.

It should be noted that this software is based heavily on DSDP 5.8, which has been described in Benson and Ye (2008).  The new
contributions here are in the two-phase approach and in the use of a presolver that identifies low-rank structure in constraint matrices.

The software is available on github under an open source MIT license. It would be relatively easy for a reader to download the software and
replicate the results in the paper.  In fact, I did download version 0.9.4 of the software and performed some spot testing.  My results
were consistent with those in the paper.

The authors tested their software on some well known benchmark problems and reported solution times and errors. These results show
that HDSDP is generally faster than DSDP 5.8 and COPT 5.0 on the selected test problems.  The problems were selected to have low-rank
constraint matrices that could be handled with the presolver.  The authors did not present the DIMACS error measures for these test
problems so it isn't possible to determine how accurately the problems were solved.

In an appendix, the authors also give results on the test problems from Hans Mittelmann's benchmark set.  Here the authors provide the
DIMACS error measures.  As with DSDP, HDSDP often produces solutions that have poor primal feasibility (DIMACS errors as large as 1.0e-4.)
In other cases, HDSDP produces solutions with large relative duality gaps as measured by DIMACS errors 5 and 6.  The authors did not
compare the solution times for the Mittelmann benchmark problems with DSDP 5.8 or COPT 5.0.

The paper could be improved significantly with more systematic
computational testing.

1. The authors have not provided computational results demonstrating the effectiveness of the homogeneous self-dual embedding in detecting
primal or dual infeasible problems.  I ran the code on the infp1, infp2, infd1, and infd2 problems in SDPLIB, and found that it failed
to correctly identify infeasibility in all four problems.  Note that these problems are confusingly named, since infp1 is actually dual
infeasible and infd1 is primal infeaible. 

2. DIMACS error measures for the selected test problems in the body of the paper should also be given in the appendix, along with errors for
the other two codes.  Run times for DSDP 5.8 and COPT 5.0 (version 6.0 is available) on the Mittelmann benchmarks should also be given in the
appendix.

3. Readers would presumably like to know whether the performance improvements on the selected test problems were due to the presolver,
the two-phase method, or both.  It would be interesting to run tests with the presolver turned on and off to determine how important it is
to the performance and accuracy of the solver.  Since the second phase is effectively the same as DSDP 5.8, it may not be necessary to
perform that comparison.

Additional Questions:
Please help ACM create a more efficient time-to-publication process: Using your best judgment, what amount of copy editing do you think this paper needs?: Light

Most ACM journal papers are researcher-oriented. Is this paper of potential interest to developers and engineers?: Yes


Referee: 2

Recommendation: Needs Major Revision

Comments:
HDSDP is an interesting contribution to the SDP solver world. The progress demonstrated merits publication. The manuscript, however, needs some improvement before it can be published. 

General remarks:
There should be more details about the software. It seems the code to HDSDP is available as Open-Source under the MIT license at http://github.com/COPT-Public/HDSDP. This is nowhere mentioned in the article. However, at the beginning of Section 5, third-party packages are mentioned. In the code, some stubs for pardiso can be found. It should be made explicitly clear what other software is needed and how everything works together.

The README on the github page has a line saying: "Number of threads available: xx. Optimizing over xx threads." 
In the article, I find nearly nothing about parallel computation. Even on Page 10, Line 518, the number of cores used is not mentioned (possibly 6, using 12 threads?). Here more details are necessary.

The number in Appendix A are peculiar. There are only two significant digits after the decimal point printed and then filled up to five digits with meaningless zeros. If a number like 1.2e-17 is shown, a few more digits should be given for 1.13e-4. The tables should be recomputed. 

Overall, the style, especially in the beginning, is somewhat repetitive. The text should be considerably shortened, for example, P2-69ff, P3-108ff. Some of the unspecific general repetitive praise, how robust etc. the software is, could be omitted. Also, the conclusion is just a repetition of the abstract and could be removed, especially as there is already a summary.

The article would benefit if a language editor would once read through it. I will state a few examples below, but there are more.

Particular remarks:
P1-24: The sentence about the availability of the pre-build binary should be changed to the source code and possibly placed elsewhere.
P2-55: delete “most”.
P2-64: maybe something more scientific than "nice"?
P3-133: “carefully” is very unspecific. 
P3-144: "solved the problem" asks for a short explanation of the implementation's outcome, taking numerical aspects into account. What precisely is solved?
P3-147 “a primal solution”
P3-151 “is” -> “are not free”
P3-154 “property” -> “properly”?
P4-158 “penalties to ensure a nonempty”
P4-159 “These attempts works well” -> “This works well”
P4-167 “the big-M method”
P6-293/294 please rephrase
P7-313 “Last” -> “Finally, “
P7-360 Be more explicit about the computational tricks and the 3rd party packages.
P10-514 give a reference (URL) and a date regarding the statement about Hans Mittelmann’s benchmark.
P11 Table 2, P12 Table 3,4: Times would be better readable if flush right.
P12 581 an explanation of why gpp250-3 fails on all solvers would be helpful. This is also true in general for the other failures.
P13-631 “enjoy” is an interesting choice here.
P13-653 delete “vast”
P13-658-661 rephrase 
P13-663 “enjoy” and “friendly” are, as “nice” before, words that convey little information to the reader.
P14-686 rephrase to, e.g., “To conclude, HDSDP is only effective for certain SDPs.” 



Additional Questions:
Please help ACM create a more efficient time-to-publication process: Using your best judgment, what amount of copy editing do you think this paper needs?: Moderate

Most ACM journal papers are researcher-oriented. Is this paper of potential interest to developers and engineers?: Yes