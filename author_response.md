**Author reponse to submission Algorithm xxxx: HDSDP: Software for Semi-definite Programming** 

We thank the reviewers for their efforts in the review process and their valuable suggestions. 

To fully address the reviewers' concerns, in the past 6 month we've been working on a new version of HDSDP with more elegant and efficient software engineering 

**Response to questions of Reviewer 1**

1. Details of the algorithm are vague and specific parameters are not provided

   We've added 

2. The convergence analysis

   [TODO: The nonsymmetric cone paper]

3. This software is based heavily on DSDP 5.8

4. DIMACS errors of reported test problems accuracy

   1e-06

5. More detailed Comparison with DSDP5.8 and COPT5.0 on Hans Mittelman's benchmark 

6. Infeasibility detection module fails to work

**Response to suggestions of Reviewer 1**

1. Efficiency of detecting infeasibility using self-dual embedding
2. DIMACS error measures for the selected test problems
3. Ablation study on the presolver 

----

**Response to questions of Reviewer 2**

1. The details of HDSDP software and third-party dependency

   By ACM TOMs

2. Multi-threading

3. Numbers in Appendix

   Zeros and 1e-17

4. Stylistic issues.

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


