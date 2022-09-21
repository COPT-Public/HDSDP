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

    <doc-data|<doc-title|Dimension-reduced Interior Point
    Method>|<doc-author|<author-data|<\author-affiliation>
      \;

      Discussion 6

      \;

      \;

      <date|>
    </author-affiliation>>>>
  </hidden>|<\hidden>
    <tit|Negative Curvature of Hessian of Potential>

    <strong|Accelerate computation of negative curvature>

    \;

    <\equation*>
      <n><rsup|2><rsub|<x><x>>\<varphi\><around*|(|<x>|)>=-<frac|\<rho\><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)>|f<around*|(|<x>|)><rsup|2>><rsup|\<top\>>+<frac|\<rho\><A><rsup|\<top\>><A>|f<around*|(|<x>|)>>+<X><rsup|-2>.
    </equation*>

    \;

    Observations

    \;

    <\itemize>
      <item><math|<n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>>>
      is dense and specially treated

      <item>Customized Lanczos solver behaves 100x more efficient

      <space|1em>than direct eigen-decomposition

      <item>Works fine for problems up to <math|10<rsup|-3>> accuracy

      <item>For higher accuracy, Lanczos might be unstable
    </itemize>
  </hidden>|<\hidden>
    <tit|Numerical difficulty>

    \;

    The spectrum of <math|\<nabla\><rsup|2>> is increasingly worse when
    optimization proceeds

    <\big-figure|<image|file:///Users/gaowenzhi/Desktop/potred/doc/discussion/fig1.png|600px|||><image|file:///Users/gaowenzhi/Desktop/potred/doc/discussion/fig2.png|600px|||>>
      Spectrum at the beginning/at the end, ploted in relative scale

      Left: min <math|-9.36\<times\>10<rsup|5>>, max
      <math|2.24\<times\>10<rsup|11>>. Right: min
      <math|-3.18\<times\>10<rsup|5>>, max <math|1.26\<times\>10<rsup|14>>
    </big-figure>

    The positive eigen-value of the potential function increases and makes
    negative hard to find.
  </hidden>|<\hidden>
    <tit|Improvement by scaling>

    \;

    <math|<n><rsup|2><rsub|<x><x>>\<varphi\><around*|(|<x>|)>=-<frac|\<rho\><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)>|<with|color|red|f<around*|(|<math-bf|x>|)><rsup|2>>><rsup|\<top\>>+<frac|\<rho\><A><rsup|\<top\>><A>|<with|color|red|f<around*|(|<math-bf|x>|)>>>+<with|color|red|<math-bf|X><rsup|-2>>>.

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<v>>>|<cell|<around*|\<langle\>|<v>,<n><rsup|2><v>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<e><rsup|\<top\>><v>=0>|<cell|>>|<row|<cell|>|<cell|<around*|\<\|\|\>|<v>|\<\|\|\>>=1>|<cell|>>>>
    </eqnarray*>

    An alternative is to use scaled Hessian\ 

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<v>\<neq\><0>>>|<cell|<around*|\<langle\>|<around*|(|<I>-<frac|<x><x><rsup|\<top\>>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>|)><X><v>,<n><rsup|2><X><around*|(|<I>-<frac|<x><x><rsup|\<top\>>|<around*|\<\|\|\>|<x>|\<\|\|\>><rsup|2>>|)><v>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|>|<cell|>>|<row|<cell|>|<cell|<around*|\<\|\|\>|<v>|\<\|\|\>>=1>|<cell|>>>>
    </eqnarray*>

    and the scaled Hessian is\ 

    <\equation*>
      -<frac|\<rho\><X><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)><rsup|\<top\>><X>|<with|color|red|f<around*|(|<math-bf|x>|)><rsup|2>>>+<frac|\<rho\><X><A><rsup|\<top\>><A><X>|<with|color|red|f<around*|(|<math-bf|x>|)>>>+<I>.
    </equation*>
  </hidden>|<\hidden>
    <tit|Further improvement>

    \;

    <\equation*>
      <n><rsup|2><rsub|<x><x>>\<varphi\><around*|(|<x>|)>=-<frac|\<rho\><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)>|f<around*|(|<x>|)><rsup|2>><rsup|\<top\>>+<frac|\<rho\><A><rsup|\<top\>><A>|f<around*|(|<x>|)>>+<X><rsup|-2>
    </equation*>

    <\itemize>
      <item>The scaled direction makes spectrum more balanced

      <item>Another way is to do diagonal scaling
      <math|diag<around*|(|<n><rsup|2><rsub|<x><x>>\<varphi\><around*|(|<x>|)>|)>>

      <item>Both solving LPs to <math|10<rsup|-4>> to <math|10<rsup|-5>> (but
      accurate eigen computation gives <math|10<rsup|-7>> to
      <math|10<rsup|-8>>)
    </itemize>

    \;

    More attempts include

    <big-table|<block|<tformat|<cwith|1|-1|1|-1|cell-tborder|1ln>|<cwith|1|-1|1|-1|cell-bborder|1ln>|<cwith|1|-1|1|-1|cell-lborder|0ln>|<cwith|1|-1|1|-1|cell-rborder|0ln>|<cwith|1|-1|1|-1|cell-halign|c>|<table|<row|<cell|Method>|<cell|Accuracy>>|<row|<cell|Lanczos>|<cell|<math|10<rsup|-4>>>>|<row|<cell|Lanczos
    + scaling (+ re-orthorgonalization)>|<cell|<math|10<rsup|-4>\<sim\>10<rsup|-5>>>>|<row|<cell|Shifted
    inverse power method >|<cell|<math|10<rsup|-5>> (solve linear
    system)>>|<row|<cell|Full eigen-decomposition>|<cell|<math|10<rsup|-7>\<sim\>10<rsup|-8>>>>>>>|Current
    attempts>
  </hidden>|<\hidden>
    <tit|One more cure: learning the support of negative curvature>

    \;

    \;

    The spectrum of <math|\<nabla\><rsup|2>> is increasingly worse when
    optimization proceeds

    \;

    <\equation*>
      <n><rsup|2><rsub|<x><x>>\<varphi\><around*|(|<x>|)>=-<frac|\<rho\><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)>|<with|color|red|f<around*|(|<math-bf|x>|)><rsup|2>>><rsup|\<top\>>+<frac|\<rho\><A><rsup|\<top\>><A>|<with|color|red|f<around*|(|<math-bf|x>|)>>>+<with|color|red|<math-bf|X><rsup|-2>>.
    </equation*>

    \;

    Intuition: <math|arg min<rsup|k><rsub|i> <around*|{|x<rsub|i>|}>> is
    <strong|not likely> to contribute to negative curvature.

    \;

    \;

    <\itemize>
      <item>Use <math|<math-bf|X><rsup|-2>> to predict the
      <with|color|red|support> of curvature <math|<v>>

      <item>Truncate <math|\<alpha\>> percent of the Hessian by
      <math|<math-bf|X><rsup|-2>>

      <item>Reducing complexity and hardness to evaluate Hessian
    </itemize>
  </hidden>|<\hidden>
    <tit|Reduced support for Hessian>

    \;

    \;

    <big-figure|<image|file:///Users/gaowenzhi/Desktop/potred/doc/discussion/fig3.png|600px|||><space|1em><image|file:///Users/gaowenzhi/Desktop/potred/doc/discussion/fig4.png|600px|||>|Spectrum
    of reduced support truncating 50% of dimension>

    \;

    <\itemize>
      <item>Retain negative curvature in general

      <item>Need more tuning
    </itemize>
  </hidden>|<\hidden>
    <tit|Conclusions and future improvement>

    \;

    <math|<n><rsup|2>=-<frac|\<rho\><n>f<around*|(|<x>|)><n>f<around*|(|<x>|)>|<with|color|red|f<around*|(|<math-bf|x>|)><rsup|2>>><rsup|\<top\>>+<frac|\<rho\><A><rsup|\<top\>><A>|<with|color|red|f<around*|(|<math-bf|x>|)>>>+<with|color|red|<math-bf|X><rsup|-2>>>.

    <\eqnarray*>
      <tformat|<table|<row|<cell|min<rsub|<v>>>|<cell|<around*|\<langle\>|<v>,<n><rsup|2><v>|\<rangle\>>>|<cell|>>|<row|<cell|<text|subject
      to>>|<cell|<e><rsup|\<top\>><v>=0>|<cell|>>|<row|<cell|>|<cell|<around*|\<\|\|\>|<v>|\<\|\|\>>=1>|<cell|>>>>
    </eqnarray*>

    <\itemize>
      <item>The customization greatly accelerate computation

      100x faster at the beginning

      <item>The spectrum of Hessian is not well-distributed\ 

      Lanczos stagnates since <math|\<lambda\><rsub|1><around*|(|\<nabla\><rsup|2>|)>\<gg\>\<lambda\><rsub|n-k><around*|(|\<nabla\><rsup|2>|)>>,
      for <math|k> of interest

      <item>Necessary to look for a method that uses the prior knowledge to
      accelerate

      We know that there exists exactly <with|color|red|one> negative
      eigen-value

      (Now experimenting on FEAST contour method)
    </itemize>

    \;
  </hidden>|<\shown>
    <tit|>

    \;

    <\equation*>
      \<varphi\><around*|(|<x>,<s>|)>=\<rho\>log<around*|(|<x><rsup|\<top\>><s>|)>-log<around*|(|<s>|)>-log<around*|(|<x>|)>
    </equation*>
  </shown>>
</body>

<\initial>
  <\collection>
    <associate|math-color|black>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|1|1>>
    <associate|auto-3|<tuple|2|?>>
    <associate|auto-4|<tuple|3|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<\surround|<hidden-binding|<tuple>|1>|>
        Spectrum at the beginning/at the end, ploted in relative scale

        Left: min <with|color|<quote|black>|font-family|<quote|rm>|<with|mode|<quote|math>|-9.36\<times\>10<rsup|5>>>,
        max <with|color|<quote|black>|font-family|<quote|rm>|<with|mode|<quote|math>|2.24\<times\>10<rsup|11>>>.
        Right: min <with|color|<quote|black>|font-family|<quote|rm>|<with|mode|<quote|math>|-3.18\<times\>10<rsup|5>>>,
        max <with|color|<quote|black>|font-family|<quote|rm>|<with|mode|<quote|math>|1.26\<times\>10<rsup|14>>>
      </surround>|<pageref|auto-1>>
    </associate>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Current
      attempts>|<pageref|auto-2>>
    </associate>
  </collection>
</auxiliary>