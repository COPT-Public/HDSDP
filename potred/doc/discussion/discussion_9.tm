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

      Discussion 9

      \;

      \;

      <date|>
    </author-affiliation>>>>
  </hidden>|<\shown>
    <tit|A Brief Overview on Implementation Progress>

    \;

    Done:

    <\itemize>
      <item>LP interface

      <item>Trust-region solver

      <item>Potential reduction framework

      Gradient-momentum type
    </itemize>

    Doing:

    <\itemize>
      <item>Lanczos negative curvature method

      <item>The support prediction method

      <item>Debugging

      <item>Abstract Ruiz-scaling
    </itemize>

    \;

    <strong|Performance>

    <\itemize>
      <item>20 to 50 times faster than raw MATLAB implementation

      <item>Higher accuracy using specialized implementation tricks
    </itemize>
  </shown>>
</body>

<\initial>
  <\collection>
    <associate|font-base-size|9>
    <associate|math-color|black>
    <associate|page-medium|paper>
  </collection>
</initial>