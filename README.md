# pcpsim

A collection of MATLAB/OCTAVE scripts to perform simulations of homogeneous
precipitation in a 2-phase alloy.

A system of differential-algebraic equations (DAEs) that describe the nucleation
and growth of precipitates and the solute balance in the alloy are solved using
the standard solvers `ode15i` (MATLAB/OCTAVE) or `daspk` (OCTAVE only).

The scripts may be used for simulating precipitate evolution during heat-treatment processes.

For more information read the documentation: [pcpsim.pdf](doc/pcpsim.pdf) 

## Contents

  - `./m`   The basic script files
  - `./doc` Documentation
  - `./FeC` An example for carbide precipitation in Fe-C alloy
  - `./FeN` Some examples for nitride precipitation in Fe-N



