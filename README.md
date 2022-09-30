# Robust_EM
## MENU
Coding files are placed under CODE repository.
R:
* Demo.Rmd demonstrates the current build of robust EM. 
* The algorithm has two part: 1. initial_cond: It generates the intial values to input into rem_core. 2. rem_core is the EM part of the algorithm.
* plot_helper and simulator are two helper files that are responsible for plotting and simulating points.

Python:
* gn_attempt is the numerical solution to the system of equations derived in constrain_struct

For Math derivation, please check the two pdf files named as REM_manuscript.pdf, and constrain_struct.pdf.

Github & Dropbox:
https://www.dropbox.com/sh/nil3ma4s0x2c2yl/AABr1mlR5ERW-N8q4dFJ4r91a?dl=0 (This is the original build. It has many helper functions, but many of them are broken and can't be used directly)
https://github.com/mkbwang/RobustEM (This is a package build with very limited features)
https://github.com/AmIACommonGuy/Robust_EM
