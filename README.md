# Grain boundary-stress fundamental zones
Computational tools to make use of grain boundary-stress orientation fundamental zones of cubic metals with Oh point group symmetry, for axisymmetric and triaxial stresses. The functions in this toolbox allow for the development of quasi-equidistant grids in disorientation, boundary plane and stress fundamental zones, to map arbitrary configurations into their symmetrically equivalent ones within the fundamental zone, and to snap configurations into their closest condition within a predefined grid. This also includes an implementation of an elastic bicrystal incompatibility stress model and calculation of damage metrics.

## Requirements
Coded in Matlab R2022a. 
Note: Function conflicts may occur if the MTEX toolbox is initialized.

## Dependencies
The 'inhull.m' function (https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull/files/inhull.m) by John D'Errico is included in this repository. 

## License
This repository is published under an MIT license ([![License: MIT](https://img.shields.io/badge/License-MIT.svg)](https://www.mit.edu/\~amini/LICENSE.md)). See `LICENSE` for full licensing information.

## Contents
- [**functions/**](./functions/): MATLAB functions.
- [**Examples/**](./Examples/): Examples for the generation and usage of fundamental zones.

