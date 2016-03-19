### Code to accompany *Quantum steering: a short review with focus on semi-definite programming*
#### Daniel Cavalcanti and Paul Skrzypczyk

This repository provides a small collection of code which implements many of the semidefinite programs presented in the review article "*Quantum steering: a short review with focus on semi-definite programming*".

All code is written in MATLAB and requires:
- [CVX](http://cvxr.com/) - a Matlab-based convex modeling framework
- [QETLAB](http://www.qetlab.com/) - A MATLAB Toolbox for Quantum Entanglement

It has been tested on Matlab R2014a, and CVX 2.1 

The code comprises the following:
- Membership, helper, and misc:
  - [NSAssemblage](https://github.com/paulskrzypczyk/steeringreview/blob/master/NSAssemblage.m): determines whether a bipartite assemblage is a valid non-signalling assemblage or not<sup>§</sup>.
  - [LHSAssemblage](https://github.com/paulskrzypczyk/steeringreview/blob/master/LHSAssemblage.m): determines whether a biparitte assemblage has an LHS model or not<sup>§</sup>.
  - [genAssemblage](https://github.com/paulskrzypczyk/steeringreview/blob/master/genAssemblage.m): generates an assemblage starting from a quantum state and a set of measurements.
  - [validPOVMs](https://github.com/paulskrzypczyk/steeringreview/blob/master/validPOVMs.m): determines whether a set of POVMs is valid or not.
  - [JMPOVMs](https://github.com/paulskrzypczyk/steeringreview/blob/master/JMPOVMs.m): determines whether a set of measurements is jointly measurable or not.
  - [genRandProjMeas](https://github.com/paulskrzypczyk/steeringreview/blob/master/genRandProjMeas.m): generates a random set of projective measurements
  - [genSinglePartyArray](https://github.com/paulskrzypczyk/steeringreview/blob/master/genSinglePartyArray.m): generates the single-party determinstic probability distributions
  - [findRadiusPolytopeInBlochSphere](https://github.com/paulskrzypczyk/steeringreview/blob/master/findRadiusPolytopeInBlochSphere.m):
  
<sup>§</sup>: These files can be used inside CVX as a means to enforce the corresponding constraint. 
