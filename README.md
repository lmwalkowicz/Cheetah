Cheetah
=======

This repo contains code + etc. for modeling starspots in photometric data (lightcurves).  

Contributors: 
Lucianne Walkowicz (Princeton Dept of Astrophysical Sciences)
Michael Thomas (Princeton 'XX --> what is your year, MT?)
Adam Finkestein (Princeton Computer Science Dept)

NEXT STEPS (AF notes from 2014-05-28):
- try > 2 spots with synth data
- range of periods
- cluster to run multiple experiments in parallel (ionic?)
- other optimization strategies instead of just L-M / nonlinear simplex...? (just try stuff built into python libraries)

OLD Plan of action:

Code Logistics:
- AF + MT find their GitHub passwords and send them to LW
- LW gives AF + MT push access to the repo
- MT pushes existing code to the repo

Lightcurve Fitting:

- Begin by fitting synthetic/artificial lightcurves generated by the code itself

- Possible to take derivatives of the analytic spot functions in Eker, e.g. using Mathematica or Maple? 

- Use python/scipy built-in optimizers (e.g. nonlinear simplex, gradient methods, etc)

- MCMC implementation - what are the confidence intervals on the fitted parameters?

- number of spots: start with 1, expand to 2 or 3. 

- Revisit issue of how to do multispot fitting (i.e. avoiding overlapping spots)

- Eventually apply model to real data, first targeting stars whose inclination is constrained to lie within a particular range. 

