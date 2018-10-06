# BCCOMICS
BCCOMICS - Baryon CDM COsMological Initial Condition Generator for Small scales.

Small-scale fluctuations in the early universe, even at redshift as high as z=200, are found strongy affected by large-scale density and streaming-velocity environments. BCCOMICS is an initial condition generator that allows study of structure formation inside a simulation box of < 4 comoving Mpc, where the simulation box can have non-zero overdensity (Delta) and streaming velocity (V_cb = mean velocity of CDM - mean velocity of baryon) as its environmental condition. This allows study of cosmic variance of e.g. first star formation under varying large-scale environments.

The main code is composed of two parts. 
(1) Realization of large-scale fluctuations (at a length resolution of 4 Mpc comoving) and transfer function calculation (bccomics_setup.m):
Because the perturbation theory studying the dual impact from large-scale Delta and V_cb on small-scale fluctuations is relatively new ([Ahn 2016](http://adsabs.harvard.edu/abs/2016ApJ...830...68A)), transfer functions from usual linear Boltzmann solvers such as CAMB do not provide the level of accuracy of this new theory. Even when Delta variance is not considered, the suppression of high-k (around k~100/Mpc) modes due to V_cb ([Tseliakhovich & Hirata 2010](http://adsabs.harvard.edu/abs/2010PhRvD..82h3520T)) is not reflected in transfer function from CAMB. "bccomics_setup.m" makes realization of enrionmental variables, and once the user chooses a specific "patch" of generically non-zero Delta and V_cb, it calculates and records the transfer function.
(2) Realization of small-scale fluctuations on a selected patch (bccomics.m):
blah blah...

## Installation Requirements


BCCOMICS is a set of scripts that run on either MATLAB(R) or gnu OCTAVE. Provide links like [this](url).

## Running

Here is an example of how to format text like source code:
```
do something here.
```

## Citing

If you use `BCCOMICS` in your work, please cite as:

```
...the BCCOMICS initial conditions generator (Ahn & Smith 2018).
```

BIBTEX:
```
@ARTICLE{2018arXiv180704063A,
   author = {{Ahn}, K. and {Smith}, B.~D.},
    title = "{Formation of First Galaxies inside Density Peaks and Voids under the Influence of Dark Matter-Baryon Streaming Velocity, I: Initial Condition and Simulation Scheme}",
  journal = {ArXiv e-prints},
archivePrefix = "arXiv",
   eprint = {1807.04063},
 keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics},
     year = 2018,
    month = jul,
   adsurl = {http://adsabs.harvard.edu/abs/2018arXiv180704063A},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
