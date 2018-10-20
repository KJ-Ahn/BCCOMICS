# BCCOMICS
BCCOMICS - Baryon CDM COsMological Initial Condition generator for Small scales.

Small-scale fluctuations in the early universe, even at redshift as high as z=200, are found strongy affected by large-scale density and streaming-velocity environments. BCCOMICS is an initial condition generator that allows study of structure formation inside a simulation box of < 4 comoving Mpc, where the simulation box can have non-zero overdensity (Delta) and streaming velocity (V_cb = mean velocity of CDM - mean velocity of baryon) as its environmental condition. This allows study of cosmic variance of e.g. first star formation under varying large-scale environments.

Currently, it only supports [enzo](http://enzo-project.org). We are inviting contributors to help port this to [Gadget](https://wwwmpa.mpa-garching.mpg.de/gadget/), [RAMSES](https://bitbucket.org/rteyssie/ramses/overview), and many other fabulous N-body+hydro simulation codes.

The main code is composed of two parts.

(1) Realization of large-scale fluctuations (at a length resolution of 4 Mpc comoving) and transfer function calculation (bccomics_setup.m):  
Because the perturbation theory studying the dual impact from large-scale Delta and V_cb on small-scale fluctuations is relatively new ([Ahn 2016](http://adsabs.harvard.edu/abs/2016ApJ...830...68A)), transfer functions from usual linear Boltzmann solvers such as CAMB do not provide the level of accuracy of this new theory. Even when Delta variance is not considered, the suppression of high-k (around k~100/Mpc) modes due to V_cb ([Tseliakhovich & Hirata 2010](http://adsabs.harvard.edu/abs/2010PhRvD..82h3520T)) is not reflected in transfer function from CAMB. "bccomics_setup.m" makes realization of environmental variables, and once the user chooses a specific "patch" of generically non-zero Delta and V_cb, it calculates and records the transfer function.

(2) Realization of small-scale fluctuations on a selected patch (bccomics.m):  
The user is asked again to choose one from already calculated set of patches, and then the corresponding transfer function is read in and used to generate 3D data of CDM and baryons. Currently, following data in simple bianry form and in enzo units are generated:  
(a) CDM particle positions (cpos1, cpos2, cpos3)  
(b) CDM particle velocities (vc1, vc2, vc3)  
(c) baryon grid density (db)  
(d) baryon grid velocities (vb1, vb2, vb3)  
(e) baryon grid thermal energy (etherm)  
(f) baryon grid kinetic+thermal energy (etot)  
(g) baryon particle positions (bpos1, bpos2, bpos3) (optional; with SPH in mind)  
(h) baryon particle velocities (vpb1, vpb2, vpb3) (optional; with SPH in mind)

*In case of MATLAB, initial conditions in hdf5 binary, ready to be used for enzo, can be generated.

(3) (OCTAVE-only) Conversion of binary data from step (2) into enzo-usable initial conditions, by "convert_enzo.py" using python+h5py. 

## Installation and Requirements

Installation of BCCOMICS:  
Either (1) clone this repo (`git clone https://github.com/KJ-Ahn/BCCOMICS.git`), or (2) download as a zip and extract its contents. Click the green "Clone or download" button and choose whichever suits you. No further installation process is required.

BCCOMICS has two main scripts that run on either [MATLAB(R)](https://www.mathworks.com/products/matlab.html) or [gnu OCTAVE](https://www.gnu.org/software/octave/). gnu OCTAVE is easily installed with its dependency by package managers (e.g. "apt", "rpm", "synaptic", "flatpak", ...) in usual linux distributions. Ask your system administrator for installation on a shared unix machine.

To install a rather recent version of OCTAVE on old linux distributions (or to install it privately on a shared machine without root priviledge), you might want to try installing [Anaconda](https://www.anaconda.com/download/), and then do `conda install -c conda-forge octave` for installing a recent OCTAVE with its dependencey without headache. 

For MATLAB, in addition to the main program, following additional packages need to be installed.  
- Image Processing Toolbox  
- Statistics and Machine Learning Toolbox

For gnu OCTAVE, in addition to the main program, following additional packages need to be installed. Use your linux distribution's package installer (e.g. `sudo apt install octave-image` for Ubuntu). Ask your system administrator for installation on a shared unix machine.  
- octave-image  
- octave-statistics (optional; if uninstalled BCCOMICS will use functions under BCCOMICS/statistics-1.3.0/)  
- octave-odepkg (optional; if uninstalled BCCOMICS will use octave functions or those under BCCOMICS/odepkg-0.8.5/)
- python
- h5py

BCCOMICS needs transfer function output files of [CAMB](https://camb.info/) (as provided in BCCOMICS/sample/CAMB_for_mode_finding), and also a redshift-ionization output file from [RECFAST](http://www.astro.ubc.ca/people/scott/recfast.html). So in general, you need to install  
- CAMB
- RECFAST

## Running
Best explained with an example. Let's assume that BCCOMICS is installed at `/home/kjahn/BCCOMICS`, whose sub-directories are `src`, `sample`, etc. (Or, assume BCCOMICS installed at `c:\Documents\BCCOMICS` on a Windows system)  
Below `$` is a linux command prompt, `>>` is either OCTAVE's or MATLAB's command prompt.

### (1) Set a work directory
First, you need a work directory under which "params.m" and "params_patch.m" exist. The name of this directory can be anything. 
Inside OCTAVE/MATLAB, you need to go to this directory by `cd` command. You also need to add the source path by `addpath` commmand. Current example "params.m" generates 151^3 unigrid patches inside (604 Mpc)^3 volume with bccomics_setup. Current example "params_setup.m" generates initial conditions for CDM and baryons with 64^3 resolution. For your own setup, copy `sample` directory to e.g. `my_params`, and modify "params.m" and "params_patch.m" which are both self-explanatory.
#### (for OCTAVE on linux machine)
```
$ octave

>> cd '/home/kjahn/BCCOMICS/sample'
>> addpath('/home/kjahn/BCCOMICS/src')
```
#### (for MATLAB on Windows machine)
```
>> cd 'c:\Documents\BCCOMICS\sample'
>> addpath('c:\Documents\BCCOMICS\src')
```

### (2) Run bccomics_setup
You will be asked to choose a patch with your desired CDM overdensity and V_cb (absolute value).
```
>> bccomics_setup
>>
```

### (3) Run bccomics
You will be asked to choose one among those patches that you have chosen in more-than-one bccomics_setup runs.
```
>> bccomics
```
After all runs are successuful, you will have initial conditions under `ICdir` with appropriate directory name.

### (4) For OCTAVE only, convert output binary ICs to enzo-readable HDF5 binary files.
```
$ python convert2enzo.py
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
