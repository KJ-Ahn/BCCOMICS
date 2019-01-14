# BCCOMICS
BCCOMICS - Baryon CDM COsMological Initial Condition generator for Small scales.

Small-scale fluctuations in the early universe, even at redshifts as high as z=200, are found to be strongly affected by large-scale density and streaming-velocity environments. BCCOMICS is an initial condition generator that allows the study of structure formation inside a simulation box of < 4 comoving Mpc, where the simulation box can have non-zero overdensity (&Delta;) and streaming velocity (V<sub>cb</sub> = mean velocity of CDM - mean velocity of baryon) as its environmental condition. This allows for the study of cosmic variance of e.g., first star formation under varying large-scale environments.

Currently, it only supports [enzo](http://enzo-project.org). We are inviting contributors to help port this to [Gadget](https://wwwmpa.mpa-garching.mpg.de/gadget/), [RAMSES](https://bitbucket.org/rteyssie/ramses/overview), and many other fabulous N-body+hydro simulation codes.

The main code is composed of two parts.

1. Realization of large-scale fluctuations (at a length resolution of 4 Mpc comoving) and transfer function calculation (`bccomics_setup.m`):  
Because the perturbation theory studying the dual impact from large-scale &Delta; and V<sub>cb</sub> on small-scale fluctuations is relatively new ([Ahn 2016](http://adsabs.harvard.edu/abs/2016ApJ...830...68A)), transfer functions from the usual linear Boltzmann solvers, such as CAMB, do not provide the level of accuracy of this new theory. Even when &Delta; variance is not considered, the suppression of high-k (around k~100/Mpc) modes due to V<sub>cb</sub> ([Tseliakhovich & Hirata 2010](http://adsabs.harvard.edu/abs/2010PhRvD..82h3520T)) is not reflected in the transfer function from CAMB. "bccomics_setup.m" makes a realization of environmental variables, and once the user chooses a specific "patch" of generically non-zero &Delta; and V<sub>cb</sub>, it calculates and records the transfer function.

2. Realization of small-scale fluctuations on a selected patch (`bccomics.m`):  
The user is asked again to choose one from a previously-calculated set of patches, and then the corresponding transfer function is read in and used to generate 3D data of CDM and baryons. Currently, the following data are generated:  
a. CDM particle positions (cpos1, cpos2, cpos3)  
b. CDM particle velocities (vc1, vc2, vc3)  
c. baryon grid density (db)  
d. baryon grid velocities (vb1, vb2, vb3)  
e. baryon grid thermal energy (etherm)  
f. baryon grid kinetic+thermal energy (etot)  
g. baryon particle positions (bpos1, bpos2, bpos3) (optional, for SPH)  
h. baryon particle velocities (vpb1, vpb2, vpb3) (optional, for SPH)  
i. baryon particle thermal energy (eptherm) (optional, for SPH)  
Data is written in [Enzo's](https://enzo.readthedocs.io) internal unit system.  
**If using MATLAB, initial conditions can be saved directly to Enzo's HDF5 format. Saving to HDF5 is not supported in OCTAVE, but a Python script for converting the native output to Enzo format is provided.**  
**Items _g_, _h_, and _i_ can be written when `baryonparticleflag=true` and `particlevelocity_accuracyflag=true` in `params_patch.m`.**

3. (OCTAVE-only) Conversion of binary data from step (2) into enzo-usable initial conditions, by "convert_enzo.py" using python+h5py. 

## Installation and Requirements

Installation of BCCOMICS:  
Either (1) clone this repo (`git clone https://github.com/KJ-Ahn/BCCOMICS.git`), or (2) download as a zip and extract its contents. Click the green "Clone or download" button and choose whichever suits you. No further installation process is required.

BCCOMICS has two main scripts that run on either [MATLAB(R)](https://www.mathworks.com/products/matlab.html) or [gnu OCTAVE](https://www.gnu.org/software/octave/). gnu OCTAVE is easily installed with its dependency by package managers (e.g. "apt", "rpm", "synaptic", "flatpak", ...) in usual linux distributions. Ask your system administrator for installation on a shared unix machine.

To install a rather recent version of OCTAVE on old linux distributions (or to install it privately on a shared machine without root priviledge), you might want to try installing [Anaconda](https://www.anaconda.com/download/), and then do `conda install -c conda-forge octave` for installing a recent OCTAVE with its dependencey without headache. 

For MATLAB, in addition to the main program, it is preferred (but not essential) that following additional packages be installed.  
- Image Processing Toolbox (optional; if uninstalled BCCOMICS will use `padarray.m` under BCCOMICS/mfiles_for_matlab/ImagProc_Tool)  
- Statistics and Machine Learning Toolbox (optional; if uninstalled BCCOMICS will use functions under BCCOMICS/mfiles_for_matlab/Stat_Tool)

For gnu OCTAVE, in addition to the main program, following additional packages need to be installed. Use your linux distribution's package installer (e.g. `sudo apt-get install octave-image` for Ubuntu). It is OK not to install following "octave-*" packages if you want, BCCOMICS is shipped with necessary function files of these packages anyways.  
- octave-image (optional; if uninstalled BCCOMICS will use `padarray.m` under BCCOMICS/mfiles_for_octave/) 
- octave-statistics (optional; if uninstalled BCCOMICS will use functions under BCCOMICS/statistics-1.3.0/)  
- octave-odepkg (optional; if uninstalled BCCOMICS will use octave functions or those under BCCOMICS/odepkg-0.8.5/)
- python
- h5py

BCCOMICS needs transfer function output files of [CAMB](https://camb.info/) (as provided in BCCOMICS/sample/CAMB_for_mode_finding), and also a redshift-ionization output file from [RECFAST](http://www.astro.ubc.ca/people/scott/recfast.html) (as provided as the file BCCOMICS/sample/output_recfast). So in general, you need to install  
- CAMB
- RECFAST

## Running
Best explained with an example. Let's assume that BCCOMICS is installed at `/home/kjahn/BCCOMICS`, whose sub-directories are `src`, `sample`, etc. (Or, assume BCCOMICS installed at `c:\Documents\BCCOMICS` on a Windows system)  
Below `$` is a linux command prompt, `>>` is either OCTAVE's or MATLAB's command prompt.

### (1) Set a work directory
First, you need a work directory under which "params.m" and "params_patch.m" exist (You should stick to this naming convention!!), and also a cosmology parameter file (as provided as BCCOMICS/sample/LCDM.m). The name of the work directory can be anything. Inside OCTAVE/MATLAB, you need to go to this directory by the `cd` command. You also need to add the source path by the `addpath` commmand.  
The current example "params.m" generates 151<sup>3</sup> unigrid patches inside a (604 Mpc)<sup>3</sup> volume with `bccomics_setup.m`. The example "params_patch.m" generates initial conditions of CDM and baryons with 64<sup>3</sup> resolution with `bccomics.m`.  
For your own setup, do a recursive copy of the `sample` directory to e.g. `my_params`, and modify "params.m" and "params_patch.m" which are both self-explanatory. The cosmology parameter file (e.g. LCDM.m) should carefully reflect parameters you used for running CAMB and RECFAST.  
**If you want to use the pre-generated (by bccomics_setup.m) gaussian random seed ("gaussseed.matbin"; stick to this naming convention!!), place it under `setupdir` specified in params.m.**
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

### (2) Run bccomics_setup (when patchidxinput_flag=false in params.m)
You will be asked to choose a patch with your desired CDM overdensity and V<sub>cb</sub> (absolute value).
```
>> bccomics_setup
... Several message outputs ...

Standard deviation of CDM overdensities (sDc) is 0.0041819
Choose CDM overdensity environment:
Input 0 for mean, 1 for overdense, 2 for underdense:
```
Well, if you are interested in overdense patch, enter 1
```
Input 0 for mean, 1 for overdense, 2 for underdense:1
What multiple of sDc away from the mean overdensity, 0? Example: for Delta_c = +1.5*sDc, Enter 1.5
Enter a floating-point number:
```
If you want 2sigma density environment for your patch, enter 2
```
Enter a floating-point number:2
CDM overdensity chosen: Delta_c = 2*sDc = 0.0083638
---------------------------------------
RMS of Vbc (rmsV) at z = 1000 is 27.528 km/s
Peak of Vbc in Maxwell-Boltzmann distribution is 22.477 km/s
Choose Vbc environment at z = 1000
Enter Vbc at z = 1000 in units of km/s: 
```
If you are interested in typical ones, from above V<sub>cb</sub>=30 km/s is a good choice, so enter 30
```
Enter Vbc at z = 1000 in units of km/s: 30
311 patches out of total 3.44295e+06 patches satisfy your chosen condition with 1% margin.
----------------One best matching patch is being found----------------
Wanted Delta_c = 0.0083638; Selected patch's Delta_c = 0.0083689
Wanted Vcb = 30 km/s; Selected patch's Vcb = 30.011 km/s
----------------Integrating----------------
1th wavenumber out of 100 is being handled.
* 1th angle outta 21 is being handled.
* 2th angle outta 21 is being handled.
...
```
Then you should wait until all wavenumbers are calculated. Be patient, it takes a while.
### (3) Run bccomics
Necessary parameters for bccomics is in "params_patch.m". It is mostly self-explanatory. Just a bit more explanation on gaussian random seed: for the final initial condition, you need to either (a) use preexisting seed or (b) generate a new seed. If (a) is the case, `olseedflag = true`, and you should specify `diroldseed` and `Noldseed`. The gaussian random seed is named as e.g. `subgaussseed512.matbin` if `Noldseed=512`, and the file has `Noldseed*Noldseed*(Noldseed/2+1)` complex values. 

After setting "params_patch.m", start octave or matlab, go to the work directory, add path of the src directory, and then run bccomics:
```
$ octave  

>> cd '/home/kjahn/BCCOMICS/sample'
>> addpath('/home/kjahn/BCCOMICS/src')
>> bccomics
```
You will be asked to choose one among those patches that you have chosen in more-than-one bccomics_setup runs.
```
Patches ordered in calculation time, from oldest(top) to newest(bottom)
-----------------------------------------------------------------------
Patch #  ix  iy  iz  Deltac/sigma(Deltac)  V_cb(km/s)  at z=1000
  1     118 102 120      1.000e+00          3.000e+01
  2      75  42 120      1.998e+00          5.003e+01
  3      40  84  79      2.966e-04          2.299e+01
Choose a patch of your interest; default is 3 if you just hit Enter below.
Enter your choice (patch #):
```
Assume you want Patch # 2 above, which is 2sigma Density peak with V_cb=50.03 km/s. Then enter 2:
```
Enter your choice (patch #):2
Patch # 2 chosen.
----- Interpolating transfer function -----
----- Convolving transfer function with random number -----
----- Calculating CDM position x -----
... (more messages)
*********** bccomics successfully ended ************
>>
```
After all runs are successuful, you will have a sub-directory under `ICdir`, named e.g. as 1.00Mpch_64_ic75_jc42_kc120. Under this direcotyr, you will have initial conditions.

### (4) For OCTAVE only, convert output binary ICs to enzo-readable HDF5 binary files.
Matlab can produce enzo initial conditions directly in step (3). But Octave still lacks the functionality to do so. The initial conditions from Octave will be files named "cpos1", "cpos2", ..., "vc1", ..., "db", "vb1", ..., "ether", "etot". Go to the initial condition directory, and run the provided python script:
```
$ cd /home/kjahn/BCCOMICS/ICs/1.00Mpch_64_ic75_jc42_kc120
$ python /home/kjahn/BCCOMICS/src_converter/convert2enzo.py
```
Then you will get enzo-readable HDF5 files named "ParticlePositions", "ParticleVelocities", "GridDensity", "GridVelocities", "GasThermalSpecEnergy", and "GasTotalSpecEnergy".

### (5) For overdense or underdense patches only, you absolutely should tweak the enzo parameter file.
Because cosmology simulation by enzo needs typically a periodic boundary condition, simulating overdense(underdense) patches requires fooling enzo or any other simulation codes (see section 3 of Ahn & Smith). For enzo, this is done by "enzo_patchcosmo.m". You first need to provide file `zglobal.dat` under the IC directory (e.g. `/home/kjahn/BCCOMICS/ICs/1.00Mpch_64_ic75_jc42_kc120`) listing "global redshifts" for enzo output dumps:
```
$ pwd
kjahn@machine:/home/kjahn/BCCOMICS/ICs/1.00Mpch_64_ic75_jc42_kc120
$ cat zglobal.dat
100
90
80
70
60
50
$
```
You may try zglobal.m in the source directory for easy creation of zglobal.dat to your taste.

Then, run octave(matlab) in this directory, and run "enzo_patchcosmo.m". It will then generate file `enzoparam_part.enzo`, which you should copy & paste to your final enzo parameter file (e.z. `localcosmo.enzo`). You'd better keep `zglobal.dat`, so that you can later match a "local redshift" in `enzoparam_part.enzo` to its corresponding "global redshift" in `zglobal.dat`. Take a look at `enzoparam_part.enzo`, and you will see non-zero curvature term, and also seemigly odd cosmological parameters and "local redshifts". But this is how you fool enzo (and similarly other simulation codes) to simulate overdense(underdense) patches.
```
$ octave
>> cd '/home/kjahn/BCCOMICS/ICs/1.00Mpch_64_ic75_jc42_kc120'
>> addpath('/home/kjahn/BCCOMICS/src')
>> enzo_patchcosmo
```
## To-Do
- Porting to Gadget, RAMSES, etc. WE WELCOME ANYONE WHO IS INTERESTED IN THIS EFFORT. DROP US A LINE, WE'LL ADD YOU AS A CONTRIBUTOR.
- Small-wavenumber tuning: For very small wavenumbers, the number of modes are too small to correctly reproduce average power spectrum. With unlucky gaussain random seed, this might cause unwated finte-box effect. May have to restrict randomness of power spectrum amplitude (not phase though) of smallest k modes.
- Nested grid: Not yet. Again, we welcome contributors on this side.

## Citing

If you use `BCCOMICS` in your work, please cite as:

```
...the BCCOMICS initial conditions generator (Ahn & Smith 2018).
```

BIBTEX:
```
@ARTICLE{0004-637X-869-1-76,
  author={Kyungjin Ahn and Britton D. Smith},
  title={Formation of First Galaxies inside Density Peaks and Voids under the Influence of Dark Matter–Baryon Streaming Velocity. I. Initial Condition and Simulation Scheme},
  journal={The Astrophysical Journal},
  volume={869},
  number={1},
  pages={76},
  url={http://stacks.iop.org/0004-637X/869/i=1/a=76},
  year={2018}
}
```
