This collection contains MATLAB scripts associated with the publication:

    Plourde A.P. and M. R. Nedimovic (Submitted). Earthquake depths, focal 
        mechanisms, and stress in the Lower St. Lawrence Seismic Zone.
        Seismological Research Letters. 

and was last updated February 3, 2021. 
 
Many of the MATLAB scripts have fairly detailed descriptions of their function
within themselves; however, this documentation is far from complete. Most of
the codes have been thoroughly tested, but I make zero guarantees that they 
function as intended. Testing was mostly performed on WINDOWS with 
MATLAB 2017a. There may be slight modifications required for operation with 
other versions of MATLAB or in UNIX or MAC environments, especially for codes 
that read or write external files.

    - Alexandre P. Plourde


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SET UP

To use these codes, I recommend adding the following lines to your MATLAB
startup.m file, modifying dr0 to your MATLAB home directory.

dr0 = 'YOUR_MATLAB_PATH\';
addpath(dr0);
addpath([dr0,'EQ_FocalMech_MT\']);
addpath([dr0,'EQ_Plotting\']);
addpath([dr0,'FocMechAPP\']);
addpath([dr0,'hypoAPP\']);
addpath([dr0,'hypoDD\']);
addpath([dr0,'hypoTD\']);
addpath([dr0,'Numerical_Methods\']);
addpath([dr0,'Signal_Processing\']);
clear dr0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    GENERAL INFO

Although many of the core functions in this package can be used indepently,
the "driver" files provided (COMING SOON) are designed to initially read in
earthquake hypocenter, station coordinates, and phase data from "hypoDD" format
event.dat, station.dat, and phase.dat files (see hypoDD manual --- 
Waldhauser [2001] USGS Open File Report) and save them in a MATLAB ".mat" file.
This file can then be reloaded to by "driver" scripts to compute waveform 
differential times, hypocenters, focal mechanisms, stress, etc., with the 
results of each analysis being appended to the ".mat" file. 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SELECTED FUNCTION INFO


% -- FocMechAPP

Estimates a focal mechanism and standard error (in degrees) using
two types of constraints P-wave first-motions and P-SV-SH absolute amplitude 
ratios. It requires a pre-defined uniform grid of focal mechanisms, which can 
be created with "createSDRgrid.m"---because this is a bit slow, I recommend 
saving the grid in a .mat file. The function computes a relative probability
for each grid-point focal mechanism, and finds the solution that minimizes the
standard error angle. See the drivers "FirstMotionAutopick.m" and 
"record_PSVSH.m" for guidance in creating input data from earthquake waveforms.


% -- hypoAPP

Estimates an earthquake hypocenter from P and S arrival times
and a 1D layered gradient velocity model. It minimizes an L1 misfit defined
as a weighted sum of absolute travel time residuals. The method uses iterative
grid-searches, each testing a 3D grid of potential hypocenters, centered around 
an initial estimate. It returns a best-fit hypocenter and 95% confidence 
intervals in each direction (E,N,Z). 


% --  hypoDD

An implementation of the double-difference earthquake relocation method
"hypoDD" by Waldhauser and Ellsworth (2000,BSSA). It also includes functions
to read/write files related to running the original hypoDD software.
 
 
% --  hypoTD

Estimates earthquake hypocenters with pick-based and waveform-based
differential travel-times, based on the "triple diffrence" or 
"double-pair double-difference" of Guo and Zhang (2017,GJI). 


% --  StressInvV2014

Estimates the tectonic stress tensor from focal mechanisms. Standard error 
angles on each focal mechanism can be input, allowing robust uncertainty 
estimates of the stress tensor. The name owes to its inspiration from the 
algorithm of Vavrycuk (2014,GJI).

%%%%%%%%%%%%%%%%%%%%%%%%%%
References

Guo, H., and H. Zhang (2017). Development of double-pair double dierence earthquake location
algorithm for improving earthquake locations. Geophysical Journal International, 208, 333-348.
doi:10.1093/gji/ggw397.

Vavrycuk V. [2014]. Iterative joint inversion for stress and fault orientations from focal
mechanisms. Geophysical Journal International, 199, 69-77, doi:10.1093/gji/ggu224.

Waldhauser F and Ellsworth WL. [2000]. A double-difference earthquake location algorithm:
Method and application to the Northern Hayward Fault, California. Bulletin of the
Seismological Society of America, 90, 1353-1368, doi:10.1785/0120000006.

    





   




