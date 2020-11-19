This .ZIP file contains MATLAB scripts associated with the publication:

    Plourde A.P. and M. R. Nedimovic (Submitted). Earthquake depths, focal 
        mechanisms, and stress in the Lower St. Lawrence Seismic Zone.
        Seismological Research Letters. 

and was last updated November 18, 2020. 
 
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

dr0 = 'C:\Users\alex\NextCloud\MATLAB\';
addpath(dr0);
addpath([dr0,'EQ_FocalMech_MT\']);
addpath([dr0,'EQ_Plotting\']);
addpath([dr0,'hypoAPP\']);
addpath([dr0,'hypoTD\']);
addpath([dr0,'hypoTD\hypoTD_functions\']);
addpath([dr0,'hypoTD\hypoTD_pieces\']);
addpath([dr0,'Numerical_Methods\']);
addpath([dr0,'Signal_Processing\']);
clear dr0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    GENERAL INFO







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SELECTED FUNCTION INFO


% -- hypoAPP

This function estimates an earthquake hypocenter from P and S arrival times
and a 1D layered gradient velocity model. It minimizes an L1 misfit defined
as a weighted sum of absolute travel time residuals. The method uses iterative
grid-searches, each testing a 3D grid of potential hypocenters, centered around 
an initial estimate. It returns a best-fit hypocenter and 95% confidence 
intervals in each direction (E,N,Z). 

    
% --  hypoTD

Guo and Zhang (2017,GJI)

% -- fmAPP



% --  StressInvV2014

Vavrycuk (2014,GJI)


    





   




