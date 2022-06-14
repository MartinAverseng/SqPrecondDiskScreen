# SqPrecondDiskScreen
Matlab Code for "Quasi-local and frequency robust preconditioners for the Helmholtz 1st-kind BIE on the disk", by François Alouges and Martin Averseng.
To run, first get FMMLIB3D to compile on your machine with working mex files.

Then, to reproduce the results of the article, run the script "launchJobs.m" after setting the variable levels to the desired value. 
(levels 1-9 is the range used in the article, requires about 32GB of RAM).
For the results of the last subsection ("Different corrections"), launch "differentCorrections.m", again setting the levels and choosing Unif = 0 or 1
according to the desired type of mesh (Unif = 1: uniform mesh, Unif = 0: graded mesh). 



