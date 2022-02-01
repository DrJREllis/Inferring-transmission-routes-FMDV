# Inferring-transmission-routes-FMDV
Code and data related to the paper 'Inferring transmission routes for foot-and-mouth disease virus within a cattle herd using approximate Bayesian computation'.

MATLAB REQUIREMENTS AND CODE:
The scripts/functions were run using Matlab version 2020b and require the Statistics and Machine Learning and Parallel Computing toolboxes. However, they can be easily adapted to run without the Parallel Computing toolbox by changing the "parfor" loop in the ParEst function to a "for" loop.

Within_Herd_Transmission_Model.m - 
Transmission_model.m - 
Direct_model.m - 
Env_model.m -

Load_parameters.m - 
Viral_profile.m -
Prior_dist_fitting.m -

Run_ABCSMC.m -
ABCSMC.m -
Model_selection.m -
