# Inferring-transmission-routes-FMDV
Code and data related to the paper 'Inferring transmission routes for foot-and-mouth disease virus within a cattle herd using approximate Bayesian computation'.

MATLAB REQUIREMENTS AND CODE:
The scripts/functions were run using Matlab version 2020b and require the Statistics and Machine Learning and Parallel Computing toolboxes. However, they can be easily adapted to run without the Parallel Computing toolbox by changing the "parfor" loop in the ParEst function to a "for" loop.

Within_Herd_Transmission_Model.m - Calls the Transmission_model function multiple times and plots transmission dynamics figures.
Transmission_model.m - Runs the transmission model with both direct and environmental transmission.
Direct_model.m - Runs the transmission model with direct transmission only.
Env_model.m -Runs the transmission model with environmental transmission only.

Load_parameters.m - Loads either the prior or posterior distributions for each parameter.
Viral_profile.m - Computes the viral profile for each animal.
Prior_dist_fitting.m - Computes prior distributions and calculates the prior probability of a particle.

Run_ABCSMC.m - Calls either the ABCSMC or Model_selection functions, either for parallel computing or on a single machine.
ABCSMC.m - Runs the ABC-SMC algorithm and saves the data after each round.
Model_selection.m - Runs the ABC-SMC for model selection algorithm and saves the data after each round.

Data folder:
SimpleWithinHostModelToo2_NF_MCMCSamples.mat
SimplePhenomModel_VirusDynOnly_MCMCSamples.mat
EnvTrans_SimpleModel21_MCMCSamples.mat
Posteriordistributions.mat
Posteriors_A1_farm1.mat
