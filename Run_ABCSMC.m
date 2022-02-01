%% Run either ABC-SMC or model selection algorithms
% Both functions automatically save data as they go

testname='test1'; %Save file will be either "ABCSMC_testname_farm#" or "ModelSel_testname_farm#"
farm=1; % farms=1-5, correspond to IP1b, IP2a, IP3b, IP4b and IP7 respectively
Nrounds=18; %Number of rounds 
Nsuc=10000; %Number of particles to generate
Nw=4; %Number of workers

%% For running on a single machine 
%This will still use parellel computing unless parfor loops are replaced
%with for loops within the functions.

%Comment out the function which is not wanted
ABCSMC(testname,farm,Nrounds,Nsuc,Nw);
% Model_selection(testname,farm,Nrounds,Nsuc,Nw);

%% For submitting to a cluster for parallel computing

sched=parcluster('MatlabCluster');
job=createCommunicatingJob(sched);
set(job,'NumWorkersRange',[1 8]);

%Comment out the function which is not wanted
createTask(job,@ABCSMC,1,{{testname,farm,Nrounds,Nsuc,Nw}}); %{'name',farm,Nrounds,Nsuc,Nw}
% createTask(job,@Model_selection,1,{{testname,farm,Nrounds,Nsuc,Nw}});

submit(job);