function [ParVp, ParTr, Ni0prior]=Load_parameters(dataset)

if dataset==0
%% Load data for priors
varload=load('Data\SimplePhenomModel_VirusDynOnly_MCMCSamples');

Vp=[varload.ParSamp{1}(:,61:62); varload.ParSamp{2}(:,61:62)];
tp=[varload.ParSamp{1}(:,63:64); varload.ParSamp{2}(:,63:64)];
lg=[varload.ParSamp{1}(:,65:66); varload.ParSamp{2}(:,65:66)];
ld=[varload.ParSamp{1}(:,67:68); varload.ParSamp{2}(:,67:68)];

Tc=[varload.ParSamp{1}(:,69); varload.ParSamp{2}(:,69)];
sTc=[varload.ParSamp{1}(:,70); varload.ParSamp{2}(:,70)];
rho=[varload.ParSamp{1}(:,71); varload.ParSamp{2}(:,71)];

varload=load('Data\EnvTrans_SimpleModel21_MCMCSamples');

alpha=[exp(varload.ParSamp{3}(:,76:79)); exp(varload.ParSamp{3}(:,76:79))];
delta=[exp(varload.ParSamp{3}(:,80:83)); exp(varload.ParSamp{3}(:,80:83))];
betaE=[exp(varload.ParSamp{3}(:,84)); exp(varload.ParSamp{3}(:,84))];

varload=load('Data\SimpleWithinHostModelToo2_NF_MCMCSamples');

beta=exp([varload.ParSamp{1}(:,end-3); varload.ParSamp{2}(:,end-3)])./log(10);

ParVp=[Vp tp lg ld Tc sTc rho];
ParTr=[beta betaE alpha delta];

Ni0prior=ones(20000,1);

else
    if dataset<6
    %Load data from posterior distributions
    load('Data\PosteriorDistributions');
    Pars=Post(:,:,dataset);
    ParVp=Pars(:,1:11);
    ParTr=Pars(:,12:end);
    Ni0prior=Ni0post(:,1,dataset);
    else if dataset==6 %Load data for A=1 scenario for farm 1
        load('Data\Fitting_parameters_A1_farm1');
        Pars=[Theta{end}; Theta{end}];
        ParVp=Pars(:,1:11);
        ParTr=Pars(:,12:end);
        Ni0prior=[Ni0(:,end); Ni0(:,end)];
        end
    end
end
