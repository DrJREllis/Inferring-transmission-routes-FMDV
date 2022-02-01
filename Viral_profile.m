function [VP, Ti, samp]=Viral_profile(ParVp,N,Ndays,dt)

tau=0:dt:Ndays;

%Locate parameters for in-host dynamics
Vp=ParVp(:,1:2);
tp=ParVp(:,3:4);
lg=ParVp(:,5:6);
ld=ParVp(:,7:8);
Tc=ParVp(:,9:10);
rho=ParVp(:,11);

% Randomly generate a sample number 
samp=randi(size(Vp,1),1,1);


%% Calculate VP
Ti=exp(mvnrnd([tp(samp,1); Tc(samp,1)],[tp(samp,2)^2 rho(samp).*tp(samp,2).*Tc(samp,2); rho(samp).*tp(samp,2).*Tc(samp,2) Tc(samp,2)^2],N));
Vpi=exp(gamrnd(Vp(samp,1),Vp(samp,2)./Vp(samp,1),N,1));
gammadi=gamrnd(ld(samp,1),ld(samp,2)/ld(samp,1),N,1);    
gammagi=gamrnd(lg(samp,1),lg(samp,2)/lg(samp,1),N,1); 

%Calculate viral profile for each animal
VP=2.*Vpi./(exp(-gammagi.*(tau-Ti(:,1)))+exp(gammadi.*(tau-Ti(:,1))));

