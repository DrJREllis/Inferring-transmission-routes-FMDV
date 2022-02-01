%% Changeable parameters
Ntrials=10000; %No. of simulations
Ndays=40; %No. of days
dt=1/24; %Time interval for calculating shedding levels

dataset=1; %Parameter dataset - 0=priors, 1-5=posterior from farms 1-5, 6= no env. size, farm 1

if dataset>0
    %Overall farm data:
    Data_Infected=[38 44 36 54 14];
    Data_Pop=[38 49 47 54 16];
    Data_age=[10 7 5 10 5];

    %Individual farm data:
    Data_clin{1}=[1 0 13 8 10 1 3 2 0 0];
    Data_clin{2}=[1 8 7 5 8 13 2];
    Data_clin{3}=[3 2 11 17 3];
    Data_clin{4}=[5 3 12 13 19 2 0 0 0 0];
    Data_clin{5}=[2 4 4 3 1];
    
    Data_clin{6}=[1 0 13 8 10 1 3 2 0 0];
    if dataset<6
    N=Data_Pop(dataset);
    else N=38;
    end
end

% Uncomment the below to choose a different herd size
N=50; %No. of animals

if dataset==6
Esize=1;
else Esize=N;
end

Insp=1; %Days between inspections
Ninsp=3; %Number of animals inspected
Nsamp=5; %Number of samples on each part of the environment


beta=1; %Probability of infection through interaction
betaE=1; %Probability of infection through the environment



%% Load parameters
%ParVp=[Vp tp lg ld Tc sTc rho]; %ParTr=[beta betaE alpha delta];
[ParVp, ParTr, Ni0prior]=Load_parameters(dataset); 
ParTr(:,1)=ParTr(:,1)*beta; ParTr(:,2)=ParTr(:,2)*betaE;

%Initial Conditions
Ni0=Ni0prior; % %Initial infected
% Uncomment the below to set different initial conditions for infected animals and environmental contamination
% Ni0=1*ones(20000,1);
E0=[0 0 0 0]; %Initial contamination: [floor, wall, trough, faeces]

if dataset==6
    dataset=1;
end
%% Run model of transmission

k=1;
while k<=Ntrials    

    [VP, Ti, samp(k)]=Viral_profile(ParVp,N,Ndays,dt); %Create the viral profile of all animals
    
    [E(:,:,k),Cs(:,k),Cfirst(k),Ni(:,k),TotInf(:,k),dailyI(:,k),InfCause(:,:,k),...
        Infectiousness(:,:,k),ProbE(:,k),ProbD(:,k)]=Transmission_model(N,Ndays,dt,ParTr(samp(k),:),Ni0(samp(k)),E0,VP,Ti,Esize); %Run the transmission model

%     [Cs(:,k),Cfirst(k),Ni(:,k),TotInf(:,k),dailyI(:,k),InfCause(:,:,k),...
%         Infectiousness(:,:,k)]=Direct_model(N,Ndays,dt,ParTr(samp(k),1),Ni0(samp(k)),VP,Ti); %Run the transmission model

%     [E(:,:,k),Cs(:,k),Cfirst(k),Ni(:,k),dailyI(:,k),InfCause(:,:,k)]=...
%       Env_model(N,Ndays,dt,ParTr(samp(k),2:end),Ni0(samp(k)),E0,VP,Ti,Esize); %Run the transmission model

            
    if dataset>0 & Cfirst(k)<Ndays/dt-Data_age(dataset)/dt
        Cs1=Cs(1/dt:1/dt:Ndays/dt,k);
        if Cfirst(k)>1/dt
        Inffit(:,k)=Cs1(ceil(Cfirst(k)*dt)-1:ceil(Cfirst(k)*dt)-1+Data_age(dataset));
        else Inffit(:,k)=[0; Cs1(1:Data_age(dataset))];
        end
    end
    if Ni(end,k)>N/2
        k=k+1;  

    end
end

%Calculate the herd generation times
xp=0:dt:Ndays;
Tgd=trapz(xp,TotInf.*xp') ./ trapz(xp,TotInf) ;
Tge=trapz(xp,sum(E,2).*xp') ./ trapz(xp,sum(E,2)) ;

%Calculate Thetas for different detection times
Esum=permute(sum(E,2),[1,3,2])./0.56;
 for k=1:Ntrials
     if Cfirst(k)-1/dt<=1
        ThetaS(k,1)=0;
        ThetaE(k,1)=0;
     else 
        ThetaS(k,1)=trapz(dt,TotInf(1:Cfirst(k)-1/dt,k))./trapz(dt,TotInf(:,k));
        ThetaE(k,1)=trapz(dt,Esum(1:Cfirst(k)-1/dt,k))./trapz(dt,Esum(:,k));
     end
     ThetaS(k,2)=trapz(dt,TotInf(1:Cfirst(k),k))./trapz(dt,TotInf(:,k));
     ThetaE(k,2)=trapz(dt,Esum(1:Cfirst(k),k))./trapz(dt,Esum(:,k));
     ThetaS(k,3)=trapz(dt,TotInf(1:Cfirst(k)+1/dt,k))./trapz(dt,TotInf(:,k));
     ThetaE(k,3)=trapz(dt,Esum(1:Cfirst(k)+1/dt,k))./trapz(dt,Esum(:,k));
     ThetaS(k,4)=trapz(dt,TotInf(1:Cfirst(k)+2/dt,k))./trapz(dt,TotInf(:,k));
     ThetaE(k,4)=trapz(dt,Esum(1:Cfirst(k)+2/dt,k))./trapz(dt,Esum(:,k));
     ThetaS(k,5)=trapz(dt,TotInf(1:Cfirst(k)+5/dt,k))./trapz(dt,TotInf(:,k));
     ThetaE(k,5)=trapz(dt,Esum(1:Cfirst(k)+5/dt,k))./trapz(dt,Esum(:,k));
     ThetaS(k,6)=trapz(dt,TotInf(1:Cfirst(k)+10/dt,k))./trapz(dt,TotInf(:,k));
     ThetaE(k,6)=trapz(dt,Esum(1:Cfirst(k)+10/dt,k))./trapz(dt,Esum(:,k));
     
     ThetaS(k,7)=trapz(dt,TotInf(1:find(Cs(:,k)>0.05*N,1)+1/dt,k))./trapz(dt,TotInf(:,k));
     ThetaE(k,7)=trapz(dt,Esum(1:find(Cs(:,k)>0.05*N,1)+1/dt,k))./trapz(dt,Esum(:,k));     
     ThetaS(k,8)=trapz(dt,TotInf(1:find(Cs(:,k)>0.1*N,1)+1/dt,k))./trapz(dt,TotInf(:,k));
     ThetaE(k,8)=trapz(dt,Esum(1:find(Cs(:,k)>0.1*N,1)+1/dt,k))./trapz(dt,Esum(:,k));     
     ThetaS(k,9)=trapz(dt,TotInf(1:find(Cs(:,k)>0.2*N,1)+1/dt,k))./trapz(dt,TotInf(:,k));
     ThetaE(k,9)=trapz(dt,Esum(1:find(Cs(:,k)>0.2*N,1)+1/dt,k))./trapz(dt,Esum(:,k));    
     ThetaS(k,10)=trapz(dt,TotInf(1:find(Cs(:,k)>0.5*N,1)+1/dt,k))./trapz(dt,TotInf(:,k));
     ThetaE(k,10)=trapz(dt,Esum(1:find(Cs(:,k)>0.5*N,1)+1/dt,k))./trapz(dt,Esum(:,k));  
 end

toc 
 %% Figures - Generate plots such as Fig. 8

IC=permute(InfCause(:,1,:),[1 3 2]);
Ndirect=sum(IC==1)';
Nenv=sum(IC==2)';
Nboth=sum(IC==3)';
Inffit=[Inffit(1,:); Inffit(2:end,:)-Inffit(1:end-1,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Subplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

subplot(1,4,1)
xp=1:size(Inffit,1)-1; Inffit=Inffit(2:end,:);
bar(xp, fliplr(Data_clin{dataset}),'blue')
hold on
errorbar(xp,flipud(median(Inffit,2)),flipud(median(Inffit,2))-fliplr(prctile(Inffit',5))',...
    fliplr(prctile(Inffit',95))'-flipud(median(Inffit,2)),'.','color','red','MarkerSize',15,'linewidth',1)
ax = gca;
ax.FontSize = 12;
xlabel('Lesion age (days)','interpreter','latex','FontSize',14);
ylabel('Number of cattle','interpreter','latex','FontSize',14); 
ylim([0 20])
      

subplot(1,4,3)   
PDirect=Ndirect./(Ndirect+Nenv+Nboth);
PEnv=Nenv./(Ndirect+Nenv+Nboth);
G=[1 2]; x=[PDirect PEnv];
boxplot(x,'PlotStyle','compact','Symbol','','Orientation','horizontal',"ColorGroup",G,'colors',['b' 'r'])
 ylim([0.6 2.4])
t = text(-0.1,0.9,'Direct','HorizontalAlignment','center','interpreter','latex','FontSize',15);
set(t,'Rotation',90);
t = text(-0.1,1.9,'Environmental','HorizontalAlignment','center','interpreter','latex','FontSize',15);
set(t,'Rotation',90);
set(gca,'yticklabel',{})
 ax = gca;
ax.FontSize = 12;
 xlabel('Proportion infected','interpreter','latex','FontSize',14);


subplot(1,4,4)
xp=0:dt:Ndays;
yyaxis left
plot(xp,median(TotInf,2),'linewidth',2) % Total infectiousness
hold on
patch([xp, flip(xp,2)],...
[prctile(TotInf',95), flip(prctile(TotInf',5),2)],...
'b','linestyle','none','FaceVertexAlpha',0.5,...
'FaceAlpha',0.3);
ylabel('Viral shedding','interpreter','latex','FontSize',14);
xlim([0 40])

yyaxis right
xp=0:dt:40;
plot(xp,nanmedian(Esum,2),'color','r','linewidth',2)
patch([xp, flip(xp,2)],...
[prctile(Esum',95), flip(prctile(Esum',5),2)],...
'r','linestyle','none','FaceVertexAlpha',0.5,...
'FaceAlpha',0.3);
ylabel('Env. contamination','interpreter','latex','FontSize',14);
xlim([0 40])
ax = gca;
ax.FontSize = 12;
xlabel('Days since first infection','interpreter','latex','FontSize',14);


subplot(1,4,2)
plot(dt:dt:Ndays,median(Ni,2),'color','r') % Number infected
hold on
plot(dt:dt:Ndays,median(Cs,2)) % Number of animals with clinical signs
hold on
xp=dt:dt:Ndays;
patch([xp, flip(xp,2)],...
[prctile(Ni',95), flip(prctile(Ni',5),2)],...
'b','linestyle','none','FaceVertexAlpha',0.5,...
'FaceAlpha',0.5,'FaceColor','red');
xlim([0 10])
ylabel('Cumulative infected','interpreter','latex','FontSize',14);
ax = gca;
ax.FontSize = 12;
xlabel('Days since first infection','interpreter','latex','FontSize',14);
hold on
patch([xp, flip(xp,2)],...
[prctile(Cs',75), flip(prctile(Cs',25),2)],...
'r','linestyle','none','FaceVertexAlpha',0.5,...
'FaceAlpha',0.5);