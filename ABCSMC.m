function []=ABCSMC(testname,farm,Nrounds,Nsuc,Nw)
%% Initialisation

dt=1/24; %Time interval for calculating shedding levels

Model=1; %Transmission: 1-both, 2-direct only, 3-environ only

%Overall farm data:
Data_Infected=[38 44 36 54 14];
Data_Pop=[38 49 47 54 16];
Data_age=[10 7 5 10 5];
Age=Data_age(farm);

%Individual farm data:
Data_clin{1}=[1 0 13 8 10 1 3 2 0 0];
Data_clin{2}=[1 8 7 5 8 13 2];
Data_clin{3}=[3 2 11 17 3];
Data_clin{4}=[5 3 12 13 19 2 0 0 0 0];
Data_clin{5}=[2 4 4 3 1];


[ParVp, ParTr,~]=Load_parameters(0); %Load data for prior distributions 
Par=[ParVp ParTr]; %Combine parameters into one matrix
[~,dpar,pardis]=Prior_dist_fitting(Par,0,'NA','NA',Model); % Fits the priors to a distribution
eta=range(Par)*0.1; %length or std dev of Perturbation kernel

N=Data_Pop(farm); %Set herd size
Esize=N; %Relative size of environment - fix to herd size

eps=5*N*Age; %Set an initial (large) value for epsilon
    

%Set a maximum number of days for the simulation to run for - This comes
%into effect if no animals show clinical signs within this time
Maxdays=100;

tic
%% First round
parfor n=1:Nw % Split each round between workers - change to for loop if parellel computing is not available   
    
    % Initialize arrays
    Dp{n}=ones(ceil(Nsuc/Nw),1); Thetap{n}=[]; i=1; %Ddaily=0;
    wp{n}=ones(ceil(Nsuc/Nw),1); beta=zeros(1,size(Par,2)); 
    
    
    while i<=ceil(Nsuc/Nw) % The number of success from each worker

        beta=Par(randi(20000),:); %Randomly select from prior data
        
        % If uniform priors are chosen, then randomly generate from the given interval
        for p=1:size(Par,2)
            if pardis(p)=="uniform"
                beta(p)=rand*dpar{p}(1)+dpar{p}(2);
            end
        end
        
        Ni00=randi(N);
        
        [VP, Ti, ~]=Viral_profile(beta(1:11),N,Maxdays,dt); %Create the viral profile of all animals 
        if Maxdays>ceil(Ti(1,2))+Age+1
          %Set the number of days depending on how long for clinical signs to appear:
            Ndays=ceil(Ti(1,2))+Age+1; 
          %Run the transmission model:
            [~,Cs,Cfirst,~,~,~,~,~]=Transmission_model(N,Ndays,dt,beta(12:21),Ni00,0,VP,Ti,Esize); 

            Cs=Cs(1/dt:1/dt:Ndays/dt); %The number of animals showing clinical signs over time
         %Determine the number of clinical signs per day since the first so that it may ve compared to IP data:
            if Cfirst>0
                fitmet=Cs(ceil(Cfirst*dt):ceil(Cfirst*dt)+Age-1); 
                fitmet=[fitmet(1) fitmet(2:end)-fitmet(1:end-1)];
            else
                fitmet=0;
            end

            Ddaily=fitmet-Data_clin{farm}; %Daily difference between simulation and IP data
            Dp{n}(i)= sum(Ddaily.^2);
        else
      % If the time for clinical signs to appear in the first animal>max no. of days then parameter set fails
            break
        end

        if Dp{n}(i)<eps(1) %If the distance metric is less than epsilon
            Thetap{n}=[Thetap{n}; beta]; %Add the parameter set to Theta
            wp{n}(i)=1; %set the weight
            Ni0p{n}(i)=Ni00; %Initially infected
            i=i+1; 
        end
        
    end
end

% Concatenate cells from the different workers:
w1=cat(1,wp{:}); D1=cat(1,Dp{:}); Theta1=cat(1,Thetap{:}); Ni01=cat(1,Ni0p{:});
w(:,1)=w1(1:Nsuc); D(:,1)=D1(1:Nsuc); Theta{1}=Theta1(1:Nsuc,:); Ni0(:,1)=Ni01(1:Nsuc);
 
%% Further rounds
for t=2:Nrounds

roundtime(t-1)=toc %Keep track of time taken for each round
    
w(:,t-1)=w(:,t-1)/sum(w(:,t-1)); %normalise the weights
eps(t)=max(median(D(:,t-1)),0.1); % Set new epsilon to be median distance from the previous round
    

save(['Data\ABCSMC_',testname,'_farm',num2str(farm),'.mat'])
tic

parfor n=1:Nw % Split each round between  workers - change to for loop if parellel computing is not available   
    
    %Initialize arrays
    Dp{n}=ones(ceil(Nsuc/Nw),1); Thetap{n}=[]; i=1; 
    Cfirst{n}=zeros(length(Data_Infected),1); 
    wp{n}=ones(ceil(Nsuc/Nw),1); K=zeros(Nsuc,1); 
    Csp{n}={}; fitmet=zeros(Age,1); 


    while i<=ceil(Nsuc/Nw)
        
        %Generate perturbarions from the perturbation kernel
        pert=randn(1,size(Par,2)).*eta;
            
        %Sample from previous round with weights and add the perturbation  
        rsamp=randsample(Nsuc,1,true,w(:,t-1));
        beta=Theta{t-1}(rsamp,:)+pert;
        Ni00=Ni0(rsamp,t-1)+randi(5)-3;
        while Ni00<1 || Ni00>N
            Ni00=Ni0(rsamp,t-1)+randi(5)-3;
        end
        
        %Calculate the probability of the parameter set in the prior distributions
        [PriorProb,~,~]=Prior_dist_fitting(Par,beta,dpar,pardis,Model);
%         PriorProb=PriorProb*geopdf(Ni00,1/(4.1+0.0049*N));
        
        % Restrict parameters where applicable by generating a new perturbation
        while beta(11)<-1 || beta(11)>1 || min(beta([3,4,9,10,12:21]))<0 || PriorProb<=0
            pert=randn(1,size(Par,2)).*eta;
                
            %Sample from previous round with weights and add the perturbation  
            rsamp=randsample(Nsuc,1,true,w(:,t-1));
            beta=Theta{t-1}(rsamp,:)+pert;
            Ni00=Ni0(rsamp,t-1)+randi(5)-3;
            
            while Ni00<1 || Ni00>N
                Ni00=Ni0(rsamp,t-1)+randi(5)-3;
            end            
            %Calculate the probability of the parameter set in the prior distributions
            [PriorProb,~,~]=Prior_dist_fitting(Par,beta,dpar,pardis,Model);
        end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [VP, Ti, ~]=Viral_profile(beta(1:11),N,Maxdays,dt); %Create the viral profile of all animals 
        if Maxdays>ceil(Ti(1,2))+Age+1
          %Set the number of days depending on how long for clinical signs to appear:
            Ndays=ceil(Ti(1,2))+Age+1; 
          %Run the transmission model:
        [~,Csp{n}{i},Cfirst{n}(i),~,~,~,~,~,~]=Transmission_model(N,Ndays,dt,beta(12:21),Ni00,0,VP,Ti,Esize); 

            Cs=Csp{n}{i}(1/dt:1/dt:Ndays/dt); %The number of animals showing clinical signs over time
         %Determine the number of clinical signs per day since the first so that it may ve compared to IP data:
            if Cfirst{n}(i)>0
                fitmet=Cs(ceil(Cfirst{n}(i)*dt):ceil(Cfirst{n}(i)*dt)+Age-1);
                fitmet=[fitmet(1) fitmet(2:end)-fitmet(1:end-1)];
            else
                fitmet=0;
            end

            Ddaily=fitmet-Data_clin{farm}; %Daily difference between simulation and IP data
            Dp{n}(i)= sum(Ddaily.^2);
        else
      % If the time for clinical signs to appear in the first animal>max no. of days then parameter set fails
            break
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if Dp{n}(i)<eps(t) %If distance metric is less than epsilon
            Thetap{n}=[Thetap{n}; beta];
            Inffitp{n}(i,:)=fitmet;  Ctimep{n}(i,:)=Ti(:,2);
            Ni0p{n}(i)=Ni00; %Initially infected
            
            %Calculate K - the probability of the parameter set being chosen from each previous parameter set
            for j=1:Nsuc
                K(j)=prod(normpdf(Theta{t-1}(j,:),Thetap{n}(i,:),eta));
            end
            wp{n}(i)=PriorProb/(sum(w(:,t-1).*K)); %Generate the weight for the parameter set
            
            i=i+1;
        end
        
    end

end

% Concattenate from the seperate workers
w1=cat(1,wp{:}); D1=cat(1,Dp{:}); Theta1=cat(1,Thetap{:}); Ni01=cat(1,Ni0p{:});

w(:,t)=w1(1:Nsuc); D(:,t)=D1(1:Nsuc);
Ni0(:,t)=Ni01(1:Nsuc); Clin=cat(2,Csp{:}); Theta{t}=Theta1(1:Nsuc,:);
Inffit=cat(1,Inffitp{:}); Ctime=cat(1,Ctimep{:});
 Cf=cat(1,Cfirst{:}); 


clear Csp Thetap Inffitp wp Cfirst Ni0p
end
    
roundtime(t)=toc


save(['Data\ABCSMC_',testname,'_farm',num2str(farm),'.mat'])

