function []=Model_selection(testname,farm,Nrounds,Nsuc,Nw)

%% Initialisation
dt=1/12; %Time interval for calculating shedding levels

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

[ParVp, ParTr,~]=Load_parameters(0); %Load data for prior distributions 
Par=[ParVp ParTr]; %Combine parameters into one matrix
[~,dpar,pardis]=Prior_dist_fitting(Par,0,'NA','NA',1); % Fits the priors to a distribution
eta=range(Par)*0.1; %length or std dev of Perturbation kernel

N=Data_Pop(farm); %Set herd size
Esize=N; %Relative size of environment - fix to herd size

eps=5*N*Data_age(farm); %Set an initial (large) value for epsilon
    

%Set a maximum number of days for the simulation to run for - This comes
%into effect if no animals show clinical signs within this time
Maxdays=100;

tic %Keep track of time taken for each round
%% First round
parfor n=1:Nw % Split each round between  workers - change to for loop if parellel computing is not available
    
    % Initialize arrays
    Dp{n}=ones(ceil(Nsuc/Nw),1); i=1; Ddaily=0;
    Cfirst=zeros(length(Data_Infected),1); ThModp{n}=cell(1,3);
    wp{n}=cell(1,3); beta=zeros(1,size(Par,2)); Ni0p{n}=cell(1,3);
    
    
    while i<=ceil(Nsuc/Nw) % The number of success from each worker

        beta=Par(randi(20000),:); %Randomly select from prior data
        Ni00=randi(N);
        
        Model=randi(3); %Randomly select a model from 1-3 uniformly
        % For single transmission route models, remove irrelevant parameters:
        if Model==2
            beta(13:size(Par,2))=NaN;
        else if Model==3
                beta(12)=NaN;
            end
        end
        
        
        [VP, Ti, ~]=Viral_profile(beta(1:11),N,Maxdays,dt); %Create the viral profile of all animals 
        if Maxdays>ceil(Ti(1,2))+Data_age(farm)+1
          %Set the number of days depending on how long for clinical signs to appear:
            Ndays=ceil(Ti(1,2))+Data_age(farm)+1; 
            
            %Run the chosen transmission model:
            if Model==1
                [~,Cs,Cfirst,~,~,~,~,~,~]=Transmission_model(N,Ndays,dt,beta(12:21),Ni00,0,VP,Ti,Esize); %Run the transmission model
            else if Model==2
                    [Cs,Cfirst,~,~,~,~,~]=Direct_model(N,Ndays,dt,beta(12),Ni00,VP,Ti);
                else if Model==3
                        [~,Cs,Cfirst,~,~,~]=Env_model(N,Ndays,dt,beta(13:21),Ni00,0,VP,Ti,Esize);
                    else print('Model number error')
                    end
                end
            end

            Cs=Cs(1/dt:1/dt:Ndays/dt); %The number of animals showing clinical signs over time
         %Determine the number of clinical signs per day since the first so that it may ve compared to IP data:
            if Cfirst>0
                fitmet=Cs(ceil(Cfirst*dt):ceil(Cfirst*dt)+Data_age(farm)-1); 
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
            ThModp{n}{Model}=[ThModp{n}{Model}; beta]; %Add the parameter set for that model to Theta
            wp{n}{Model}=[wp{n}{Model}; 1]; %set the weight
            Ni0p{n}{Model}=[Ni0p{n}{Model}; Ni00]; %Initially infected
            i=i+1; 
        end
        
    end
end

% Concatenate cells from the different workers:
w1=cat(1,wp{:}); D1=cat(1,Dp{:}); Ni01=cat(1,Ni0p{:});  ThMod1=cat(1,ThModp{:});

D(:,1)=D1(1:Nsuc); %Ni0(:,1)=Ni01(1:Nsuc); %w(:,1)=w1(1:Nsuc); ThMod{1}=ThMod1(1:Nsuc,:); Theta{1}=Theta1(1:Nsuc,:);
for k=1:3
    ThMod{1}{k}=cat(1,ThMod1{:,k});
    w{1}{k}=cat(1,w1{:,k});
    Ni0{1}{k}=cat(1,Ni01{:,k});
end

%% Further rounds
for t=2:Nrounds

roundtime(t-1)=toc %Keep track of time taken for each round
    
%Normalise the weights:
wsum=sum([w{t-1}{1}; w{t-1}{2}; w{t-1}{3}]);
for k=1:3
    w{t-1}{k}=w{t-1}{k}./wsum;
end
eps(t)=max(median(D(:,t-1)),0.1); % Set new epsilon to be median distance from the previous round
    

save(['Data\ModelSel_',testname,'_farm',num2str(farm),'.mat'])
tic

parfor n=1:Nw % Split each round between  workers - change to for loop if parellel computing is not available
    
    %Initialize arrays
    Dp{n}=ones(ceil(Nsuc/Nw),1); i=1; Ddaily=0; 
    Cfirst{n}=zeros(length(Data_Infected),1); Ndays=Maxdays;
    Csp{n}={};  fitmet=zeros(Data_age(farm),1);
     ThModp{n}=cell(1,3); wp{n}=cell(1,3);  Ni0p{n}=cell(1,3);
     K=ones(Nsuc,1);
    
    while i<=ceil(Nsuc/Nw)
        
        %Generate particle perturbations from the perturbation kernel
        pert=randn(1,size(Par,2)).*eta;
            
        %Sample from model previous round with weights and perturb
        Model=randsample(3,1,true,[sum(w{t-1}{1}) sum(w{t-1}{2}) sum(w{t-1}{3})]);
        Model=mod(Model+(rand<0.2)*randi(2)-1,3)+1;
        
        while  sum(w{t-1}{Model})==0 %If an 'extinct' model is chosen then try again
            Model=randsample(3,1,true,[sum(w{t-1}{1}) sum(w{t-1}{2}) sum(w{t-1}{3})]);
            Model=mod(Model+(rand<0.2)*randi(2)-1,3)+1;
        end
            
    %Sample particle for the selected model from previous round with weights and add the perturbation  
        rsamp=randsample(1:size(ThMod{t-1}{Model},1),1,true,w{t-1}{Model});
        beta=ThMod{t-1}{Model}(rsamp,:)+pert;
        Ni00=Ni0{t-1}{Model}(rsamp)+randi(5)-3; 
        while Ni00<1 || Ni00>N
            Ni00=Ni0{t-1}{Model}(rsamp)+randi(5)-3;
        end
        
        %Calculate the probability of the parameter set in the prior distributions
        [PriorProb,~,~]=Prior_dist_fitting(Par,beta,dpar,pardis,Model);
        
    % Restrict parameters where applicable by generating a new perturbation (redo the above)
        while beta(11)<-1 || beta(11)>1 || min(beta([3,4,9,10,12:21]))<0 || PriorProb<=0    
                pert=randn(1,size(Par,2)).*eta;
                
                Model=randsample(3,1,true,[sum(w{t-1}{1}) sum(w{t-1}{2}) sum(w{t-1}{3})]);
                Model=mod(Model+(rand<0.2)*randi(2)-1,3)+1;

                while  sum(w{t-1}{Model})==0
                    Model=randsample(3,1,true,[sum(w{t-1}{1}) sum(w{t-1}{2}) sum(w{t-1}{3})]);
                    Model=mod(Model+(rand<0.2)*randi(2)-1,3)+1;
                end

                rsamp=randsample(1:size(ThMod{t-1}{Model},1),1,true,w{t-1}{Model});
                beta=ThMod{t-1}{Model}(rsamp,:)+pert;
                Ni00=Ni0{t-1}{Model}(rsamp)+randi(5)-3; 
                while Ni00<1 || Ni00>N
                    Ni00=Ni0{t-1}{Model}(rsamp)+randi(5)-3;
                end

            [PriorProb,~,~]=Prior_dist_fitting(Par,beta,dpar,pardis,Model);
        end
        
        [VP, Ti, ~]=Viral_profile(beta(1:11),N,Maxdays,dt); %Create the viral profile of all animals  
        if Maxdays>ceil(Ti(1,2))+Data_age(farm)+1
            %Set the number of days depending on how long for clinical signs to appear:
            Ndays=ceil(Ti(1,2))+Data_age(farm)+1; 
            
            %Run the chosen transmission model:
            if Model==1
                [~,Csp{n}{i},Cfirst{n}(i),~,~,~,~,~,~]=Transmission_model(N,Ndays,dt,beta(12:21),Ni00,0,VP,Ti,Esize); 
            else if Model==2
                    [Csp{n}{i},Cfirst{n}(i),~,~,~,~,~]=Direct_model(N,Ndays,dt,beta(12),Ni00,VP,Ti);
                else if Model==3
                        [~,Csp{n}{i},Cfirst{n}(i),~,~,~]=Env_model(N,Ndays,dt,beta(13:21),Ni00,0,VP,Ti,Esize);
                    else print('Model number error')
                    end
                end
            end
            %The number of animals showing clinical signs over time
            Cs=Csp{n}{i}(1/dt:1/dt:Ndays/dt);
         %Determine the number of clinical signs per day since the first so that it may ve compared to IP data:
            if Cfirst{n}(i)>0
                fitmet=Cs(ceil(Cfirst{n}(i)*dt):ceil(Cfirst{n}(i)*dt)+Data_age(farm)-1);
                fitmet=[fitmet(1) fitmet(2:end)-fitmet(1:end-1)];
            else
                fitmet=0;
            end

            Ddaily=fitmet-Data_clin{farm};%Daily difference between simulation and IP data
            Dp{n}(i)= sum(Ddaily.^2);
        else
      % If the time for clinical signs to appear in the first animal>max no. of days then parameter set fails
            break
        end


        if Dp{n}(i)<eps(t) %If distance metric is less than epsilon
            ThModp{n}{Model}=[ThModp{n}{Model}; beta]; %Add the parameter set for that model to Theta
            Ni0p{n}{Model}=[Ni0p{n}{Model}; Ni00]; %Initially infected

            
            %Calculate K - the probability of the parameter set being chosen from each previous parameter set
            for j=1:size(ThMod{t-1}{Model},1)
                K(j)=prod(normpdf(ThMod{t-1}{Model}(j,:),beta,eta),'omitnan'); 
            end
            
            %Calculate 'S' from the ABC-SMC model selection algorithm
            S=(sum(w{t-1}{Model})*0.8+(1-sum(w{t-1}{Model}))*0.1)/sum(w{t-1}{Model});
            
            %Calculate the particle weight
            wp{n}{Model}=[wp{n}{Model}; (1/3)*PriorProb/(S*sum( w{t-1}{Model}.*K(1:size(w{t-1}{Model},1))))]; 
            i=i+1;
        end
        
    end

end

% Concattenate from the seperate workers
w1=cat(1,wp{:}); D1=cat(1,Dp{:}); Ni01=cat(1,Ni0p{:}); ThMod1=cat(1,ThModp{:});

 D(:,t)=D1(1:Nsuc);
 Clin=cat(2,Csp{:});  Cf=cat(1,Cfirst{:}); 

for k=1:3
    ThMod{t}{k}=cat(1,ThMod1{:,k});
    w{t}{k}=cat(1,w1{:,k});
    Ni0{t}{k}=cat(1,Ni01{:,k});
end


clear   Csp Inffitp  wp  Cfirst Ni0p ThMod1 ThModp
end

%Normalise the final weights
wsum=sum([w{t}{1}; w{t}{2}; w{t}{3}]);
for k=1:3
    w{t}{k}=w{t}{k}./wsum;
end

roundtime(t)=toc


save(['Data\ModelSel_',testname,'_farm',num2str(farm),'.mat'])
