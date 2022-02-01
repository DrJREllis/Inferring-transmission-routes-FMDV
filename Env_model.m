function [E,Cs,Cfirst,Ni,dailyI,InfCause]=Env_model(N,Ndays,dt,ParTr,Ni0,E0,VP,Ti,Esize)


%% Parameters from data
Nt=Ndays/dt; %Number of time steps

%Contamination parameters
betaE=ParTr(1,1);
alpha=ParTr(1,2:5);
delta=ParTr(1,6:9);

    %% Initialization
%Initialize vectors and matrices
Ni=zeros(Nt,1); V=zeros(N,Nt); State=zeros(N,Nt); E=zeros(Nt+1,4);
E(1,:)=E0; State(1:Ni0,1)=ones(1,Ni0); C=zeros(N,Nt); InfCause=zeros(N,2); 
dailyI=zeros(Nt,1); NewI=(1:Ni0); dailyI(1)=1; ProbE=NaN(Nt,1); 

%% Time loop
for j=1:Nt
    %Identify and count susceptible and infected animals
    S=find(State(:,j)==0); Ns=size(S,1);
    Infd=find(State(:,j)==1); Ni(j)=size(Infd,1);
    State(:,j+1)=State(:,j);


    %For newly infected animals, calculate their viral titre in each compartment and their 
    %shedding profile from day of infection to the end of the simulation.
    if dailyI(j)>0
        for i=NewI'
            V(i,j:Nt+1)=VP(i,1:Nt-j+2);
            C(i,j+ceil(Ti(i,2)/dt):Nt)=1;
        end
    end
    
    
    %Calculate environmental contamination - 0.56 accounts for unit conversion
    E(j+1,:)=E(j,:)- delta*dt.*E(j,:) + alpha.*dt*sum(max(0,log(0.56*V(:,j))))./Esize;
    
    % New infections    
    if Ni(j)<N %If there are still uninfected animals in the herd

        %Environmental transmission (homogeneous for each type of contamination)
        ProbE(j)=1-exp(-betaE.*dt*sum( (E(j+1,:)+E(j,:))./2, 2)); %Probability of transmission
        Einf=rand(Ns,1)<ProbE(j); % Generate new infections of susceptible animals
        
        % Record the route of infection and time 
        InfCause(S,1)=Einf; 
        InfCause(S,2)=j*dt;
        State(S,j+1)=State(S,j)+Einf; 
        
        NewI=find(State(:,j+1)-State(:,j)==1 );
        dailyI(j+1)=size(NewI,1);

    end

end

dailyI= sum(reshape(dailyI(1:Nt),1/dt,Nt*dt)); %Find how many cows are infected on each day

%Calculate the number of cows showing clinical signs and when they first occur
Cs=sum(C(:,:),1);
if sum(Cs)>0
    Cfirst=find(Cs>=1,1);
else 
    Cfirst=nan;
end




