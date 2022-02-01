function [Cs,Cfirst,Ni,TotInf,dailyI,InfCause,Infectiousness]=Direct_model(N,Ndays,dt,ParTr,Ni0,VP,Ti)

%% Parameters from data
Nt=Ndays/dt; %Number of time steps

%Contamination parameters
beta=ParTr(1,1);

    %% Initialization
%Initialize vectors and matrices
Ni=zeros(Nt,1); V=zeros(N,Nt); Infectiousness=zeros(N,Nt+1); State=zeros(N,Nt); 
State(1:Ni0,1)=ones(1,Ni0); C=zeros(N,Nt); InfCause=zeros(N,2); 
dailyI=zeros(Nt,1); NewI=(1:Ni0); dailyI(1)=1; 
ProbD=NaN(Nt,1);

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
        %Infectiousness is assumed to be proportional to log titre
        Infectiousness(NewI,:)=max(0,log(V(NewI,:)));
    end
    
    
    % New infections    
    if Ni(j)<N %If there are still uninfected animals in the herd
        
        %Direct transmission
        ProbD(j)=1-exp(-beta.*dt.*(mean((Infectiousness(:,j+1)+Infectiousness(:,j))./2) )); %Probability of transmission
        Dinf=rand(Ns,1)<ProbD(j); % Generate new infections of susceptible animals
         
        % Record the route of infection and time 
        InfCause(S,1)=Dinf; 
        InfCause(S,2)=j*dt;
        State(S,j+1)=State(S,j)+Dinf; 
        
        NewI=find(State(:,j+1)-State(:,j)==1 );
        dailyI(j+1)=size(NewI,1);

    end


end

dailyI= sum(reshape(dailyI(1:Nt),1/dt,Nt*dt)); %Find how many cows are infected on each day
TotInf=sum(Infectiousness(:,:)); %Calculate total infectiousness

%Calculate the number of cows showing clinical signs and when they first occur
Cs=sum(C(:,:),1);
if sum(Cs)>0
    Cfirst=find(Cs>=1,1);
else 
    Cfirst=nan;
end


