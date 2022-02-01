function [PriorProb,dpar,pardis]=priordistfitting(Par,beta,dpar,pardis,Model)

%Model input: 1-combined model, 2-direct only transmission, 3-env only transmission

%If distributions have not yet been assigned:
if size(pardis,2)~=size(Par,2)
    clear pardis
% Fitted prior distributions - can choose between normal, gev or gamma for best fit.
% These prior distributions are used to generate the prior probabilities to calculate particle weights
% and are only meant to capture the shape of the distribution of the prior data.
% Change to "uniform" for uninformative priors.
pardis(1)="gev";
pardis(2)="gamma";
pardis(3)="normal";
pardis(4)="gev";
pardis(5)="gev";
pardis(6)="gev";
pardis(7)="gev";
pardis(8)="kernel";
pardis(9)="normal";
pardis(10)="gev";
pardis(11)="kernel";
pardis(12)="gamma";
pardis(13)="gamma";
pardis(14)="gamma";
pardis(15)="gamma";
pardis(16)="gev";
pardis(17)="gamma";
pardis(18)="normal";
pardis(19)="normal";
pardis(20)="normal";
pardis(21)="normal";
end

% For single transmission route models, remove irrelevant parameters
if Model==2
    pardis(13)="zero";
else if Model==3
    pardis(12)="zero";
    end
end

priorprob=size(Par,2);
PriorProb=0;

if beta==0 %If a particle isn't given then save the prior distributions 
    clear dpar
    for i=1:size(Par,2)

        parmin=min(Par(:,i));
        parmax=max(Par(:,i));
        parint=(parmax-parmin)/100;


        if pardis(i)=="gamma"
                dpar{i}=gamfit(Par(:,i));          
                
            else if pardis(i)=="gev" 
                dpar{i}=gevfit(Par(:,i));

            else if pardis(i)=="normal" 
                [a,b]=normfit(Par(:,i));   
                dpar{i}=[a,b];
                
            else if pardis(i)=="kernel"
                    dpar{i} = fitdist(Par(:,i),'Kernel','Kernel','normal');

            else if pardis(i)=="uniform"
                    dpar{i}=[parmax-parmin, parmin];
            end
            end
            end
            end
        end
        if Model==2
            dpar{13}=0;
        else if Model==3
            dpar{12}=0;
            end
        end
    end  
else %Otherwise calculate the prior probability of selecting the given particle:
    for ii=1:size(beta,2)
        if pardis(ii)=="gamma"
            priorprob(ii)=gampdf(beta(ii),dpar{ii}(1),dpar{ii}(2)); 
            else if pardis(ii)=="normal"
                priorprob(ii)=normpdf(beta(ii),dpar{ii}(1),dpar{ii}(2)); 
            else if pardis(ii)=="gev" 
                priorprob(ii)=gevpdf(beta(ii),dpar{ii}(1),dpar{ii}(2),dpar{ii}(3));    
            else if pardis(ii)=="kernel"
                priorprob(ii)=pdf(dpar{ii},beta(ii));  
            else if pardis(ii)=="uniform"
                    if beta(ii)>dpar{ii}(2) && beta(ii)<=sum(dpar{ii})
                        priorprob(ii)=1/dpar{ii}(1);
                    else priorprob(ii)=0;
                    end
            else if pardis(ii)=="zero"
                priorprob(ii)=1;
            end
            end
            end
            end
            end
        end
    end
    
    if Model==2
        priorprob(13:end)=1;
    else if Model==3
            priorprob(12)=1;
        end
    end
    PriorProb=prod(priorprob);
end