%% SphPruning
%
% Prune a set of multi-objective vectors using spherical coordinates.
%
%%
%%
% J  [OUT] : The objective Vector. J is a matrix with as many rows as
%            trial vectors in X and as many columns as objectives.
% X   [IN] : Decision Variable Vector. X is a matrix with as many rows as
%            trial vector and as many columns as decision variables.
% Dat [IN] : Parameters defined in spMODEparam.m and updated in spMODE.m
%
% PFront [OUT] : The Pruned Front. PFront is a matrix with as many rows
%                as vectors in PSet and as many columns as objectives.
% PSet   [OUT] : The Pruned Set. PSet is a matrix with as many rows as
%                vectors and as many columns as decision variables.
% Dat    [OUT] : Updated parameters
%
%%
%% Beta version 
% Copyright 2006 - 2014 - CPOH  
%
% Predictive Control and Heuristic Optimization Research Group
%      http://cpoh.upv.es
%
% ai2 Institute
%      http://www.ai2.upv.es
%
% Universitat Politècnica de València - Spain.
%      http://www.upv.es
%
%%
%% Author
% Gilberto Reynoso-Meza
% http://www.researchgate.net/profile/Gilberto_Reynoso-Meza
% http://www.mathworks.es/matlabcentral/fileexchange/authors/289050
%%
%% For new releases and bug fixing of this Tool Set please visit:
% 1) http://cpoh.upv.es/en/research/software.html
% 2) Matlab Central File Exchange
%%
%% Overall Description
% This code implements the spherical pruning described in:
%
% Gilberto Reynoso-Meza, Javier Sanchis, Xavier Blasco, Miguel Martínez.
% Design of Continuous Controllers Using a Multiobjective Differential
% Evolution Algorithm with Spherical Pruning. Applications of Evolutionary 
% Computation. LNCS Volume 6024, 2010, pp 532-541.
% 
% It can be used with any set of solutions. Dominance is guaranted IIF
% the input set has only Pareto Optimal solutions.
%
function [PFront PSet Dat]=SphPruning(J,X,Dat)

% Reading variables from Dat
Nobj   = Dat.NOBJ;     % Number of objectives.
Nvar   = Dat.NVAR;     % Number of variables.
Nres   = Dat.NRES;     % Number of costraints.
Xpop   = (size(J,1));  % Population size.
Alphas = Dat.Alphas;   % Number of arcs in the grid.

% Updating variables
UpperLim                = max([Dat.AnchorsY;J]); % Look for extremes.
LowerLim                = min(Dat.AnchorsY); 

Dat.ExtremesPFront(1,:) = max([Dat.AnchorsY;J]); % Update Extremes
Dat.ExtremesPFront(2,:) = min(Dat.AnchorsY);

% Initalization
Normas = zeros(Xpop,1);

% To be sure we have numbers to work on.
for n=1:size(UpperLim,2)
    if UpperLim(1,n)==LowerLim(1,n)
        UpperLim(1,n)=UpperLim(1,n)+1;
    end
    if isnan(UpperLim(1,n)) || isnan(LowerLim(1,n))
        disp('alarma NAN');
        pause;
    elseif isinf(UpperLim(1,n)) || isinf(LowerLim(1,n))
        disp('alarma INF');
        pause;
    end
    
end

    %% CALCULATING SPHERICAL COORDINATES

Arcs = HypCar2HypSph(J(:,1:Nobj),Dat.ExtremesPFront(2,1:Nobj),Dat.ExtremesPFront(1,1:Nobj));

if size(Arcs,1)==1
    MaxArc=Arcs+1;
    MinArc=Arcs;
else
    MaxArc=max(Arcs);
    MinArc=min(Arcs);
end

for i=1:size(Arcs,2)
    if MaxArc(1,i)==MinArc(1,i)
        MaxArc(1,i)=MaxArc(1,i)+1;
    end
end

    %% DEFINING THE SIGHT RANGE
    
Sight = (MaxArc-MinArc)./Alphas;
Sight = (1./Sight);

    %% ASSIGNING SPHERICAL SECTORS TO POPULATION
    
for xpop=1:Xpop
    Arcs(xpop,:)=ceil(Sight(1,:).*Arcs(xpop,:));
end

    %% COMPUTING NORMS
    
for xpop=1:Xpop
    if strcmp(Dat.Norm,'euclidean')
        Normas(xpop,1) = ...
            norm((J(xpop,1:Nobj)-Dat.ExtremesPFront(2,1:Nobj)) ./ ...
            (Dat.ExtremesPFront(1,1:Nobj)-Dat.ExtremesPFront(2,1:Nobj)),2);
    elseif strcmp(Dat.Norm,'manhattan')
        Normas(xpop,1) = ...
            norm((J(xpop,1:Nobj)-Dat.ExtremesPFront(2,1:Nobj)) ./ ...
            (Dat.ExtremesPFront(1,1:Nobj)-Dat.ExtremesPFront(2,1:Nobj)),1);
    elseif strcmp(Dat.Norm,'infinite')
        Normas(xpop,1) = ...
            norm((J(xpop,1:Nobj)-Dat.ExtremesPFront(2,1:Nobj)) ./ ...
            (Dat.ExtremesPFront(1,1:Nobj) - ...
            Dat.ExtremesPFront(2,1:Nobj)),Inf);
    elseif strcmp(Dat.Norm,'physical')
        Normas(xpop,1)=PhyIndex(J(xpop,:),Dat);
    elseif strcmp(Dat.Norm,'custom')
        Normas(xpop,1)=CustomNorm(J(xpop,:),X(xpop,:),Dat);
    end
end


    %% PRUNING
    
k = 0;
Dominancia = zeros(Xpop,1);
PFront     = zeros(Xpop,Nobj+Nres);
PSet       = zeros(Xpop,Nvar);
Arcos      = zeros(Xpop,Nobj-1);


for xpop=1:Xpop
    Dominado=Dominancia(xpop,1);

    if strcmp(Dat.Norm,'physical')
       if Normas(xpop,1)>Dat.PhyIndexMax
           Dominado=1;
           Dominancia(xpop,1)=1;
       end
    end
    
    if Dominado==0    
        for compara=1:Xpop
            if (Arcs(xpop,:)==Arcs(compara,:) )
                if (xpop~=compara)
                    if Normas(xpop,1)>Normas(compara,1)
                        Dominancia(xpop,1)=1;
                        break;
                    elseif Normas(xpop,1)<Normas(compara,1)
                        Dominancia(compara,1)=1;
                    elseif Normas(xpop,1)==Normas(compara,1)
                        if compara>xpop
                            Dominancia(compara,1)=1;
                        end
                    end
                end
            end
        end
    end
    
    if Dominancia(xpop,1)==0
        k=k+1;
        PFront(k,:)=J(xpop,:);
        PSet(k,:)=X(xpop,:);
        Arcos(k,:)=Arcs(xpop,:);
    end
end

if k==0
    if strcmp(Dat.Norm,'physical')
         Acomoda=sortrows([Normas,J,X]);
         PFront=Acomoda(1,2:Nobj+Nres+1);
         PSet=Acomoda(1,Nobj+Nres+2:end);
    else
        PFront=J;
        PSet=X;
    end
elseif k>Dat.PFrontSize
         Acomoda=sortrows([Normas,J,X]);
         PFront=Acomoda(1:Dat.PFrontSize,2:Nobj+Nres+1);
         PSet=Acomoda(1:Dat.PFrontSize,Nobj+Nres+2:end);    
else
    PFront=PFront(1:k,:);
    PSet=PSet(1:k,:);
end


%%

%% Computing Spherical Coordinates
function Arcs=HypCar2HypSph(X,LowerLim,UpperLim)

Xpop  = size(X,1);
Nobj  = size(X,2);
Narcs = Nobj-1;
Arcs  = zeros(Xpop,Narcs);

for xpop = 1:Xpop
    X(xpop,:)=1+(X(xpop,:)-LowerLim)./(UpperLim-LowerLim);
    Arcs(xpop,1)=atan2(X(xpop,Nobj),X(xpop,Nobj-1));
    
    if Narcs>1
        for narcs=2:Narcs
            Arcs(xpop,narcs) = ...
                atan2(norm(X(xpop,(Nobj-narcs+1):Nobj)), ...
                X(xpop,Nobj-narcs));
        end
    end
    
    for narcs=1:Narcs
        if isnan(Arcs(xpop,narcs))==1
            Arcs(xpop,narcs)=0;
        end
        if isinf(Arcs(xpop,narcs))==1
            Arcs(xpop,narcs)=0;
        end
    end
end

Arcs=(360/(2*pi))*(Arcs);

%%
%% Release and bug report:
%
% June 2014: Initial release