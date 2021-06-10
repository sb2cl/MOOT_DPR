%% spMODEparam II
% Generates the required parameters to run the spMODE-II optimization 
% algorithm.
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
% This code implements a version of the multi-objective differential
% evolution algorithm with spherical pruning based on preferences 
% (spMODE-II, second version of the spMODE algorithm) described in:
%
% Gilberto Reynoso-Meza. Controller Tuning by Means of Evolutionary 
% Multiobjective Optimization: a Holistic Multiobjective Optimization 
% Design Procedure. PhD. Thesis (2014), Universitat Politècnica de 
% València. Url: http://hdl.handle.net/10251/38248. Supervisors: Javier 
% Sanchis and Xavier Blasco.
%
%%

dbstop if error
%% Variables regarding the multi-objective problem

spMODEDat.NOBJ = 3;                   % Number of decision objectives.

spMODEDat.NRES = 1;                   % Number of constraints plus
                                      % non-decision objectives. That is,
                                      % variables that are not used to
                                      % built the spherical grid, but to
                                      % calculate the physical index.
                                      % 1: number of oscillations
                                      % 2:sigma>antisigma

spMODEDat.NVAR = 7;                   % Number of decision variables

spMODEDat.FieldD  = [0.1 10;...   %4 pa RBS anti-sigma (1/min)
                   1 15;...  % CNa Copy number anti-sigma
                   0.1 20;...% ph RBS CHS enzyme (1/min)
                   1 15;...% CNh Copy number CHS enzyme
                   0.01  20; %kc binding rate sigma.Asigma complex (1/min)
                   0.01 10000;...%kd20 dissociation constant sigma-promoter (molecules)
                   log(2)/100 log(2)/30];   % Growth rate mu(1/min) [2h 45min]

spMODEDat.Initial = spMODEDat.FieldD; % Initialization bounds

spMODEDat.mop  =...
    str2func('CostFunction2_pathway_Anti');         % Cost Function

spMODEDat.CostProblem  = '';     % Problem Instance; In case you have 
                                      % several "subcases" for the 
                                      % "CostFunction"

%
%%
%% Variables regarding the optimization algorithm (Differential Evolution)

spMODEDat.Xpop = 50;                 % Population size

spMODEDat.SubXpop=50;                % SubPopulation size

spMODEDat.ScalingFactor = 0.5;       % Scaling factor

spMODEDat.CRrate= 0.5;               % Croosover Probability

spMODEDat.Recombination='binomial';  % binomial or lineal                                           
%
%%
%% Variables regarding convergence improving (elitism)

spMODEDat.CarElite = [];             % Solutions from the Approximated
                                     % Pareto front in a generation
                                     % to be merged with the population
                                     % in the evolution process. If empty,
                                     % a default value is used.
%%
%% Variables regarding spreading (spherical pruning)

spMODEDat.Strategy='SphP';           % 'Push' for a basic Dominance-based 
                                     % selection; 'SphP' for the spherical
                                     % pruning;
                                     
spMODEDat.Alphas=...                 % Number of Arcs (Strategy='SphP').   
    10*spMODEDat.NOBJ^2;             % Number of Arcs (Strategy='SphP'). alfa = 10.Nobj^(Nobj-1)

spMODEDat.Norm='physical';           % Norm to be used in Strategy='SphP';
                                     % It could be 'euclidean','manhattan',
                                     % 'infinite', 'physical' or 'custom'.
                                     % When using "custom" the user needs
                                     % to define his/her own custom
                                     % function to calculate the norm with
                                     % the format:
                                     % IndexesOUT=...
                                     %   CustomNorm(Front,Set,spMODEDat)
                                     
spMODEDat.PFrontSize=10*spMODEDat.NOBJ;             % Maximum Pareto optimal solutions required
%
%%
%% Variables regarding pertinency (Global Physical Programming)
% The following values for each design objective (decision objective +
% non-decision objective) and constraint shall be defined:
%
% Physical Matrix Definition.
% HD Highly Desirable
% D  Desirable
% T  Tolerable
% U  Intolerable
% HU Highly intolerable
%                           HD  -  D -  T  -  U  - HU
 spMODEDat.PhyMatrix{1} = [0.0    10  100  200   400;... %J1: g/L
                           0.0    10   25   30    40;... %J2:relative error %
                           0     0.5   0.7    1   10;... %J3:sigma oscillations  restricción 0     0.5   0.7    1   10;
                           0     0     0    0.5    1];   %2:sigma value
                       
                           %20    1200  1800 2000  2000;...%itae time min
% The above is based on _Global Physical Programming_ and _Physical
% Programming_; both are used to state preferences when dealing with 
% several objectives. For more see the following:

% J. Sanchis, M. Martínez, X. Blasco, G. Reynoso. Modelling preferences in
% multi-objective engineering design. Engineering Applications of 
% Artificial Intelligence. Vol. 23, num. 8, pp. 1255 - 1264, 2010.
%
% and
%
% A. Messac. Physical programming: effective optimization for computational
% design. AIAA Journal 34 (1), 149 – 158, 1996
                      
%%%spMODEDat.PhyIndexMax=...            % Tolerable vector is used as default.
    %PhyIndex(spMODEDat.PhyMatrix{1}(:,3)',spMODEDat);
spMODEDat.PhyIndexMax=PhyIndex([100, 30, 0.7, 0],spMODEDat);

   %spMODEDat.PhyIndexMax=PhyIndex([100, 30, 1700, 0.7, 0],spMODEDat);
% When using the tolerable vector as default, the algorithm will look to
% approximate a Pareto front approximation with only Tolerable values.

%%
%% Regarding Constraint Handling (different for objectives bound)
%
% Constraints could be defined as additional objectives. For an example
% please refer to:
%
% G. Reynoso-Meza, X. Blasco, J. Sanchis, M. Martínez. Multiobjective 
% optimization algorithm for solving constrained single objective problems.
% Evolutionary Computation (CEC), 2010 IEEE Congress on. 18-23 July 2010
%
% In such case, we encourage to define a Pertinency vector (above) or the
% PhyMatrix cell (also above) accordingly.
%
%%
%% Execution Variables

spMODEDat.MAXGEN =1e3;              % Generation bound

spMODEDat.MAXFUNEVALS = 1e4;        % Function evaluations bound

load 'xnominal.mat'
spMODEDat.PobInitial=xnominal;      % [] Initial population (if any) 

spMODEDat.SaveResults='yes';         % Write 'yes' if you want to 
                                    % save your results after the
                                    % optimization process;
                                    % otherwise, write 'no';

spMODEDat.Plotter='no';            % 'yes' if you want to see some
                                    % a graph at each generation.

spMODEDat.SeeProgress='yes';        % 'yes' if you want to see some
                                    % information at each generation.
                                    
spMODEDat.Plot_Cost_Func='yes';     %  'yes' if you want to plot all the variables
                                    % time courses in the cost function
%
%%
%% Put here the variables required by your code (if any).
% These variables could be used in your cost function file and/or your
% custom norm.
%
%
%
%% Release and bug report:
%
% June 2014: Initial release

