%% Multi-objective Differential Evolution Algorithm with Spherical Pruning based on preferences (spMODE-II) - Beta version -
%% 
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
%
% http://www.researchgate.net/profile/Gilberto_Reynoso-Meza
%
%%
%% For new releases and bug fixing of this Tool Set please visit:
%
% 1) http://cpoh.upv.es/en/research/software.html
%
% 2) Matlab Central File Exchange
%
%%
%% Tutorial Description
%
% In this tutorial, basic problems are solved using the spMODE-II algorithm,
% which is a version of the multi-objective differential evolution algorithm
% with spherical pruning based on preferences described in:
%
% *G. Reynoso-Meza*. _Controller Tuning by Means of Evolutionary 
% Multiobjective Optimization: a Holistic Multiobjective Optimization 
% Design Procedure._PhD. Thesis (2014), Universitat Politècnica de 
% València. Url: http://hdl.handle.net/10251/38248. Supervisors: Javier 
% Sanchis and Xavier Blasco.
% 
% Basic features of the algorithm are:
%
% # Improving Convergence by using an external file to store pertinent 
% solutions and include them in the evolutionary process.
% # Improving Spreading by using a spherical pruning mechanism based on 
% preferences.
% # Improving Pertinency of solutions with a preferences handling mechanism
% based on global physical programming. This features enables the algorithm
% to deal efficiently with constrained and many-objective optimization
% instances.
% # Size control mechanism to retain a fixed amount of solutions (the most
% preferable solutions according to the set of preferences).
%
%% 
%% Scripts and functions listing
% # RunTutorial.m  - Runs the tutorial.
% # Tutorial.m     - The Tutorial script.
% # spMODEparam.m  - Script to built the struct required for the
% optimization.
% # spMODE.m       - The optimization algorithm.
% # SphPruning.m   - The spherial pruning mechanism.
% # CostFunction.m - Cost function definition.
%
%%
%% Basic Example
%
% Run the spMODEparam file to build the variable "spMODEDat" with the
% variables required for the optimization.

spMODEparam;

spMODEDat

%%
% All of them are self-explainend in the spMODEparam.m script. Some default
% values will be modified in order to improve clarity on the features of
% this algorithm.
%
% The problem to solve is a known benchmark problem (DTLZ2). Now, run the
% algorithm:
%% BI-OBJECTIVE BASE CASE EXAMPLE: execution with the original sp-MODE algorithm.

disp('Running base case example:')

clear spMODEDat

% Default parameters are used
spMODEparam;

% Now we define a clasical p-norm (preferences are not included) fromt the
% original sp-MODE algorithm.
spMODEDat.Norm='euclidean';

% We modify the default values of the number of arcs in the spherical grid
% and the quantity of solutions in the approximated Pareto front that we
% are seeking (to improve clarity).
spMODEDat.Alphas=100;
spMODEDat.PFrontSize=100;

% Execute the algorithm
OUT=spMODEII(spMODEDat)

% Ploting an usual run without preferences (base case).
plot(OUT.PFront(:,1),OUT.PFront(:,2),'xb'); grid on;


%% FIRST EXAMPLE: basic set of preferences.

clear spMODEDat

% Default parameters are used
spMODEparam;

% We use, instead a p-norm, the physical programming index.

spMODEDat.Norm='physical';

% Now we define our set of preferences (mandatory) in the format:
% J-ith objective : HD  -  D -  T  -  U  - HD
% HD Highly Desirable
% D  Desirable
% T  Tolerable
% U  Untolerable
% HU Highly Untolerable

 spMODEDat.PhyMatrix{1} = [0.0  0.20  0.70  0.8  2.0   10;
                           0.0  0.60  0.90  1.2  2.0   10];   
                       
% We modify the default values of the number of arcs in the spherical grid
% and the quantity of solutions in the approximated Pareto front that we
% are seeking (to improve clarity).

spMODEDat.Alphas=100;
spMODEDat.PFrontSize=100;

% Execute the algorithm

OUT1=spMODEII(spMODEDat)

% Ploting the Tolerable, Desirable and Highly Desirable vectors.

figure;
plot(0.8,1.2,'ob'); hold on; grid on;
plot(0.7,0.9,'ob');
plot(0.2,0.6,'ob');

% Ploting an usual run without preferences (base case).
plot(OUT.PFront(:,1),OUT.PFront(:,2),'xb');
% Ploting an usual run taking into account preferences.
plot(OUT1.PFront(:,1),OUT1.PFront(:,2),'dr');

%% SECOND EXAMPLE: size control mechanism.

clear spMODEDat

% Default parameters are used
spMODEparam;

% The following means we are looking just for 30 solutions in the
% pertinent Pareto front approximation; Solutions with better physical
% index will be prefered.

spMODEDat.PFrontSize=30;

% Now we define our set of preferences (mandatory) in the format:
% J-ith objective : HD  -  D -  T  -  U  - HD
% HD Highly Desirable
% D  Desirable
% T  Tolerable
% U  Untolerable
% HU Highly Untolerable

 spMODEDat.PhyMatrix{1} = [0.0  0.20  0.70  0.8  2.0   10;
                           0.0  0.60  0.90  1.2  2.0   10];   

% We modify the default values of the number of arcs in the spherical grid
% (to improve clarity).

spMODEDat.Alphas=100;                       

% Execute the algorithm

OUT2=spMODEII(spMODEDat)

% Ploting the Tolerable, Desirable and Highly Desirable vectors.
figure;
plot(0.8,1.2,'ob'); hold on; grid on;
plot(0.7,0.9,'ob');
plot(0.2,0.6,'ob');
% Ploting an usual run without preferences (base case).
plot(OUT.PFront(:,1),OUT.PFront(:,2),'xb');
% Ploting an usual run taking into account preferences ans size control. 
% Solutions with higher GPP will be discarded (less preferable solutions).
plot(OUT2.PFront(:,1),OUT2.PFront(:,2),'dr');

%% THIRD EXAMPLE: two set of preferences.

clear spMODEDat

% Default parameters are used
spMODEparam;

% The following means that two preference conditions are stated.

spMODEDat.PhyMatrix{1} = [0.0  0.5  1.0  1.2  2.0 10;
                          0.0  0.4  0.5  1.2  2.0 10];

spMODEDat.PhyMatrix{2} = [0.0  0.4  0.5  1.2  2.0 10;
                          0.0  0.5  1.0  1.2  2.0 10];  


% The following means that solutions with more than 2 objectives values in 
% the tolerable region are not allowed. It is valid for all instances.

spMODEDat.PhyIndexMax=PhyIndex([0.5 1.0],spMODEDat);                  

% or

spMODEDat.PhyIndexMax=PhyIndex([1.0 0.5],spMODEDat);                  
% In both cases, the numerical value for "Maximum Tolerability" is the same

% We modify the default values of the number of arcs in the spherical grid
% and the quantity of solutions in the approximated Pareto front that we
% are seeking (to improve clarity).

spMODEDat.Alphas=100;
spMODEDat.PFrontSize=100;

% Execute the algorithm
OUT3=spMODEII(spMODEDat)

% Ploting the Tolerable, Desirable and Highly Desirable vectors for both 
% instances.
figure;

plot(1.2,1.2,'ob'); hold on; grid on;
plot(1.0,0.5,'ob');
plot(0.5,0.4,'ob');

plot(1.2,1.2,'sb'); hold on; grid on;
plot(0.5,1.0,'sb');
plot(0.4,0.5,'sb');

% Ploting an usual run without preferences (base case).
plot(OUT.PFront(:,1),OUT.PFront(:,2),'xb');
% Ploting an usual run taking into account preferences and size control. 
% Solutions with higher GPP will be discarded (less preferable solutions).
plot(OUT3.PFront(:,1),OUT3.PFront(:,2),'dr');

%% Tri-OBJECTIVE BASE CASE EXAMPLE: execution with the original sp-MODE algorithm.

disp('Running base case example:')

clear spMODEDat

% Default parameters are used
spMODEparam;

% We define a tri-objective problem
spMODEDat.NOBJ = 3;

% Now we define a clasical p-norm (preferences are not included) fromt the
% original sp-MODE algorithm.
spMODEDat.Norm='euclidean';

% We modify the default values of the number of arcs in the spherical grid
% and the quantity of solutions in the approximated Pareto front that we
% are seeking (to improve clarity). Remember, if m is the number of
% objectives, and alpha the number of partitions in each dimension, that 
% means that the spherical grid have m^(alpha-1) sectors.

spMODEDat.Alphas=15;
spMODEDat.PFrontSize=225;

% Execute the algorithm
OUT=spMODEII(spMODEDat)

% Ploting an usual run without preferences (base case).
figure;
plot3(OUT.PFront(:,1),OUT.PFront(:,2),OUT.PFront(:,3),'xb'); grid on;

%% FOURTH EXAMPLE: Maximum tolerability mechanism

clear spMODEDat

% Default parameters are used
spMODEparam;

% We define a tri-objective problem
spMODEDat.NOBJ = 3;

% Now we define our set of preferences (mandatory) in the format:
% J-ith objective : HD  -  D -  T  -  U  - HD
% HD Highly Desirable
% D  Desirable
% T  Tolerable
% U  Untolerable
% HU Highly Untolerable

spMODEDat.PhyMatrix{1} = [0.0  0.05  0.1  0.4  1.0   10;
                          0.0  0.30  0.4  0.6  1.0   10;
                          0.0  0.50  0.8  0.9  1.0   10]; 

% The following means that solutions with more than 2 objectives values in 
% the tolerable region are not allowed:

spMODEDat.PhyIndexMax=PhyIndex([0.4, 0.6, 0.8],spMODEDat);

% you could also use:

spMODEDat.PhyIndexMax=PhyIndex([0.4, 0.4, 0.9],spMODEDat);

% or

spMODEDat.PhyIndexMax=PhyIndex([0.1, 0.6, 0.9],spMODEDat);

% We modify the default values of the number of arcs in the spherical grid
% and the quantity of solutions in the approximated Pareto front that we
% are seeking (to improve clarity). 

spMODEDat.Alphas=30;
spMODEDat.PFrontSize=150;

% Execute the algorithm

OUT4=spMODEII(spMODEDat)

% Ploting the Tolerable, Desirable and Highly Desirable vectors for both 
% instances.
figure; 
plot3(0.40,0.6,0.9,'ob'); hold on; grid on;
plot3(0.10,0.4,0.8,'ob');
plot3(0.05,0.3,0.5,'ob');

% Ploting an usual run without preferences (base case).
plot3(OUT.PFront(:,1),OUT.PFront(:,2),OUT.PFront(:,3),'xb');
% Ploting an usual run taking into account preferences ans size control. 
% Solutions with higher GPP will be discarded (less preferable solutions)
% but with, at most, 2 tolerable values.
plot3(OUT4.PFront(:,1),OUT4.PFront(:,2),OUT4.PFront(:,3),'dr');

%% FIFTH EXAMPLE: Constraints

clear spMODEDat

% Default parameters are used
spMODEparam;

% The following means that it will be used a spherical grid in two
% dimensions, but all objectives and all constraints are used to calculate
% the preferability of a solution and that information will be used in the
% pruning mechanism and selection mechanism.

spMODEDat.NOBJ = 2;
spMODEDat.NRES = 1;

% Now we define our set of preferences (mandatory) in the format:
% J-ith objective : HD  -  D -  T  -  U  - HD
% HD Highly Desirable
% D  Desirable
% T  Tolerable
% U  Untolerable
% HU Highly Untolerable

spMODEDat.PhyMatrix{1} = [0.00  0.05  0.1  0.4  1.0   10;
                          0.00  0.30  0.4  0.6  1.0   10;
                          0.90  0.90  0.9  0.9  1.0   10]; 

% The following means that solutions with more than 3 objectives values in 
% the tolerable region are not allowed:

spMODEDat.PhyIndexMax=PhyIndex([0.4, 0.6, 0.9],spMODEDat);

% We modify the default values of the number of arcs in the spherical grid
% and the quantity of solutions in the approximated Pareto front that we
% are seeking (to improve clarity). 

spMODEDat.Alphas=30;
spMODEDat.PFrontSize=150;

% Execute the algorithm

OUT5=spMODEII(spMODEDat)

% Ploting the Tolerable, Desirable and Highly Desirable vectors for both 
% instances.
figure;
plot3(0.40,0.6,0.9,'ob'); hold on; grid on;
plot3(0.10,0.4,0.9,'ob');
plot3(0.05,0.3,0.9,'ob');

% Ploting an usual run without preferences (base case).
plot3(OUT.PFront(:,1),OUT.PFront(:,2),OUT.PFront(:,3),'xb');
% Ploting an usual run taking into account preferences.
plot3(OUT5.PFront(:,1),OUT5.PFront(:,2),OUT5.PFront(:,3),'dr');

%% Release and bug report:
%
% June 2014: Initial release