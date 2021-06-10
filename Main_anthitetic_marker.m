%   Naringerin Metabolic pathway, Anthitetic controller and 
%   QdoR biosensor model.
%   Parameters: structure contains all rates and constants
%   Updated 24/02/2021 by Yadira Boada, Alejandro Vignoni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set to 1 if this script is being used for the perturbation loop
PLOT_ON = 0; %Put 1 to get plots.
PLOT_DETAILED = 0;
AHL_EXTERNAL = 5000 ;   %nM

Variance = 0;
Stdeviation =0;
Ncell = 1;
ODinitial = 0.001;
ODmax = 120;
p = parameters(Ncell,Stdeviation,ODmax);
Cellinitial = ODinitial*p.Vext*p.OD_to_cells;  

%AHLe concentration. Vector optional
ahle_nM = [3, 2500];    %Induction in the lab [nM]

nM = 1e-9;  %nM in Molarity
to_molecules = p.Vext*p.nA*nM;
ahle0 = ahle_nM.*to_molecules; %vector

%System size
NumberStates = 16; 

%Closed loop gain
p.phc = 1000*6.5096e-04; %= 0.8454*0.0005*1.54; %b*dms;    % translation rate  [1/min] [1.5424 - 3.0848] from our rates calculator  %RBS of the constitutive promoter
p.ph = 15.6230; % = 0.8454*12*1.54; %RBS of the inducible promoter
     
% Open loop
%p.phc =  6.5493;% = 0.96655* 4.4*1.54;      %Open loop gain. Only comparative plots
%p.ph = 0;

%Enzymes (molecules). Max values from each enzyme range.
p.TAL = 20*1.6e5;  
p.CL4 = 15*4.32e5;  
p.CHI = 10*3.54e5;
p.F3H = 2.81;
p.FLS = 5.84;
MAL_PERCENT = 0.4; %40% of Malonyl

X = OUT15.PSet;
Titer_before = zeros(length(X),8); %g/L
Titer_after  = zeros(length(X),8); 
Productivity_before = zeros(length(X),8); %g/(Lh)
Productivity_after = zeros(length(X),8);

for xpop=1:size(X,1) %size(X,1) calcula la cantidad total de candidatos.

    % Decision variables & parameters   
    p.pa = X(xpop,1); % RBS anti-sigma
    p.CNa = X(xpop,2); % Copy number anti-sigma
    p.ph = X(xpop,3); % RBS CHS enzyme
    p.CNh = X(xpop,4); % Copy number CHS enzyme
    p.kc = X(xpop,5); % binding rate sigma.Asigma complex 1/(molecule.min)
     p.k_c = p.kdc*p.kc;
    p.kd20 = X(xpop,6); % dissociation constant sigma-promoter (molecules)
    p.mu = X(xpop,7); % Growth rate 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0 Null initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
tfin = 60*16;     %simulation time
step = 0.1;
tspan = 0:step:tfin-step;
options = odeset('AbsTol',1e-8,'RelTol',1e-6);      % for ode function 

p.Mal3 = 0;             %Input: 3 Malonyl-CoA
  p.Mal30 = p.Mal3;     %Input: 3 Malonyl-CoA
Initial = [zeros(1,NumberStates-2) Cellinitial 0];  %ini conditions[species, cells, ahle]
p.Size = length(Initial)-1;
[t0,x0] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

%% 1 Adding Malonyl  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial = x0(end,1:end); %Initial conditions
tfin = 60*8;             %Tiempo de simulacion (min)
tspan = 0:step:tfin-step;

p.Mal3 = 1.17e3;     %Mean amount Malonyl-CoA=3.54e-5 (M) from table
                     %Maximum Malonyl-CoA=3.09e-3  - Minimum 4.05e-7 (M)from table
  p.Mal31 = p.Mal3;     %
[t1,x1] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

%% 2 Adding ahle 1st. Closing the loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Mal32 = p.Mal3;     %Input: 3 Malonyl-CoA

%AHLe=AHLe(1)
Initial = [x1(end,1:end-1) ahle0(1)]; %Initial conditions
tfin = 60*10; % Tiempo de simulacion (min)
tspan = 0:step:tfin-step;
[t2,x2] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

%% 3 Adding ahle 2nd  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Mal33 = p.Mal3;     %Input: 3 Malonyl-CoA

%AHLe=AHLe(2)
Initial = [x2(end,1:end-1) ahle0(2)]; %Initial conditions

tfin = 60*47; % Tiempo de simulacion (min)
tspan = 0:step:tfin-step;
[t3,x3] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

% %% 3 Continue closed loop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.Mal33 = p.Mal3;     %Input: 3 Malonyl-CoA
% Initial = x2(end,1:end); %Initial conditions
% 
% tfin = 60*42; % Tiempo de simulacion (min)
% tspan = 0:step:tfin-step;
% [t3,x3] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

%% 4 Malonyl Perturbation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Mal3 = MAL_PERCENT*p.Mal3;     %Input: 3 Malonyl-CoA
  p.Mal34 = p.Mal3;     %Input: 3 Malonyl-CoA
Initial = x3(end,1:end); %Initial conditions
tfin = 60*36; %31
tspan = 0:step:tfin-step;
[t4,x4] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

%% Dynamics together
sigma = [x1(:,1); x2(:,1); x3(:,1); x4(:,1)]; 
asigma = [x1(:,2); x2(:,2); x3(:,2); x4(:,2)];  
sa_complex = [x1(:,3); x2(:,3); x3(:,3); x4(:,3)]; 
luxR = [x1(:,4); x2(:,4); x3(:,4); x4(:,4)]; 
ahl = [x1(:,5); x2(:,5); x3(:,5); x4(:,5)]; 
chs = [ x1(:,6); x2(:,6); x3(:,6); x4(:,6)];  
qdoR = [x1(:,7); x2(:,7); x3(:,7); x4(:,7)]; 
Ltyrosine = [x1(:,8); x2(:,8); x3(:,8); x4(:,8)]; 
pC_acid = [x1(:,9); x2(:,9); x3(:,9); x4(:,9)]; 
p_CoA = [x1(:,10); x2(:,10); x3(:,10); x4(:,10)]; 
nchalcone = [ x1(:,11); x2(:,11); x3(:,11); x4(:,11)]; 
naringenin = [x1(:,12); x2(:,12); x3(:,12); x4(:,12)]; 
dykaempferol = [x1(:,13); x2(:,13); x3(:,13); x4(:,13)]; 
kaempferol = [x1(:,14); x2(:,14); x3(:,14); x4(:,14)]; 
ahle = [x1(:,16); x2(:,16); x3(:,16); x4(:,16)];
od = [x1(:,15); x2(:,15); x3(:,15); x4(:,15)]./(p.OD_to_cells*p.Vext);
malonyl = [p.Mal31.*ones(length(t1),1); p.Mal32.*ones(length(t2),1);...
           p.Mal33.*ones(length(t3),1); p.Mal34.*ones(length(t4),1)];

time = [t1; t2+t1(end); t3+t2(end)+t1(end); t4+t3(end)+t2(end)+t1(end)]./60;


%% Metabolites production mg/L  

%molecular_weight=[L-ty,p-Co acid,p-CoA, Malonyl,N chal,Naringenin,Dikaem,Kaemp ];   %g/mol
molecular_weight=[181.19, 164.047, 913.67, 853.6, 272.25, 272.25, 288.25, 286.23]';   %g/mol
titer = (([x4(end,8:10), malonyl(end), x4(end,11:14)]'.*molecular_weight.*p.OD_to_cells.*ODmax)./p.nA).*1e3;
productivity = titer./(time(end)/(60*24));% mg /(L*days)

%Saving data
Production_after(xpop,:) = production;
Titer_after(xpop,:) = titer;
%total_product = [production, total_product];

if p.Mal31==p.Mal34
    mal_percent = 0;
else
    mal_percent = (1-p.Mal34/p.Mal31)*100;
end
total_mal = [mal_percent, total_mal];

end



%%
% if (PLOT_ON)
figure1= figure('Color',[1 1 1]); 
subplot(211)
c = categorical({'L-ty','p-Co acid','p-CoA','Malonyl','N chal','Naringenin','Dikaem','Kaemp'});
barfig = bar(c,Titer_after(1,:)); ylabel('Production mg/L');
set(gca,'YGrid','on','YMinorTick','on','YScale','log');
ylim([0.0001 1500]);
subplot(212)
barfig = bar(c,production_y4(1:end-1,:)); ylabel('Production mg/L');
set(gca,'YGrid','on','YMinorTick','on','YScale','log');
ylim([0.0001 1500]);
% 
% T_points=round([length(t1), length(t1)+length(t2),length(t1)+length(t2)+length(t3)]);
% figure20 = figure('Color',[1 1 1]); 
% subplot(221)
% plot(time,naringenin,'-d','MarkerIndices',T_points,'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',1); 
% hold on;
% plot(time,NAringenin,'-d','MarkerIndices',T_points,'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',1); 
% ylabel({'Naringenin';'(molecules)'});
% 
% subplot(222)
% plot(time,chs,'-d','MarkerIndices',T_points,'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',1);
% hold on;
% plot(time,CHS,'-d','MarkerIndices',T_points,'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',1);
% ylabel('CHS','Interpreter', 'Latex', 'FontSize',12);  xlabel('time (min)');
% 
% subplot(223)
% plot(kaempferol,chs,'Linewidth',1);
% hold on;
% plot(KAempferol,CHS,'Linewidth',1);
% xlabel('Kaempferol'); ylabel('CHS'); 
% subplot(224)
% tranf = 1e3*272.25*p.OD_to_cells.*ODmax./p.nA;
% plot(naringenin(24010:end)*tranf,chs(24010:end),'Linewidth',1);
% hold on;
% plot(NAringenin(24010:end)*tranf,CHS(24010:end),'Linewidth',1);
% xlabel('Naringenin'); ylabel('CHS'); 
% end