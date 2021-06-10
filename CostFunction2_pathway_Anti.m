function  [J,X,TITER,PRODUCTIVITY,YIELD,MOLECULEx2,MOLECULEx3,MOLECULEx4] = CostFunction2_pathway_Anti(X,Dat)
%[J,X,TITER,PRODUCTIVITY,YIELD,MOLECULEx2,MOLECULEx3,MOLECULEx4]
%% Load experimental data from spMODEparam (execute first) 

%% MOO
if strcmp(Dat.Plot_Cost_Func,'yes')
    Blue_pink = ['172e12','44732c','5d8740','779a57','91ae70','acc28a','c7d6a6','e9c9bc','d7aea6','c29494','ac7c84','936577','78506c','5a3d64','362c5f'];
   
   
    Blue_pink15 = hex2rgb(Blue_pink);
    mycolormap = [Blue_pink15;0 0 0];
    f1 = figure('Color',[1 1 1],'Colormap',mycolormap); 
    f2 = figure('Color',[1 1 1],'Colormap',mycolormap); 
    f3 = figure('Color',[1 1 1],'Colormap',mycolormap); 

    
    %Blue_pink = ['395CAF','4366B1','4C6FB3','5579B6','5E82B8','678CBB','7095BD','7A9FC0','84A8C3','8FB1C6','9CBAC9','AAC2CD','BBCAD1','D0D0D7','EDD3DF'];  
%     Blue_pink = ['5A72B0', '587BAF', '5784AF', '588CB0', '5B94B1',
%     '619BB2', '69A2B3', '73A9B5', '7FAFB7', '8CB5B9', '9CBABC', 'ACBEBF', 'BEC2C3', 'd2C6C7', 'E6C8CB'];
    
    %Blue_pink = ['5A72B0', '5C7AAe', '6083AC', '648BAA', '6992A7',
    %'709AA4', '78A1A0', '81A89C', '8CAF97', '98B592', 'A7BB8C', 'B7C085', 'C9C47C', 'DFC773', 'F8C867'];
    %Blue_pink = ['7d5d98', '7d6a99', '7d759a', '7d809b', '7f8b9c', '82959e', '879ea0', '8da8a2', '96b1a5', 'a0b9a8', 'adc1ac', 'bdc8b0', 'd0ceb4', 'e6d2ba', 'ffd6c0'];    
      
    %Save data
    %molec_weight=[L-ty,p-CoAcid, p-CoA, Malonyl,N chal,Naringenin,Dikaem,Kaemp ];   %g/mol
    molecularWEIGHT=[181.19, 164.047, 913.67, 853.6, 272.25, 272.25, 288.25, 286.23];   %g/mol
    TITER= zeros(size(X,1),8); %mg/L 8metabolites
    PRODUCTIVITY= zeros(size(X,1),8); %g/L/days)
    MOLECULEx2=zeros(size(X,1),14); %x2(end only) NO AHLe nor OD
    MOLECULEx3=zeros(size(X,1),14); %x3(end only) NO AHLe nor OD
    MOLECULEx4=zeros(size(X,1),14); % x4(end only) NO AHLe nor OD
    YIELD =zeros(size(X,1),2); %yield of Naringenin
end

% General parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

p.Mal3 = 0;  %Input: 3 Malonyl-CoA
  p.Mal30 = p.Mal3;  

Initial = [zeros(1,NumberStates-2) Cellinitial 0];  %ini conditions[species, cells, ahle]

p.Size = length(Initial)-1;
[t0,x0] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

%% 1 Adding Malonyl  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfin = 60*8;             %Tiempo de simulacion (min)
tspan = 0:step:tfin-step;

Initial = [x0(end,1:end-2) Cellinitial 0]; %Initial conditions

p.Mal3 = 1.17e3;     %(nM) Mean amount Malonyl-CoA=3.54e-5 (M) from table
                     %Maximum Malonyl-CoA=3.09e-3  - Minimum 4.05e-7 (M)from table
  p.Mal31 = p.Mal3;     %
[t1,x1] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

%% 2 Adding ahle 1st. Closing the loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Mal32 = p.Mal3;     %Input: 3 Malonyl-CoA

tfin = 60*10; % Tiempo de simulacion (min)
tspan = 0:step:tfin-step;

%AHLe=AHLe(1)
Initial = [x1(end,1:end-2) Cellinitial ahle0(1)]; %Initial conditions
[t2,x2] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

%% 3 Adding ahle 2nd  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Mal33 = p.Mal3;     %Input: 3 Malonyl-CoA

tfin = 60*47; % Tiempo de simulacion (min)
tspan = 0:step:tfin-step;

%AHLe=AHLe(2)
Initial = [x2(end,1:end-1) ahle0(2)]; %Initial conditions
[t3,x3] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);

%% 4 Malonyl Perturbation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Mal3 = MAL_PERCENT*p.Mal3;     %Input: 3 Malonyl-CoA
  p.Mal34 = p.Mal3;     %Input: 3 Malonyl-CoA

tfin = 60*36; %31
tspan = 0:step:tfin-step;

Initial = x3(end,1:end); %Initial conditions
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
population = [x1(:,15); x2(:,15); x3(:,15); x4(:,15)];
od = [x1(:,15); x2(:,15); x3(:,15); x4(:,15)]./(p.OD_to_cells*p.Vext);
malonyl = [p.Mal31.*ones(length(t1),1); p.Mal32.*ones(length(t2),1);...
           p.Mal33.*ones(length(t3),1); p.Mal34.*ones(length(t4),1)];

time = [t1; t2+t1(end); t3+t2(end)+t1(end); t4+t3(end)+t2(end)+t1(end)]./60;%hours


%% Optimization
   
%Naringenin avegage in the last 1/4 of time t3 (before Ma perturbation)
element = 3*length(t3)/4;
Naringenin_x3 = mean(x3((element:end), 12));
Naringenin_x4 = mean(x4((element:end), 12));

% Naringenin's Titer (mg/L) without perturbation 
 molecular_weight=272.257; % molecular weight Naringenin (g/mol)
 J1 = abs(1000-(Naringenin_x3*molecular_weight*p.OD_to_cells*od(end))/p.nA*1e3); %(mg/L)

% Relative Error (%) before and after perturbation
 J2 = 100*(abs(Naringenin_x3 - Naringenin_x4)/Naringenin_x3);

% Sigma`s oscillations  before Ma perturbation
SIGMA = x3(:,1)';
zero_mean = SIGMA-mean(SIGMA);
zero_clipped = zeros(size(SIGMA));
zero_cross = 0;
for j = 2:length(zero_mean)
    if zero_mean(j)>=0
        zero_clipped(j) = 1;
    end
    zero_cross=zero_cross+(zero_clipped(j)-zero_clipped(j-1))^2;
end
J3 = zero_cross/2;  %oscillations

%Itae of chs. settling time of CHS after perturbation (min)
%  itae = cumsum(abs(x4(:,6)-x4(end,6)).*t4,1);
%  n=itae(end)*0.98;
%  [val,idx]=min(abs(itae-n));
%  J3 = t4(idx);

%% Constraints

%1  Sigma's molecules after adding AHL
sigma_t3t4 = [x3(:,1); x4(:,1)]';
asigma_t3t4 = [x3(:,2); x4(:,2)]';

if min(sigma_t3t4)>=3000
vector_size = length(sigma_t3t4);
vector = zeros(1,vector_size);
%1:asigma>sigma (NO), 0.5:asigma=sigma (NO), 0:asigma<sigma (YES)
    for j=1:vector_size
        vector(1,j) = (sign(asigma_t3t4(j)-sigma_t3t4(j))+1)/2; 
    end
    Value = sum(vector);
else
    Value = 1; 
end
    
J(xpop,:)=[J1, J2, J3, Value];

%% Plots

if strcmp(Dat.Plot_Cost_Func,'yes')
    figure(f1); 
        
    subplot(321,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,sigma.*p.molecules_to_M,'LineWidth',2);
    ylabel({'\sigma_{20}';'(uM)'});
    hold on;
    subplot(323,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,asigma.*p.molecules_to_M,'LineWidth',2);
    ylabel({'anti-\sigma_{20}';'(uM)'});
    hold on;
    subplot(325,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,sa_complex.*p.molecules_to_M,'LineWidth',2);
    ylabel({'(\sigma \cdot a\sigma)complex';'(uM)'}); xlabel('Time (h)')
    hold on;
    subplot(326,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap);  
    %naringeninTempo = (naringenin*molecular_weight./p.nA).*(population./p.Vext);
    %plot(time,naringeninTempo,'LineWidth',2);
    %ylabel({'Naringenin';'g L^{-1}'}); xlabel('Time (h)')
    plot(time,naringenin.*p.molecules_to_M,'LineWidth',2);
    ylabel({'Naringenin';'(uM)'}); xlabel('Time (h)')
    hold on;
    subplot(324,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,chs.*p.molecules_to_M,'LineWidth',2);
    ylabel({'CHS enzyme';'(uM)'});
    hold on;
    subplot(322,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,malonyl.*p.molecules_to_M,'LineWidth',2);
    ylabel({'Malonyl-CoA';'(uM)'});
    hold on;
    %%
    figure(f2);
    %if xpop == size(X,1)
        subplot(421,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
        plot(time,od,'LineWidth',2);
        ylabel({'OD_{600}'});xlabel('Time (h)'); hold on;
    %end  
    subplot(423,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap);  
    plot(time,Ltyrosine.*p.molecules_to_M,'LineWidth',2);
    ylabel({'Ltyrosine';'(uM)'}); hold on; 
    subplot(425,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,pC_acid.*p.molecules_to_M,'LineWidth',2);
    ylabel({'pCoumaroyl acid';'(uM)'}); hold on;
    subplot(427,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap);  
    plot(time,p_CoA.*p.molecules_to_M,'LineWidth',2);
    ylabel({'pCoumaroyl-CoA';'(uM)'});hold on;
    
    subplot(422,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,nchalcone.*p.molecules_to_M,'LineWidth',2);
    ylabel({'N chalcone';'(uM)'}); hold on;
    subplot(424,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,naringenin.*p.molecules_to_M,'LineWidth',2);
    ylabel({'Naringenin';'(uM)'}); hold on;
    subplot(426,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap);  
    plot(time,dykaempferol.*p.molecules_to_M,'LineWidth',2);
    ylabel({'Dykaempferol';'(uM)'}); hold on;
    subplot(428,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,kaempferol.*p.molecules_to_M,'LineWidth',2);
    ylabel({'Kaempferol';'(uM)'});  hold on; 
    
    %%  
    figure(f3);

    subplot(421,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,ahl.*p.molecules_to_M,'LineWidth',2);  ylabel({'AHL';'(uM)'}); 
    hold on;
    subplot(423,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,luxR.*p.molecules_to_M,'LineWidth',2); ylabel({'LuxR';'(uM)'}); hold on;
    subplot(425,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap);  
    plot(time,sigma.*p.molecules_to_M,'LineWidth',2); ylabel({'\sigma_{20}';'(uM)'}); hold on;
    subplot(427,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,chs.*p.molecules_to_M,'LineWidth',2); ylabel({'CHS enzyme';'(uM)'}); hold on;
    
    subplot(422,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,kaempferol.*p.molecules_to_M,'LineWidth',2); ylabel({'Kaempferol';'(uM)'}); hold on;
    subplot(424,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,qdoR.*p.molecules_to_M,'LineWidth',2);ylabel({'QdoR';'(uM)'}); hold on;
    subplot(426,'xlim',[0 100],'XGrid','on','YGrid','on','Colormap',mycolormap,'ColorOrder',mycolormap); 
    plot(time,asigma.*p.molecules_to_M,'LineWidth',2); ylabel({'anti-\sigma_{20}';'(uM)'}); hold on;

    %%
    
    %Integral of Ltyrosine
    time_feedback = [t3; t4+t3(end)]; %min
    cells_feedback = [x3(:,15); x4(:,15)]; %cells
    Lty_pop = p.KLt*trapz(time_feedback,cells_feedback);%molecules
    gramsNa = (naringenin(end)*population(end)/p.nA)*molecularWEIGHT(6);
    gramsLty = (Lty_pop/p.nA)*molecularWEIGHT(1);
    yieldNa =  gramsNa/gramsLty;%g g^-1
    
    %Metabolites production mg/L    
    titer = (([x4(end,8:10), malonyl(end,1), x4(end,11:14)].*molecularWEIGHT.*p.OD_to_cells.*od(end))./p.nA).*1e3; 
    %productivity = (titer./(time(end)/(60*24))).*0.001;% g /(L*days)
    productivity = (titer./(time_feedback(end)/60));% mg /(L*hours)
       
    %Saving data
    PRODUCTIVITY(xpop,:) = productivity;
    TITER(xpop,:) = titer;
    MOLECULEx2(xpop,:) = x2(end,1:end-2);
    MOLECULEx3(xpop,:) = x3(end,1:end-2);
    MOLECULEx4(xpop,:) = x4(end,1:end-2);
    YIELD(xpop,:) = [yieldNa, Lty_pop/(p.nA*p.Vext)];
end
end

end