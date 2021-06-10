%   These are the parameters of the Anthitetic circuit Sigma-Anti-sigma and
%   LuxR protein.
%   Update 05/08/2019 by Yadira Boada

function [p] = parameters(Ncell,sd, ODmax)
%%%%%%%%%%%%%%%%%%%%%%%%  General parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    p.Ncell = Ncell;
    p.doubling = 82;          % doubling time [min]
    p.mu = log(2)/p.doubling;        % growth rate [1/min]
    p.nA = 6.023e23;                          % Avogadro's number: # particles/mol
    p.Vcell =  1.1e-15;                   % typical volume of E. coli (liters). Source: Bionumbers
    p.Vext = 4e-3;                   %culture medium volume [l] in a plate reader.  microfluidic device = 1e-9 liters
                                     % From Solvej Siedler, Novel biosensors based on flavonoid
    p.OD_to_cells = 8e11;               % cells/liter for OD=1, Agilent, E. coli Cell Culture
    p.cellmax = ODmax*p.Vext*p.OD_to_cells;     % maximum number of cells
    p.molecules_to_M = 1/(p.Vcell*p.nA)*1e6; % uM Conversion factor from number of particles to concentration (1 nMolar=nanomols/liter)

for k = 1:Ncell
    p.w(k) = 0;   % PERTURBATION
    %p.pN_luxI = 17;                   % plasmid number  pBR322 (15-20 copies/cell)
    %Copy number 
    p.CN(k) = 10 + sd*randn(1);   % plasmid number  pACYC184 (10 copies/cell)
    %Copy number antisigma
    p.CNa(k) = 10 + sd*randn(1);   % plasmid number  pACYC184 (10 copies/cell)
    %Copy number CHS
    p.CNh(k) = 10 + sd*randn(1);
    
    %Sigma
    p.dms(k) = log(2)/3 + sd*randn(1);             % degradation rate mRNA [1/min]
    p.ds(k) = 0.0003 + sd*randn(1);                 % degradation rate  [1/min]. Rapid degradation: RR Burgess, doi: 10.1006/rwgn.2001.1192
    p.ps(k) =  3 + sd*randn(1); %b*dms;    % translation rate  [1/min] [2.9801 - 5.9603] from our rates calculator
    p.ks(k) = 2*0.99 + sd*randn(1); %200/b;          % transcription rate [1/min] [0.99338 - 9.9338] from our rates calculator
    
    p.kd20(k) = 1000 + sd*randn(1);        % dissociation cte to promoter [molecules], Annunziata 2017
    
    %Anti-sigma
    p.alpha(k) = 0.01 + sd*randn(1);             % basal expresion pLux
    p.dma(k) = log(2)/3 + sd*randn(1);           % mRNA degradation rate [1/min]
    p.da(k) = 0.0003 + sd*randn(1);               % protein degradation rate [1/min]
    p.pa(k) =  0.8*3.96 + sd*randn(1);               % translation rate  [1/min] [3.9648 - 7.9295]
    p.ka(k) = 1.5*1.32 + sd*randn(1);              % transcription rate [1/min] [1.3216 - 13.2159] from our rates calculator
    
    %cI
    p.kdcI(k) = 30;
    p.kd_lamcI(k) = 1000;                                % dissociation cte to promoter [molecules]. Strong promoter
    p.phc_cI(k) = 6.5096e-04;
    p.ph_cI(k) = 0.001*1.54;
    
    %SigmaComplex
    p.kdc = 0.01 + sd*randn(1);                  % dissociation constant [molecules] Annunziata 2017 an orthogonal multi-input
    p.k_c = 1.8e-3 + sd*randn(1);                % [1/min] Annunziata 2017 an orthogonal multi-input
    %p.kc = p.k_c/p.kdc + sd*randn(1);            % binding rate sigma to anti-sigma [min^-1 molecules^-1]
    p.dc = 0.001 + sd*randn(1);                   % degradation rate [1/min] Annunziata 2017 an orthogonal multi-input
    
    %LuxR
   
    p.beta1(k) = 0.01 + sd*randn(1);                % basal expresion p20
    p.dmR(k) = log(2)/3 + sd*randn(1);              % mRNA degradation rate  [1/min]
    p.dR(k) = 0.02 + sd*randn(1);                   % degradation rate [1/min]
    p.pR(k) = 2.34 + sd*randn(1);                  % translation rate  [1/min] [1.1749 - 2.3499] from our rates calculator
    p.kR(k) = 2*0.39 + sd*randn(1);                   % transcription rate [1/min] [0.39164 - 3.9164] from our rates calculator
    
    % Monomer LuxR.AHL
    p.kd1(k) = 100 + sd*randn(1);                   % dissociation constant of R to A [nM], Urbanowski etal. 2004
    p.k_1(k) = 10 + sd*randn(1);                   % unbinding rate LuxR to AHL [1/min]
    p.k1(k) = p.k_1(k)/p.kd1(k) + sd*randn(1);            % binding rate LuxR to AHL [1/min]
    p.dRA(k) = log(2)/5 + sd*randn(1);              % degradation rate of (LuxR.A) [1/min]. Buchler et al. 2004 Monomer half-life is just few minutes.
    
    % Dimer (R.A)2
    p.kd2(k) = 20 + sd*randn(1);                   %dissociation cte (LuxR.A) to (LuxR.A) [nM], Buchler et al. 2003
    %Koren, R. & Hammes, G. G. (1976) Biochemistry 15, 1165�1171.
    %Northrup, S. H. & Erickson, H. P. (1992) Proc. Natl. Acad. Sci. USA 89, 3338�3342.
    p.k_2(k) = 1 + sd*randn(1);                     % dissociation rate [1/min]
    p.k2(k) = p.k_2(k)/p.kd2(k) + sd*randn(1);               % binding rate LuxR to AHL [1/min]
    p.kdlux(k) = 600 + sd*randn(1);           % dissociation cte (LuxR.A)2 to promoter [nM], Bucler et al [1 1000]nM
    
    %Sigma dimer
    p.kds(k) = 1000 + sd*randn(1);
    
    %CHS (389 amino acids)
    p.beta(k) = 1.5 + sd*randn(1);  
    p.dmh(k) = log(2)/3 + sd*randn(1);              % degradation rate mRNA [1/min]
    p.dh(k) = 0.0003 + sd*randn(1);                 % degradation rate  [1/min]. 
    
    %p.ph(k) =   1.557*2*1.54 + sd*randn(1);        %Open loop gain. Only comparative plots
    p.phc(k) =  2*1.54 + sd*randn(1); % Re-written in COST.m   CONSTITUTIVE translation rate   
    p.ph(k) = 60*1.54 + sd*randn(1);                %translation rate [1/min] [1.5424 - 3.0848] from our rates calculator  %RBS of the constitutive promoter
    p.kh(k) = 0.02* 48*7.5*0.51 + sd*randn(1);   %%%CHANGING 1/10        %   transcription rate [1/min] [0.51414 - 5.1414] from our rates calculator
        
    %QdoR (842 bp)
    p.dmq(k) = log(2)/3 + sd*randn(1);             % degradation rate mRNA [1/min]
    p.dq(k) = 0.0003 + sd*randn(1);                % degradation rate  [1/min]. 
    p.pq(k) =  1.2*2.13 + sd*randn(1); %b*dms;         % translation rate  [1/min] [2.1378 - 4.2755] from our rates calculator
    p.kq(k) = 0.71 + sd*randn(1); %200/b;          % transcription rate [1/min] [0.71259 - 7.1259] from our rates calculator
    p.kdq = 150 + sd*randn(1);                     % kd = 1 nM dissociation constant complex (Q.Kae)2 from promoter PqdoI.  'Novel biosensors based on 
                                                   %flavonoid-responsive transcriptional regulators introduced into E. coli'
    
    %Dissociation constant Kaempferol from QdoR. kdq>kdk
    p.kdk = 75 + sd*randn(1);              %kd^n = 5 uM.  Ref.: Novel biosensors based on 
                                                   %flavonoid-responsive transcriptional regulators introduced into E.
    
    %AHL and AHLe
    p.D = 2;        %kinetic rate of AHL external transport [1/min] across the cell membrane, calculated
    
    p.dA(k) = 0.0004 + sd*randn(1);    %[0.05 0.03 min^-1]Degradation from Bionumbers online
    p.dAe(k) = 0.0000481 + sd*randn(1);          % Horswill et al., 2007  %0.0164, Degradation rate for AHL. From Kaufmann etal. 2005. Similar to You etal. Nature, 2004
    %[0.05 0.03 min^-1]Degradation from Bionumbers online
    %0.000282 Degradation rate for external AHL. Fitted using half-life of 180 minutes, from Englmann etal. 2007
    % In Kauffmann & Sartorio, 2005 they use 0.0018
    
    %TAL,L-tyrosine
    p.tyrosine0= 2e6;       % [nM]From 2,89e-5M Wild type e.coli table. Refs:3mM in the culture
    p.KLt = p.tyrosine0; %[nM]2e6;%/9.05;   %This is 1.75 mM L-tyrosine initial in culture (other Refs 3mM)
    p.KmLt(k) = 1.9e-5*1e9+ sd*randn(1);  %[1.9E-5M, 1.6E-4]M 
    p.catp(k) = 0.02*60 + sd*randn(1);    % 0.02- 4.32 (1/sec)
    
    %4CL, p-Coumeric acid
    p.catA(k) = 8.20e-3*60+ sd*randn(1);           %[8.20e-3, 8.87e1] 1/sec
    p.KmP(k) =  1.4e-5*1e9+ sd*randn(1);            %[1.40E-5, 0.000432]M
       
    %CHS,p-CoA   %%%% haciendo esta 10 veces o 20 mas rapida
    p.catNc(k) = 4*0.007*60+ sd*randn(1);    %[0.007,0.042]	 s^(-1)
    p.Ka(k) =  1/(1.0e-6*1e9) + sd*randn(1);    %[1.00E-6, 3.50E-2]M
    p.Kb(k) =  1/(1.0e-6*1e9) + sd*randn(1);    %malonyl-CoA  https://www.brenda-enzymes.org/enzyme.php?ecno=2.3.1.74
    
    %CHI, Naringenin chalcone
    p.catN(k)= 0.07*60+ sd*randn(1);    %[0.07,23]	 s^(-1)
    p.KmNc(k) = 2.8e-5*1e9+ sd*randn(1); %[2.40E-06, 3.54E-04]M   https://www.brenda-enzymes.org/enzyme.php?ecno=5.5.1.6
    
   %F3H, Naringenin  %%% Esta cambiada
    p.catD(k)= 5.8*0.5*1*60 + sd*randn(1);     %[1.00E+00	3.90E+00]s^(-1)
    p.KmN(k) = 100000*5e-6*1e9 + sd*randn(1);   %[5.00E-06, 2.18E-04]M https://www.brenda-enzymes.org/enzyme.php?ecno=1.14.11.11
                                                % this value is 2293 times less effective
   %FLS, Dihydrokaempferol
    p.catK(k)= 0.1*1*60 + sd*randn(1);       %6.6 (1/sec) https://www.brenda-enzymes.org/enzyme.php?ecno=1.14.20.6
    p.KmD(k) = 10e-6*1e9 + sd*randn(1);    %[1.00E-06,5.84E-05]M    https://www.brenda-enzymes.org/enzyme.php?ecno=1.14.20.6

end