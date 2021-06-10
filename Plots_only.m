%load OUT_spMODE_20210115T043652.mat;
load OUTPATH;
OUT15 = OUT;

%Blue_pink = ['395CAF','4366B1','4C6FB3','5579B6','5E82B8','678CBB','7095BD','7A9FC0','84A8C3','8FB1C6','9CBAC9','AAC2CD','BBCAD1','D0D0D7','EDD3DF'];
Blue_pink = ['172e12','44732c','5d8740','779a57','91ae70','acc28a','c7d6a6','e9c9bc','d7aea6','c29494','ac7c84','936577','78506c','5a3d64','362c5f'];
   
Blue_pink15 = hex2rgb(Blue_pink);
mycolormap = [Blue_pink15;0 0 0];

%pfront = ["Titer error(mg/L)","% error before/after pert.","ITAE of CHS",...
       %"# Oscillation","Sigma molecule"];

pfront = ["Titer error(mg/L)","% error before/after pert.",...
       "# Oscillation","Sigma molecule"];
%%
figure; 
mycolormap = colormap(Blue_pink15);
index = OUT15.Param.NOBJ+OUT15.Param.NRES;
for j=1:index
    subplot(2,3,j);
    for k=1:size(OUT.PFront,1)
        plot(k,OUT15.PFront(k,j),'o','MarkerEdgeColor',[73/255,0,106/255],'MarkerFaceColor',mycolormap(k,:)); 
        hold on;
    end
    if index==4
    plot(Jnominal4(j),'kv','MarkerFaceColor','k');
    else
    plot(Jnominal(j),'kv','MarkerFaceColor','k');
    end
    title(pfront(j));
end
colormap(gcf,mycolormap)
%
str = ["pa","CNa","ph",...
       "CNh","kc","kd20","mu growth"];
   
figure; 
for j=1:7
    subplot(2,4,j);
    for k=1:size(OUT.PSet,1)
        plot(k,OUT15.PSet(k,j),'o','MarkerEdgeColor',[73/255,0,106/255],'MarkerFaceColor',mycolormap(k,:));
        hold on;
    end
    plot(xnominal(j),'kv','MarkerFaceColor','k')
    title(str(j))
end

%% LEVEL DIAGRAM
% load OUT_spMODE_20210113T151106.mat;
% OUT13 = OUT;
% load OUT_spMODE_20210114T052202.mat;
% OUT14 = OUT;
% load OUT_spMODE_20210115T043652.mat;
% OUT15 = OUT;
% load OUT_spMODE_20210126T131404.mat;
% OUT26 = OUT;
% load OUT_spMODE_20210211T022357.mat;
% OUT11 = OUT;
% 
% % Load and create Concepts 1, 2 and 3
% %load data
% conceptCreate(OUT13.PFront(:,1:4),OUT13.PSet,'concept1')

% 
% % Calculate the 2-norm normalizing the front.
% % Bounds for Pareto front normalization
% max_all=max([concept1.maxpf,concept2.maxpf,concept3.maxpf,concept4.maxpf,concept5.maxpf]);
% min_all=min([concept1.minpf,concept2.minpf,concept3.minpf,concept4.minpf,concept5.minpf]);
% 


%% LEVEL DIAGRAM
conceptCreate(OUT15.PFront(:,1:3),OUT15.PSet,'concept3');
conceptCreate(Jnominal4(:,1:3),xnominal,'concept4');

% Calculate the 2-norm normalizing the front.
% Bounds for Pareto front normalization
max_all=max([concept3.maxpf,concept4.maxpf]);
min_all=min([concept3.minpf,concept4.minpf]);
bounds= [max_all; min_all];
%bounds=[concept3.maxpf;concept3.minpf];
basicNorm('ld1','concept3',bounds,2);
basicNorm('ld1','concept4',bounds,2);

ldDraw('ld1','concept3');
hold on;
ldDraw('ld1','concept4');
ldChangeMarker(ld1,'v','concept4')

%c=winter(concept3.nind); %colormap
% Creating an array with sizes for each point
s5=(20:5:20+5*(concept3.nind-1))';

%Ordering colors and sizes wrt parameter x1
%[nil0,idx0]=sort(concept3_data(:,1));
%c2(idx0,:)=c;
ldChangeColor(ld1,Blue_pink15,'concept3')
ldChangeColor(ld1,[0 0 0],'concept4')

%Ordering sizes wrt 3rd Objective: oscillations of sigma
[nil,idx]=sort(concept3_data(:,3));
%
%%%Creating s5 for sizes
[uv,~,idX] = unique(nil);
ne = accumarray(idX(:),1);
s5=[20.*ones(1,ne(1)), 65.*ones(1,ne(2)), 175.*ones(1,ne(3)) 260.*ones(1,ne(4))]';

s6(idx,1)=s5;
ldChangeSize(ld1,s6,'concept3')



%% Metabolites
%load('OUTPATH.mat');

metabolite=["L-ty","p-Co Acid", "p-CoA", "Malonyl","N chal","Naringenin","Dikaem","Kaemp"];
pathway=["Sigma","Antisigma", "Complex SA", "LuxR","AHLi","CHS","QdoR",...
          "L-ty","p-Co Acid", "p-CoA","N chal","Naringenin","Dikaem","Kaemp"];
metabolite2 = ["L-ty", "p-CoA", "Malonyl","Naringenin"];


%%
f1 = figure('Color',[1 1 1],'Colormap',Blue_pink15);

Titers= [1000-OUT.PFront(:,1), TITER(1:end-1,6)];
%Titers= [TITER(1:end-1,6)];
subplot(1,1,1);
b1 = bar(Titers,'grouped','FaceColor','flat')%Naringenin
for k = 1:size(Titers,2)
    b1(k).CData = Blue_pink15;
end

applyhatch_pluscolor(f1, 'x', 0, [1 0 1 0], jet(4));
hold on;

%nominal
plot([1:(size(TITER,1))],TITER(end,6).*ones(1,size(TITER,1)),'-v','LineWidth',2,...
    'Color',[0 0 0]);


%% scatter3d
sizeOscillation =zeros(size(OUT.PFront,1),1);

f2 = figure('Color',[1 1 1],'Colormap',Blue_pink15);
for j=1:size(OUT.PFront,1)
    if OUT.PFront(j,3)==0.5
        sizeOscillation(j)=20;
    elseif OUT.PFront(j,3)==1
        sizeOscillation(j)=75;
    elseif OUT.PFront(j,3)==1.5
        sizeOscillation(j)=120;
    elseif OUT.PFront(j,3)==2
        sizeOscillation(j)=185;
    else OUT.PFront(j,3)==2.5
        sizeOscillation(j)=270;
    end
scatter(OUT.PFront(j,1),OUT.PFront(j,2),sizeOscillation(j),'MarkerEdgeColor','k',...
        'MarkerFaceColor',Blue_pink15(j,:),'MarkerFaceAlpha',0.8)
hold on;
end
scatter(Jnominal4(:,1),Jnominal4(:,2),80,'v','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0.25 0.25 0.25]); 
grid on;
xlabel('Titer target error (mg/L)'); 
ylabel('Production before-after perturbation(%)');
%zlabel('Antithetic controller(# oscillations)');




























