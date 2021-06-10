% Naringerin Metabolic, Anthitetic controller and QdoR biosensor model.
% Updated 17/09/2019 Alejandro Vignoni, Yadira Boada

function [dxdt] = model(t,x,p)

Ncell = p.Ncell;
Size = p.Size; 

for k = 1:Ncell
    m = Size*(k-1)+1;      

%% Genetic model
%x1 = Sigma (act)    
c1 = p.ps(k)*p.CN(k)*p.ks(k)./( p.dms(k)+p.mu );
dxdt(m,1)=c1*( p.alpha(k) +(1-p.alpha(k))* x(m+4)^2./( p.kdlux(k)*(p.kd2(k)*p.CN(k)./x(m+3))^2 + x(m+4)^2) )-...
          p.k_c(k)./p.kdc(k)*x(m)*x(m+1) + p.k_c(k)*x(m+2) - (p.ds(k)+p.mu)*x(m);
      
      
%x2 = AntiSigma 
c2 =p.kdq(k)*p.CNa(k);
c3 =p.kdk(k)+ x(m+13); %Kae=x(m+13)
dxdt(m+1,1)= p.pa(k)*p.CNa(k)*p.ka(k)./(p.dma(k)+p.mu)*( p.alpha(k) +...
             (1-p.alpha(k))*(c2^2*c3^2)./( c2^2*c3^2 + (p.kdk(k)*x(m+6))^2))-...
             p.k_c(k)./p.kdc(k)*x(m)*x(m+1) +p.k_c(k)*x(m+2)-(p.da(k)+p.mu)*x(m+1);

%x3 = Complex Sigma-Antisigma
dxdt(m+2,1) = p.k_c(k)/p.kdc(k)*x(m)*x(m+1) - p.k_c(k)*x(m+2)-(p.dc(k)+p.mu)*x(m+2);

%x4 = LuxR
dxdt(m+3,1)= p.pR(k)*p.CN(k)*p.kR(k)./(p.dmR(k)+p.mu) - (p.dR(k)+p.mu)*x(m+3); 

%x5 = AHLint
dxdt(m+4,1) = p.D*p.Vcell/p.Vext*x(Size*Ncell+1)-p.D*x(m+4)- (p.dA(k)+ p.mu)*x(m+4);
           
%x6 = CHS
c6 = p.CNh(k)*p.kh(k)./( p.dmh(k)+p.mu );
dxdt(m+5,1) = c6*p.beta(k)*p.phc(k)+ c6*p.ph(k)*( p.alpha(k) + ...
              (1-p.alpha(k))*x(m)^2./( p.kd20*(p.kds(k)*p.CNh(k)) + x(m)^2 ) )-...
              (p.dh(k)+p.mu)*x(m+5);

%x7 = AHLext
%x8 = QdoR
dxdt(m+6,1) =  p.pq(k)*p.CN(k)*p.kq(k)/(p.dmq(k)+p.mu)-(p.dq(k)+p.mu)*x(m+6); 

%% Metabolic model
%x9 = L-tyrosine
%dxdt(m+7,1) =  -p.catp(k)*p.TAL*x(m+7)./(p.KmLt(k)+x(m+7))-p.mu*x(m+7); 
dxdt(m+7,1) =  p.KLt - p.catp(k)*p.TAL*x(m+7)./(p.KmLt(k)+x(m+7))-p.mu*x(m+7); 

%x10 = p-Coumeric acid
dxdt(m+8,1) = p.catp(k)*p.TAL*x(m+7)./(p.KmLt(k)+x(m+7)) -...
              p.catA(k)*p.CL4*x(m+8)./(p.KmP(k)+x(m+8))-p.mu*x(m+8); 

%x11 = p-CoA, CHS=x(m+5)
dxdt(m+9,1) = p.catA(k)*p.CL4*x(m+8)./(p.KmP(k)+x(m+8)) -...
              p.catNc(k)*x(m+5)*(x(m+9)*p.Mal3./( 1/(p.Ka(k)*p.Kb(k))+x(m+9)/p.Kb(k)+p.Mal3/p.Ka(k)+x(m+9)*p.Mal3))-...
              p.mu*x(m+9); 
                    
%x12 = Naringenin chalcone
dxdt(m+10,1) = p.catNc(k)*x(m+5)*(x(m+9)*p.Mal3./( 1/(p.Ka(k)*p.Kb(k))+x(m+9)/p.Kb(k)+p.Mal3/p.Ka(k)+x(m+9)*p.Mal3))-...
               p.catN(k)*p.CHI*x(m+10)./(p.KmNc(k)+x(m+10))- p.mu*x(m+10);
           
%x13 = Naringenin
dxdt(m+11,1) = p.catN(k)*p.CHI*x(m+10)./(p.KmNc(k)+x(m+10))-...
               p.catD(k)*p.F3H*x(m+11)./(p.KmN(k)+x(m+11)) - p.mu*x(m+11);

%x14 = Dihydrokaempferol
dxdt(m+12,1) = p.catD(k)*p.F3H*x(m+11)./(p.KmN(k)+x(m+11))-...
               p.catK(k)*p.FLS*x(m+12)./(p.KmD(k)+x(m+12)) - p.mu*x(m+12);
                      
%x15 = Kaempferol
dxdt(m+13,1) = p.catK(k)*p.FLS*x(m+12)./(p.KmD(k)+x(m+12)) - p.mu*x(m+13);

%x16 = Number of cells
dxdt(m+14,1) = p.mu*x(m+14)*(1-x(m+14)/p.cellmax);

end
%AHLext
    m=1;  %x(m+4)=AHL of Cell1
    dxdt((Size*Ncell+1),1) = -p.D*x(m+14)*p.Vcell/p.Vext*x(Size*Ncell+1) +...
                              p.D*x(m+14)*sum(x(m+4:Size:Size*Ncell)) - p.dAe*x(Size*Ncell+1);
                              %dxdt((Size*Ncell+1),1) = -p.D*x(m+14)*p.Vcell/p.Vext*x(Size*Ncell+1) +...
                              %p.D*sum(x(m+4:Size:Size*Ncell)) - p.dAe*x(Size*Ncell+1);
end