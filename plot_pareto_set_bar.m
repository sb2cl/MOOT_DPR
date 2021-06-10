f1 = figure('Color',[1 1 1],'Colormap',Blue_pink15);
str = ["p_a","C_{Na}","p_h",...
       "C_{Nh}","k_c","k_{d20}","\mu"];

subplot(1,7,1);
b1 = barh(OUT.PSet(:,1),'grouped','FaceColor','flat')
 for k = 1:15
     b1(k).CData = Blue_pink15;
 end
xlabel(str(1),'Interpreter','tex')
ylabel('Pareto solution','Interpreter','tex')
 
 subplot(1,7,2);

b1 = barh(OUT.PSet(:,2),'grouped','FaceColor','flat')
 for k = 1:15
     b1(k).CData = Blue_pink15;
 end
  xlabel(str(2),'Interpreter','tex')
 
 subplot(1,7,3);
b1 = barh(OUT.PSet(:,3),'grouped','FaceColor','flat')
 for k = 1:15
     b1(k).CData = Blue_pink15;
 end
 xlabel(str(3),'Interpreter','tex')
 
 N = 4;
 subplot(1,7,N);
b1 = barh(OUT.PSet(:,N),'grouped','FaceColor','flat')
 for k = 1:15
     b1(k).CData = Blue_pink15;
 end
 xlabel(str(N),'Interpreter','tex')
 
  N = 5;
 subplot(1,7,N);
b1 = barh(OUT.PSet(:,N),'grouped','FaceColor','flat')
 for k = 1:15
     b1(k).CData = Blue_pink15;
 end
 xlabel(str(N),'Interpreter','tex')
 
   N = 6;
 subplot(1,7,N);
b1 = barh(OUT.PSet(:,N),'grouped','FaceColor','flat')
 for k = 1:15
     b1(k).CData = Blue_pink15;
 end
 xlabel(str(N),'Interpreter','tex')

   N = 7;
 subplot(1,7,N);
b1 = barh(OUT.PSet(:,N),'grouped','FaceColor','flat')
 for k = 1:15
     b1(k).CData = Blue_pink15;
 end
 xlabel(str(N),'Interpreter','tex')
 
%applyhatch_pluscolor(f1, 'x', 0, [1 0 1 0], jet(4));
%hold on;
