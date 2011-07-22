%set(gcf,'PaperUnits','centimeters'
%'Color',[1 1 1],'colormap',colormap(jet));
% Can be replaced by my personal map:
% Search for number of string matches per line.  
clear all

%addpath('C:\Users\mostafa\Matlab\matlab_Mostafa\function');

%addpath('/Applications/MATLAB74/functions');
Rearth = 6371.2;
boltzmann=1.3807e-23;
pmass= 1.6726e-27;
echarge=1.6022e-19;
y = 0;
shiftm=46.0;


x_min=-50;x_max=0.;
y_min=-5.;y_max=25.
z_min=-5.;z_max=5.


time_min=0.0;time_max=60.;

ener_min = 0;ener_max = 3.0




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTENS = 'gsmb';
% RUNE = 'March0808';
% RUNEVENT = strcat(RUNE,'event');
% RUN = strcat(RUNE,EXTENS);
% TIME1 = '1240'
% TIME2 = '1250'
% TIME3 = '1300'

PARTICLE = '00001'

DIR_OUT=strcat('./plots/');

DIR_DATA=strcat('./particles/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
file_part=strcat(DIR_DATA,'particle.',PARTICLE);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pngfile=strcat(DIR_OUT,'particle.',PARTICLE,'.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%initialize plot 1 axes limits
left1 = .1;bottom1 = .70;width1 = .7;height1 = .20;           
left2 = .1;bottom2 = .42;width2 = .7;height2 = .20;           
left3 = .1;bottom3 = .10;width3 = .7;height3 = .20;           


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %
%                  Read in particle...
%                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
icounter = 0;
int8 istart=0
fid = fopen(file_part);
BB = fscanf(fid, '%lf %lf %lf %lf %lf %lf %lf %lf %lf ', [9, inf]);
BB=BB';
[nx,ny]=size(BB);
for i=1:nx  
    t_part(i)               = BB(i,2)/60.;
    x_part(i)               = -BB(i,3);
    y_part(i)               = -BB(i,4);
    z_part(i)               = BB(i,5);
    en_part(i)               = BB(i,9);
end





% for i=1:istart
%     x_part_k(i) = x_part_k(istart)
%     y_part_k(i) = y_part_k(istart)
%     z_part_k(i) = z_part_k(istart)
% end
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set (gcf,...
    'PaperUnits','inches','PaperPosition',[0. 0. 11. 8.],...
    'PaperOrientation','portrait','PaperType','usletter');
FontName=['Helvetica'];


zero = [-20000.0 0.0 ; 20000.0 0.0];


delete(figure(1))
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Bx component     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(4,1,1);              

% AXES-1
axes('Position',[left1,bottom1,width1,height1])
plot(x_part,z_part,'-k','LineWidth',1);

hold on
plot(zero(:,1),zero(:,2),'--k','linewidth',0.5);

% set(gca,'xdir','reverse','xticklabel',[],'TickDir','out','XMinorTick','on', 'YMinorTick','on')
axis([x_min x_max z_min  z_max]);
                      
%title(PARTICLE,'FontSize',16,'FontWeight','Bold','color','k')
ylabel('Z (R_E)','FontSize',12);

%set(gca,'xdir','reverse','ydir','reverse','xticklabel',[],'TickDir','out','XMinorTick','on', 'YMinorTick','on',...

set(gca,'xdir','reverse','TickDir','out','XMinorTick','on', 'YMinorTick','on',...
                          'FontSize',12,...
                      'xtick',[-120.0,-115.0,-110.0,-105.0,-100.0,-95.0,-90.0,...
                      -85.0,-80.0,-75.0,-70.0,-65.0,-60.0,-55.0,-50.0,-45.0,-40.0,-35.0,...
                      -30.0,-25.0,-20.0,-15.0,-10.0,-5.0,0.0])
                  
title=strcat('Particle # ',PARTICLE);
posx=(x_min - x_max)/2. ;
posy=z_max+2;
text(posx,posy,title,'FontSize',20,'FontWeight','Bold','color','r','HorizontalAlignment','center')
                           
                           
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AXES-2
axes('Position',[left2,bottom2,width2,height2])

plot(x_part,y_part,'-k','LineWidth',1);
hold on;

plot(zero(:,1),zero(:,2),'--k','linewidth',0.5);
axis([x_min x_max y_min  y_max]);
%set(gca,'FontName','Helvetica','FontSize',16,'FontWeight','Bold')
xlabel('X (R_E)','FontSize',12);
ylabel('Y (R_E)','FontSize',12);
%ylabel('B_{y} (nT)','FontSize',12);

% set(gca,'xticklabel',[],'TickDir','out','XMinorTick','on', 'YMinorTick','on',...
%                           'FontSize',12,'xtick',[4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,...
%                           14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0])
set(gca,'xdir','reverse','ydir','reverse','TickDir','out','XMinorTick','on', 'YMinorTick','on',...
                          'FontSize',12,...
                      'xtick',[-120.0,-115.0,-110.0,-105.0,-100.0,-95.0,-90.0,...
                      -85.0,-80.0,-75.0,-70.0,-65.0,-60.0,-55.0,-50.0,-45.0,-40.0,-35.0,...
                      -30.0,-25.0,-20.0,-15.0,-10.0,-5.0,0.0])


                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
% AXES-3
axes('Position',[left3,bottom3,width3,height3])

plot(t_part,en_part,'-k','LineWidth',1);
hold on;
axis([time_min time_max ener_min  ener_max]);

xlabel('Time (Minutes)','FontSize',12);
ylabel('Energy (keV)','FontSize',12);
 set(gca,'TickDir','out','XMinorTick','on', 'YMinorTick','on',...
     'FontSize',12)
 %'xtick',[0.0,1.0,2.0,3.0,4.0])
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


print ('-r144','-dpng',pngfile)
% saveas(gcf,pngfile, 'png')
% saveas(gcf,jpgfile, 'png')


break

