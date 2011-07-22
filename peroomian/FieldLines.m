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


x_min=-120;x_max=0.;
z_min=-5.;z_max=5.

iFile_End = 88



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTENS = 'gsmb';
% RUNE = 'March0808';
% RUNEVENT = strcat(RUNE,'event');
% RUN = strcat(RUNE,EXTENS);
% TIME1 = '1240'
% TIME2 = '1250'
% TIME3 = '1300'

RUN = '001'

DIR_OUT=strcat('./plots/');

DIR_DATA=strcat('./flines/');

pngfile=strcat(DIR_OUT,'fieldlines.',RUN,'.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%initialize plot axes limits
left1 = .1;bottom1 = .10;width1 = .7;height1 = .60; 
delete(figure(1))
zero = [-20000.0 0.0 ; 20000.0 0.0];
set (gcf,...
    'PaperUnits','inches','PaperPosition',[0. 0. 11. 8.],...
    'PaperOrientation','portrait','PaperType','usletter');
FontName=['Helvetica'];
axes('Position',[left1,bottom1,width1,height1])
plot(zero(:,1),zero(:,2),'--k','linewidth',0.5);
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

iFile=1
while( iFile <= iFile_End)
    
    ppart=num2str(iFile,'%3.3i');
    file_fline=strcat(DIR_DATA,'fline.',ppart,'.out')
    
    
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %
%                  Read in field line...
%                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
icounter = 0;
int8 istart=0
fid = fopen(file_fline);
BB = fscanf(fid, '%lf %lf %lf ', [3, inf]);
BB=BB';
[nx,ny]=size(BB)
for i=1:nx  
    
    x_part(i)               = -BB(i,1);
    z_part(i)               = BB(i,3);
    
end
fclose(fid);

plot(x_part,z_part,'-k','LineWidth',1);

iFile = iFile + 1
clear x_part;
clear z_part;
end


axis([x_min x_max z_min  z_max]);
                      
%title(PARTICLE,'FontSize',16,'FontWeight','Bold','color','k')
xlabel('X (R_E)','FontSize',12);
ylabel('Z (R_E)','FontSize',12);

%set(gca,'xdir','reverse','ydir','reverse','xticklabel',[],'TickDir','out','XMinorTick','on', 'YMinorTick','on',...

set(gca,'xdir','reverse','TickDir','out','XMinorTick','on', 'YMinorTick','on',...
                          'FontSize',12,...
                      'xtick',[-120.0,-110.0,-100.0,-90.0,...
                      -80.0,-70.0,-60.0,-50.0,-40.0,...
                      -30.0,-20.0,-10.0,0.0])
                  
%                   'xtick',[-120.0,-115.0,-110.0,-105.0,-100.0,-95.0,-90.0,...
%                       -85.0,-80.0,-75.0,-70.0,-65.0,-60.0,-55.0,-50.0,-45.0,-40.0,-35.0,...
%                       -30.0,-25.0,-20.0,-15.0,-10.0,-5.0,0.0])
                  
% title=strcat('Particle # ',PARTICLE);
% posx=(x_min - x_max)/2. ;
% posy=z_max+2;
% text(posx,posy,title,'FontSize',20,'FontWeight','Bold','color','r','HorizontalAlignment','center')
%                            
                           
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


print ('-r144','-dpng',pngfile)
% saveas(gcf,pngfile, 'png')
% saveas(gcf,jpgfile, 'png')


break

