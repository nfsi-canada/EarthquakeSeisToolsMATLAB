% 2019-11-13
% This piece produces a 3D figure showing the true hypocenters, estimated 
% hypocenters, and stations

figure(30288)
clf
hold on

plot3(ex,ey,ez,'ko','MarkerFaceColor',[0.5,0.5,0.5])
plot3(exF,eyF,ezF,'ko','MarkerFaceColor','b')
plot3([ex';exF'],[ey';eyF'],[ez';ezF'],'r--')

plot3(sx,sy,zeros(ns,1),'kv','MarkerFaceColor','r','MarkerSize',7)

set(gca,'ZDir','reverse')
grid on
axis equal


xlabel('E (km)')
ylabel('N (km)')
zlabel('Depth (km)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Repeat zoomed in on events

figure(30289)
clf
hold on

plot3(ex,ey,ez,'ko','MarkerFaceColor',[0.5,0.5,0.5])
plot3(exF,eyF,ezF,'ko','MarkerFaceColor','b')
plot3([ex';exF'],[ey';eyF'],[ez';ezF'],'r--')

set(gca,'ZDir','reverse')
grid on
axis equal


xlim([min([ex;exF])-0.5, max([ex;exF])+0.5])
ylim([min([ey;eyF])-0.5, max([ey;eyF])+0.5])
zlim([min([ez;ezF])-0.5, max([ez;ezF])+0.5])

xlabel('E (km)')
ylabel('N (km)')
zlabel('Depth (km)')
