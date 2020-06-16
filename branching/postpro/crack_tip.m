clc;clear all;close all;

%iplot = 3;%1=std|2=marigo|3=pfczm
for iplot = 1:3

dtime      = 10^6; format = '%f %f %f'; Cr = 2119; % m/sec % thickness = 1.0; %mm

if iplot == 1 
    fileName   = 'sigma1b05/branching-std_tips.dat'; [time,x,y] = textread(fileName,format); [vA] = computeVelocity(time,x,y);
    fileName   = 'sigma1b025/branching-std_tips.dat'; [time,x,y] = textread(fileName,format); [vB] = computeVelocity(time,x,y);
end

if iplot == 2 
    fileName   = 'sigma1b05/branching-marigo_tips.dat'; [time,x,y] = textread(fileName,format); [vA] = computeVelocity(time,x,y);
    fileName   = 'sigma1b025/branching-marigo_tips.dat'; [time,x,y] = textread(fileName,format); [vB] = computeVelocity(time,x,y);
end

if iplot == 3 
    fileName   = 'sigma1b05/branching-implicit_tips.dat'; [time,x,y] = textread(fileName,format); [vA] = computeVelocity(time,x,y);
    fileName   = 'sigma1b025/branching-implicit_tips.dat'; [time,x,y] = textread(fileName,format); [vB] = computeVelocity(time,x,y);
end

%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',2.0);set(0, 'DefaultAxesFontSize',30);

%----------------------------------------------------------------------
%---- Compare load deformation -------------------
%----------------------------------------------------------------------
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot(dtime*vA.t,vA.svel/1000,'b-');
p2 = plot(dtime*vB.t,vB.svel/1000,'r-');

plot([0 80],[0.6*Cr 0.6*Cr],'k--','LineWidth',2.5);
text(75,0.63*Cr,'$0.6C_R$','HorizontalAlignment','center','interpreter','latex','FontSize',24);

xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30);
ylabel('velocity [m/s]','interpreter','latex','FontSize',30);

leg1=legend([ p1 p2 ],{'$b = 0.50$ mm','$b = 0.25$ mm'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 24);
 set(leg1,'Box','off');

%--
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';ax.YAxis.TickLabelFormat='%,.0f';

xlim([0 80]); ylim([0 1500]); yticks([0 300 600 900 1200 1500]);

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

if iplot == 1; fileName = 'branching-velocity-std'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); end;
if iplot == 2; fileName = 'branching-velocity-marigo'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); end;
if iplot == 3; fileName = 'branching-velocity-wu'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); end;
end

%---------------------------------- Functions ------------------------------------
function [velocity] = computeVelocity(time,x,y)

% not every step is used
interval = 1;

time = time(1:interval:length(time));
x    = x(1:interval:length(x));

% remove duplicated values 
%[C,IA,IC] = unique(x); time = time(IA); x    = C; y    = y(IA); %(unique x values)
[C,IA,IC] = unique(time); x = x(IA); time    = C; y    = y(IA); %(unique time values)

timeCount = length(time);
dLength   = zeros(timeCount,1); % incremental crack lengths
cLength   = zeros(timeCount,1); % total crack lengths 
dvelocity  = zeros(timeCount,1); % 
svelocity = zeros(timeCount,1); % smooth velocity

for i=2:timeCount
    dLength(i) = sqrt( (x(i)-x(i-1))^2 + (y(i)-y(i-1))^2 );
end

for i=2:timeCount
    cLength(i) = sum(dLength(1:i));
end

for i=4:length(time)-4
    %velocity(i) = (x(i)-x(i-1))/(time(i)-time(i-1));
    dvelocity(i) = (x(i+1)-x(i))/(time(i+1)-time(i));
    p = polyfit([time(i-3) time(i-2) time(i-1) time(i) time(i+1) time(i+2) time(i+3)],[x(i-3) x(i-2) x(i-1) x(i) x(i+1) x(i+2) x(i+3)],1);
    svelocity(i)=p(1);
end

velocity.t = time;
velocity.dvel = dvelocity;
velocity.svel = svelocity;

end