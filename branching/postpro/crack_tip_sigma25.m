clc;clear all;close all;

dtime      = 10^6; format = '%f %f %f'; Cr = 2119; % m/sec % thickness = 1.0; %mm

fileName   = 'sigma1b025/branching_tips.dat'; [time,x,y] = textread(fileName,format); [vA] = computeVelocity(time,x,y);
fileName   = 'sigma25/branching_tips.dat'; [time,x,y] = textread(fileName,format); [vB] = computeVelocity(time,x,y);

%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',1.5);set(0, 'DefaultAxesFontSize',20);

%----------------------------------------------------------------------
%---- Compare load deformation -------------------
%----------------------------------------------------------------------
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot(dtime*vA.t,vA.svel/1000,'b-','LineWidth',1.5);
p2 = plot(dtime*vB.t,vB.svel/1000,'r-','LineWidth',1.5);

plot([0 80],[0.6*Cr 0.6*Cr],'k--','LineWidth',2.5);
text(75,0.62*Cr,'$0.6v_R$','HorizontalAlignment','center','interpreter','latex','FontSize',20);

xlabel('time [$\mu$s]','interpreter','latex','FontSize',20);
ylabel('velocity [m/s]','interpreter','latex','FontSize',20);

leg1=legend([ p1 p2 ],{'$\sigma = 1.0$ MPa','$\sigma = 2.5$ MPa'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 18);
 set(leg1,'Box','off');

%--
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';ax.YAxis.TickLabelFormat='%,.0f';

xlim([0 80]); ylim([0 1800]); %yticks([0 200 400 800 1000 1200 1400]);

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[15 15];fig.PaperPosition = [0. 0. 15 15];%fig.PaperPositionMode = 'auto';

fileName = 'branching-velocity-sigma2'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); 


%---------------------------------- Functions ------------------------------------
function [velocity] = computeVelocity(time,x,y)

% not every step is used
interval = 50;

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