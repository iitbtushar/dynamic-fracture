% computation of Damage Dissipation ratio [Bleyer, IJF 2017]

clc;clear all; close all;

iplot = 1; %1=> v/Cr || 2=> G/Gf || 3 => energy

thick = 5;%thickness of plate

formatSpec1 = '%f %f %f %f'; formatSpec2 = '%f %f %f'; 

dtime = 1e6; xtol = 0.5; vtol = 10.0;%mm/s

fileName1 = 'dats/d25-energies.dat'; fileName2 = 'dats/d25_tips.dat'; 
fileName3 = 'dats/d30-energies.dat'; fileName4 = 'dats/d30_tips.dat';

%======================================================================
% energy-input
fileName   = fileName1; fileID     = fopen(fileName,'r'); A = textscan(fileID, formatSpec1, 'HeaderLines', 20); A = cell2mat(A); time1 = A(:,1); ef = A(:,2); es = A(:,3); ek = A(:,4); 

% tip-location-input
fileName   = fileName2; [time,x,y] = textread(fileName,formatSpec2);

% get smooth velocity and energy
interval=30; [vA] = computeVelocityEnergy(time,x,y, time1,ef,es,ek, interval); 
x0 = 8; L = 70; 
xmin = x0+xtol; xmax = L-xtol; [vvA] = removeEdges(vA,x0,xmin,xmax);

%======================================================================
fileName   = fileName3; fileID     = fopen(fileName,'r'); A = textscan(fileID, formatSpec1, 'HeaderLines', 20); A = cell2mat(A); time1 = A(:,1); ef = A(:,2); es = A(:,3); ek = A(:,4); 
fileName   = fileName4; [time,x,y] = textread(fileName,formatSpec2);
interval = 30; [vB] = computeVelocityEnergy(time,x,y, time1,ef,es,ek,interval); 
x0 = 8; L = 70;
xmin = x0+xtol; xmax = L-xtol; [vvB] = removeEdges(vB,x0,xmin,xmax);



%materials
E = 4.39e3; nu = 0.37; Gf = 0.0933; rho = 1170*10^(-12); ft = 40; %Epoxy

mu = E/(2*(1+nu)); Cs = sqrt(mu/rho); 
Cr = (0.862 + 1.14* nu)/(1+nu) * Cs; % Rayleigh wave speed %Cr = 920e3 mm/sec;


%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',2.0);set(0, 'DefaultAxesFontSize',30);

%---v/Cr--------
%======================================================================
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot((vvA.t-vvA.tp)*dtime, vvA.svel/Cr,'-b');
p2 = plot((vvB.t-vvB.tp)*dtime, vvB.svel/Cr,'-r');


ylim([0.0 1.1]); 
ylabel('crack tip velocity $\hat{v}/C_R$','interpreter','latex','FontSize',30); 

% % %crack branching
% patch('vertices', [64 0; 65 0; 65 5.0; 64 5.0],'faces', [1, 2, 3, 4],'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
% patch('vertices', [57 0; 58 0; 58 5.0; 57 5.0],'faces', [1, 2, 3, 4],'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);

set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';ax.YAxis.TickLabelFormat='%,.1f';

xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30); 

leg1=legend([p1 p2 ],{'$a_0 = 25$ mm','$a_0 = 30$ mm'},...
    'NumColumns',1,'Location','northwest','interpreter','latex','FontSize', 30);
set(leg1,'Box','off');

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

fileName = 'jap-vCr'; 
%print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');



%---G/Gf------------
%======================================================================
figure(2); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot((vvA.t-vvA.tp)*dtime, vvA.dEds/Gf,'-b');
p2 = plot((vvB.t-vvB.tp)*dtime, vvB.dEds/Gf,'-r');

%xlim([-5 100]); ylim([0.0 3.5]); 
xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30); 
ylabel('energy dissipation rate $\mathcal{G}/G_f$','interpreter','latex','FontSize',30); 

%----------------------------------
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';ax.YAxis.TickLabelFormat='%,.1f';


leg1=legend([p1 p2 ],{'$a_0 = 25$ mm','$a_0 = 30$ mm'},...
    'NumColumns',1,'Location','northwest','interpreter','latex','FontSize', 30);
set(leg1,'Box','off');

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

fileName = 'jap-GGf'; 
%print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');


% %---Energy------------
% %======================================================================
% figure(3); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));
% 
% p1 = plot((vvA.t-vvA.tp)*dtime, vvA.es*thick,'-b');
% p2 = plot((vvB.t-vvB.tp)*dtime, vvB.es*thick,'-r');
% 
% 
% q1 = plot((vvA.t-vvA.tp)*dtime, vvA.ef*thick,'--b');
% q2 = plot((vvB.t-vvB.tp)*dtime, vvB.ef*thick,'--r');
% 
% %xlim([-5 100]);
% xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30); 
% ylabel('$\Psi_s$, $\Psi_c$ [mJ]','interpreter','latex','FontSize',30); 
% 
% %----------------------------------
% set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
% 
% leg1=legend([p1 p2 ],{'$d = 25$ mm','$d = 30$ mm'},...
%     'NumColumns',1,'Location','northeast','interpreter','latex','FontSize', 18);
% set(leg1,'Box','off');
% 
% %---
% fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';
% 
% fileName = 'jap-strain-energy'; 
% %print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');
% 
% %======================================================================
% figure(4); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));
% 
% p1 = plot((vvA.t-vvA.tp)*dtime, vvA.ek*thick,'-b');
% p2 = plot((vvB.t-vvB.tp)*dtime, vvB.ek*thick,'-r');
% 
% %xlim([-5 100]); 
% xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30); 
% ylabel('$\mathcal{K}$ [mJ]','interpreter','latex','FontSize',30); 
% 
% %----------------------------------
% set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
% 
% leg1=legend([p1 p2 ],{'$d = 25$ mm','$d = 30$ mm'},...
%     'NumColumns',1,'Location','northwest','interpreter','latex','FontSize', 18);
% set(leg1,'Box','off');
% 
% %---
% fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';
% 
% fileName = 'jap-kinetic-energy'; 
% %print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');



%---------------------------------- Functions ------------------------------------
%========================================================================
function [updatedVel] = removeEdges(vA,x0,xmin,xmax)
% get rid of start/end
for i = 1:length(vA.t); if vA.s(i)> xmin; i1 =i; break; end; end;
for i = 1:length(vA.t); if vA.s(i)> xmax; i2 =i; break; end; end;

vA.s = vA.s(i1:i2); vA.t = vA.t(i1:i2); 

vA.dEds = vA.dEds(i1:i2); vA.svel = vA.svel(i1:i2);

vA.es = vA.es(i1:i2); vA.ef = vA.ef(i1:i2); vA.ek = vA.ek(i1:i2);

%pre-append zero
vA.s = [0.0; x0; vA.s]; vA.t = [0.0; vA.t(1); vA.t]; vA.dEds = [0.0; 0.0; vA.dEds]; vA.svel = [0.0; 0.0; vA.svel]; 

vA.es = [vA.es(1); vA.es(2); vA.es]; vA.ef = [0.0; vA.ef(1); vA.ef]; vA.ek = [0.0; vA.ek(1); vA.ek];

% update
updatedVel.s = vA.s; updatedVel.t = vA.t; updatedVel.dEds = vA.dEds; updatedVel.svel = vA.svel; 
updatedVel.es = vA.es; updatedVel.ef = vA.ef; updatedVel.ek = vA.ek; 
updatedVel.tp = vA.t(1); %time of propagation
end


%========================================================================
function [velocity] = computeVelocityEnergy(time,x,y,time1,ef,es,ek, interval)

tol_velocity = 1.0; %tolerance in velocity to compute dE/dx

% not every step is used
%interval = 30; %30
time = time(1:interval:length(time)); x    = x(1:interval:length(x));

% remove duplicated values 
[C,IA,IC] = unique(time); time    = C; x = x(IA); y    = y(IA); %(unique time values)
[C,IA,IC] = unique(time1); time1    = C; ef = ef(IA); es = es(IA); ek = ek(IA); %(unique time values)

timeCount = length(time);
dLength   = zeros(timeCount,1); % incremental crack lengths
cLength   = zeros(timeCount,1); % total crack lengths 
svelocity = zeros(timeCount,1); % smooth velocity

eF = zeros(timeCount,1); % fracture energy
eS = zeros(timeCount,1); % strain energy
eK = zeros(timeCount,1); % kinetic energy
eT = zeros(timeCount,1); % strain energy + kinetic energy
dEdx = zeros(timeCount,1); %dE/dx
dE1dx = zeros(timeCount,1); %dE1/dx

% compute crack length
for i=2:timeCount
    dLength(i) = sqrt( (x(i)-x(i-1))^2 + (y(i)-y(i-1))^2 );
    cLength(i) = sum(dLength(1:i));
    x(i) = cLength(i);%arc length
end

% smooth velocity
for i=10:length(time)-10
    dvelocity(i) = (x(i+1)-x(i))/(time(i+1)-time(i));
    p = polyfit([time(i-7) time(i-6) time(i-5) time(i-4) time(i-3) time(i-2) time(i-1) time(i) time(i+1) time(i+2) time(i+3) time(i+4) time(i+5) time(i+6) time(i+7)],...
        [x(i-7) x(i-6) x(i-5) x(i-4) x(i-3) x(i-2) x(i-1) x(i) x(i+1) x(i+2) x(i+3) x(i+4) x(i+5) x(i+6) x(i+7)],1);
    svelocity(i)=p(1);
end

% interpolated energy
eF = interp1(time1,ef,time); eS = interp1(time1,es,time); eK = interp1(time1,ek,time);
eT = eK + eS;

% compute dE/dx
for i=4:length(time)-4
    if svelocity(i)> tol_velocity
    q = polyfit([x(i-3) x(i-2) x(i-1) x(i) x(i+1) x(i+2) x(i+3)],[eF(i-3) eF(i-2) eF(i-1) eF(i) eF(i+1) eF(i+2) eF(i+3)],1);
    dEdx(i)=q(1);
    q = polyfit([x(i-3) x(i-2) x(i-1) x(i) x(i+1) x(i+2) x(i+3)],[eT(i-3) eT(i-2) eT(i-1) eT(i) eT(i+1) eT(i+2) eT(i+3)],1);
    dE1dx(i)=q(1);  
    end
end


velocity.t = time; %time
velocity.s = x; %space
%velocity.dvel = dvelocity; %discrete velocity
velocity.svel = svelocity; %smooth velocity

velocity.ef = eF; %crack surface energy
velocity.es = eS; %strain energy
velocity.ek = eK; %kinetic energy
velocity.dEds = dEdx; %dE/dx (based on kinetic + strain energy)
velocity.dE1ds = dE1dx; %dE/dx (based on crack surface energy)

end