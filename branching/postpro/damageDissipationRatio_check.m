% computation of Damage Dissipation ratio [Bleyer, IJF 2017]

clc;clear all; close all;

iplot = 1; %1=> v/Cr || 2=> G/Gf || 3 => energy

formatSpec1 = '%f %f %f %f'; formatSpec2 = '%f %f %f'; 

dtime = 1e6; xtol = 1e-1; x0 = 50.0; L = 100; %length of specimen 

fileName1 = 'sigma25-random/mesh1/wu-energies.dat'; fileName2 = 'sigma25-random/mesh1/wu_tips.dat'; 
fileName3 = 'sigma25-random-reduceMass/wu-energies.dat'; fileName4 = 'sigma25-random-reduceMass/wu_tips.dat'; 

%======================================================================
% energy-input
fileName   = fileName1; fileID     = fopen(fileName,'r'); A = textscan(fileID, formatSpec1, 'HeaderLines', 20); A = cell2mat(A); time1 = A(:,1); ef = A(:,2); es = A(:,3); ek = A(:,4); 

% tip-location-input
fileName   = fileName2; [time,x,y] = textread(fileName,formatSpec2);

% get smooth velocity and energy
[vA] = computeVelocityEnergy(time,x,y, time1,ef,es,ek); 
L = 100; xmin = x0+xtol; xmax = L-xtol; [vvA] = removeEdges(vA,x0,xmin,xmax);

%======================================================================
fileName   = fileName3; fileID     = fopen(fileName,'r'); A = textscan(fileID, formatSpec1, 'HeaderLines', 20); A = cell2mat(A); time1 = A(:,1); ef = A(:,2); es = A(:,3); ek = A(:,4); 
fileName   = fileName4; [time,x,y] = textread(fileName,formatSpec2);
[vB] = computeVelocityEnergy(time,x,y, time1,ef,es,ek); 
L = 100; xmin = x0+xtol; xmax = L-xtol; [vvB] = removeEdges(vB,x0,xmin,xmax);



% material properties
E = 3.2e4; nu = 0.2; ft = 12; Gf = 3e-3; rho = 2450*10^(-12);%dynamic-branching

mu = E/(2*(1+nu)); Cs = sqrt(mu/rho); 
Cr = (0.862 + 1.14* nu)/(1+nu)* Cs; % Rayleigh wave speed %Cr = 3172e3 mm/sec;


%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',2.0);set(0, 'DefaultAxesFontSize',30);

%---v/Cr--------
%======================================================================
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot(vvA.s, vvA.svel/Cr,'-k');
p2 = plot(vvB.s, vvB.svel/Cr,'-r');

% p20 = plot(vvB.s(1:10:end), vvB.svel(1:10:end)/Cr,'o','MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','b');
% p30 = plot(vvC.s(1:10:end), vvC.svel(1:10:end)/Cr,'^','MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','b');

%ylim([-0.05 1.0]); 
ylabel('crack tip velocity $\hat{v}/C_R$','interpreter','latex','FontSize',30); 

% % %crack branching
% patch('vertices', [64 0; 65 0; 65 5.0; 64 5.0],'faces', [1, 2, 3, 4],'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
% patch('vertices', [57 0; 58 0; 58 5.0; 57 5.0],'faces', [1, 2, 3, 4],'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);

set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
%xlim([-0.5 L]);
xlabel('crack length $l$ [mm]','interpreter','latex','FontSize',30); 

% leg1=legend([p1 p2 p3 p4 p5 p6],{'$U^\ast = 0.035$ mm','$U^\ast = 0.040$ mm','$U^\ast = 0.060$ mm',...
%     '$U^\ast = 0.08$ mm','$U^\ast = 0.10$ mm','$U^\ast = 0.14$ mm'},...
%     'NumColumns',3,'Location','northwest','interpreter','latex','FontSize', 18);
% set(leg1,'Box','off');

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

%fileName = 'zhou-vCr'; 
%print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');



%---G/Gf------------
%======================================================================
figure(2); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot(vvA.s, vvA.dEds/Gf,'-k');
p2 = plot(vvB.s, vvB.dEds/Gf,'-r');

%xlim([-0.5 L]); ylim([-1.0 25]); 
xlabel('crack length $l$ [mm]','interpreter','latex','FontSize',30); 
ylabel('energy dissipation rate $\mathcal{G}/G_f$','interpreter','latex','FontSize',30); 

%----------------------------------
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');

% leg1=legend([p1 p2 p3 p4 p5 p6],{'$U^\ast = 0.035$ mm','$U^\ast = 0.040$ mm','$U^\ast = 0.060$ mm',...
%     '$U^\ast = 0.08$ mm','$U^\ast = 0.10$ mm','$U^\ast = 0.14$ mm'},...
%     'NumColumns',3,'Location','northwest','interpreter','latex','FontSize', 18);
% set(leg1,'Box','off');

ylim([0 9]);

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

fileName = 'zhou-GGf'; 
%print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');


%---Energy------------
%======================================================================
figure(3); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot(vvA.t*dtime, vvA.es,'-k');
p2 = plot(vvB.t*dtime, vvB.es,'-r');

q1 = plot(vvA.t*dtime, vvA.ef,'--k');
q2 = plot(vvB.t*dtime, vvB.ef,'--r');

%xlim([-0.5 32]); %ylim([-1.0 25]); 
xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30); 
ylabel('$\Psi_s$, $\Psi_c$ [mJ]','interpreter','latex','FontSize',30); 

%----------------------------------
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');

% leg1=legend([p1 p2 p3 p4 p5 p6],{'$U^\ast = 0.035$ mm','$U^\ast = 0.040$ mm','$U^\ast = 0.060$ mm',...
%     '$U^\ast = 0.080$ mm','$U^\ast = 0.100$ mm','$U^\ast = 0.140$ mm'},...
%     'NumColumns',3,'Location','northeast','interpreter','latex','FontSize', 18);
% set(leg1,'Box','off');

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

%fileName = 'zhou-strain-energy'; 
%print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');

%======================================================================
figure(4); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot(vvA.t*dtime, vvA.ek,'-k');
p2 = plot(vvB.t*dtime, vvB.ek,'-r');

%xlim([-0.5 32]); %ylim([-1.0 25]); 
xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30); 
ylabel('$\Psi_k$ [mJ]','interpreter','latex','FontSize',30); 

%----------------------------------
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');

% leg1=legend([p1 p2 p3 p4 p5 p6],{'$U^\ast = 0.035$ mm','$U^\ast = 0.040$ mm','$U^\ast = 0.060$ mm',...
%     '$U^\ast = 0.080$ mm','$U^\ast = 0.100$ mm','$U^\ast = 0.140$ mm'},...
%     'NumColumns',3,'Location','northwest','interpreter','latex','FontSize', 18);
% set(leg1,'Box','off');

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

%fileName = 'zhou-kinetic-energy'; 
%print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');



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
end


%========================================================================
function [velocity] = computeVelocityEnergy(time,x,y,time1,ef,es,ek)

tol_velocity = 1.0; %tolerance in velocity to compute dE/dx

% not every step is used
interval = 15; %30
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