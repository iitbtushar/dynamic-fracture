% computation of Damage Dissipation ratio [Bleyer, IJF 2017]

clc;clear all; close all;

iplot = 2; %1=> v/cR || 2=> G/Gf

dtime = 1e6; xtol = 1e-3; x0 = 50; L = 100; %length of specimen

fileName1 = 'sigma1b05-KE/wu-energies.dat'; fileName2 = 'sigma1b05-KE/wu_tips.dat';  
fileName3 = 'sigma25-random/mesh1/wu-energies.dat'; fileName4 = 'sigma25-random/mesh1/wu_tips.dat';


% material properties
E = 3.2e4; nu = 0.2; ft = 12; Gf = 3e-3; rho = 2450*10^(-12);%dynamic-branching

mu = E/(2*(1+nu)); Cs = sqrt(mu/rho); 
Cr = (0.862 + 1.14* nu)/(1+nu)* Cs; % Rayleigh wave speed %Cr = 3172e3 mm/sec;

%=============== Sigma0 = 1MPa =======================================
% energy-input
fileName   = fileName1; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f';
A          = textscan(fileID, formatSpec, 'HeaderLines', 2); A          = cell2mat(A); % convert cells to matrices
% time % fracture energy  % strain energy % kinetic energy
time1 = A(:,1); ef = A(:,2); es = A(:,3); ek = A(:,4); 


% tip-location-input
fileName   = fileName2; format = '%f %f %f';
[time,x,y] = textread(fileName,format);

% get smooth velocity and energy
[vA] = computeVelocityEnergy(time,x,y, time1,ef,es,ek);

% get rid of start/end
xmin = x0+xtol; xmax = L-xtol; 
for i = 1:length(vA.t); if vA.s(i)> xmin; i1 =i; break; end; end;
for i = 1:length(vA.t); if vA.s(i)> xmax; i2 =i; break; end; end;
vA.s = vA.s(i1:i2); vA.t = vA.t(i1:i2);
vA.dEds = vA.dEds(i1:i2); vA.svel = vA.svel(i1:i2);

%pre-append zero
vA.s = [0.0; x0; vA.s]; vA.t = [0.0; vA.t(1); vA.t]; vA.dEds = [0.0; 0.0; vA.dEds]; vA.svel = [0.0; 0.0; vA.svel]; 

%=============== Sigma0 = 2.5MPa =======================================
% energy-input
fileName   = fileName3; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f';
A          = textscan(fileID, formatSpec, 'HeaderLines', 2); A          = cell2mat(A); % convert cells to matrices
% time % fracture energy  % strain energy % kinetic energy
time1 = A(:,1); ef = A(:,2); es = A(:,3); ek = A(:,4); 

% tip-location-input
fileName   = fileName4; format = '%f %f %f';
[time,x,y] = textread(fileName,format);

% get smooth velocity and energy
[vB] = computeVelocityEnergy(time,x,y, time1,ef,es,ek);

% get rid of start/end
xmin = x0+xtol; xmax = 93; 
for i = 1:length(vB.t); if vB.s(i)> xmin; i1 =i; break; end; end;
for i = 1:length(vB.t); if vB.s(i)> xmax; i2 =i; break; end; end;
vB.s = vB.s(i1:i2); vB.t = vB.t(i1:i2);
vB.dEds = vB.dEds(i1:i2); vB.svel = vB.svel(i1:i2);

%pre-append zero
vB.s = [0.0; x0; vB.s]; vB.t = [0.0; vB.t(1); vB.t]; vB.dEds = [0.0; 0.0; vB.dEds]; vB.svel = [0.0; 0.0; vB.svel]; 

%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',2.0);set(0, 'DefaultAxesFontSize',30);

%-- v/cR ----
%-----------------------------------------------------------------------
if iplot == 1
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot(vA.s, vA.svel/Cr,'-b');
p2 = plot(vB.s, vB.svel/Cr,'-r');

% %crack branching
% sigma0 = 1.0
patch('vertices', [63.5 0; 65, 0; 65 0.8; 63.5 0.8],'faces', [1, 2, 3, 4],'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
%sigma0 = 2.5
patch('vertices', [54 0; 56 0; 56 0.8; 54 0.8],'faces', [1, 2, 3, 4],'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 
patch('vertices', [74 0; 76 0; 76 0.8; 74 0.8],'faces', [1, 2, 3, 4],'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);


xlabel('crack length $l$ [mm]','interpreter','latex','FontSize',30); 
ylabel('crack tip velocity $\hat{v}/C_R$','interpreter','latex','FontSize',30); 
ylim([0 0.8]); xlim([0 L]);
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f'; ax.YAxis.TickLabelFormat='%,.1f';

leg1=legend([ p1 p2 ],{'$\sigma^\ast = 1.0$ MPa','$\sigma^\ast = 2.5$ MPa'},'NumColumns',1,'Location','northwest','interpreter','latex','FontSize', 22);
set(leg1,'Box','off');

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

fileName = 'branching1-vCr'; 
%print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');
end


%-- G/Gf ----
%-----------------------------------------------------------------------
if iplot == 2
figure(2); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p1 = plot(vA.s, vA.dEds/Gf,'-b');
p2 = plot(vB.s, vB.dEds/Gf,'-r');
q0 = plot([0 L],[2 2],'k--', 'LineWidth',2.0);
q0 = plot([0 L],[4 4],'k--', 'LineWidth',2.0);

% %crack branching
% sigma0 = 1.0
patch('vertices', [63.5 0; 65, 0; 65 7; 63.5 7],'faces', [1, 2, 3, 4],'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
%sigma0 = 2.5
patch('vertices', [54 0; 56 0; 56 7; 54 7],'faces', [1, 2, 3, 4],'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 
patch('vertices', [74 0; 76 0; 76 7; 74 7],'faces', [1, 2, 3, 4],'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);


ylim([0 8.0]); xlim([0 L]);
ylabel('energy dissipation rate $\mathcal{G}/G_f$','interpreter','latex','FontSize',30); 
xlabel('crack length $l$ [mm]','interpreter','latex','FontSize',30); 
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');

leg1=legend([p1 p2 ],{'$\sigma^\ast = 1.0$ MPa','$\sigma^\ast = 2.5$ MPa'},'NumColumns',1,'Location','northwest','interpreter','latex','FontSize', 22);
set(leg1,'Box','off');

ax=gca;ax.XAxis.TickLabelFormat='%,.0f'; ax.YAxis.TickLabelFormat='%,.1f';

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

fileName = 'branching1-GGf'; 
%print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');
end


%---------------------------------- Functions ------------------------------------
function [velocity] = computeVelocityEnergy(time,x,y,time1,ef,es,ek)

tol_velocity = 1.0; %tolerance in velocity to compute dE/dx

% not every step is used
interval = 30; %30
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

velocity.psiF = eF; %crack surface energy
velocity.psiS = eS; %strain energy
velocity.psiK = eK; %kinetic energy
velocity.dEds = dEdx; %dE/dx (based on kinetic + strain energy)
velocity.dE1ds = dE1dx; %dE/dx (based on crack surface energy)

end
