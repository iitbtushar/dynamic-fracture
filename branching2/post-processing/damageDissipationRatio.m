% computation of Damage Dissipation ratio [Bleyer, IJF 2017]

clc;clear all; close all;

iplot = 2; %2=> sigma3;||3=> sigma6

iprintData = 0; %1= to print out the text data

dtime = 1e6; xtol = 1e-3; x0 = 50; L = 100; %length of specimen

if iplot == 1; fileName1 = 'loading2/wu05-ft-energies.dat'; fileName2 = 'loading2/wu05-ft_tips.dat'; char1 = 'L2a'; end
if iplot == 2; fileName1 = 'loading2/wu3-ft-energies.dat'; fileName2 = 'loading2/wu3-ft_tips.dat'; char1 = 'L2b'; end
if iplot == 3; fileName1 = 'loading2/wu6-ft-energies.dat'; fileName2 = 'loading2/wu6-ft_tips.dat'; char1 = 'L2c'; end

R = xlsread('researchers/bobaru-results.xlsx',char1);


% material properties
E = 72e3; nu = 0.22; Gf = 3.8e-3; rho = 2440*10^(-12); ft = 28.65; %Bobaru(2015)

mu = E/(2*(1+nu)); Cs = sqrt(mu/rho); 
Cr = (0.862 + 1.14* nu)/(1+nu)* Cs; % Rayleigh wave speed %Cr = 3172e3 mm/sec;

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

% get rid of infinity at start/end
dEdx = zeros(length(vA.t),1); %dE/dx
dE1dx = zeros(length(vA.t),1); %dE1/dx
for i = 1:length(vA.t)
    if vA.s(i)> x0+xtol && vA.s(i)< L-xtol
        dE1dx(i) = vA.dE1ds(i); dEdx(i) = vA.dEds(i);
    end  
end

if iprintData == 1
% a set of data for table in paper
fout = fopen('bobaru-energy-rate-computation.csv','w');
fprintf(fout,'t,x,y,psiF,psiS,psiK,psiS+psiK,a, svel, dE/da, -dE1/da\n');
for i = 1:length(vA.t)
fprintf(fout,'%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n',vA.t(i), vA.s(i), 25, vA.psiF(i), vA.psiS(i), vA.psiK(i),...
    (vA.psiS(i) + vA.psiK(i)), vA.s(i), vA.svel(i), dEdx(i), -dE1dx(i));
end
fclose(fout);
end


%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',1.5);set(0, 'DefaultAxesFontSize',30);


%-----------------------------------------------------------------------
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

yyaxis left;
p1 = plot(vA.s, vA.svel/Cr,'-b','LineWidth',1.5);

% if iplot == 2; aa =[19.6 19.6]; bb = [0 0.6]; end;
% if iplot == 3; aa =[9.8 9.8]; bb = [0 0.6]; end;
% if (iplot == 2) ||(iplot == 3);p00 = plot(aa,bb,'--k','LineWidth',1.5); end;

ylabel('crack tip velocity $\hat{v}/C_R$','interpreter','latex','FontSize',30); 
set(gca,'XMinorTick','on','YMinorTick','on');
ylim([0 0.8]);

yyaxis right;
q1 = plot(vA.s, dEdx/Gf,'r-','LineWidth',1.5); 
q0 = plot([0 100],[2 2],'k--', 'LineWidth',2.0);

%if (iplot == 2) ||(iplot == 3);  patch('vertices', [0 2; 100 2; 100 3; 0 3],'faces', [1, 2, 3, 4],'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);end;


% %crack branching
if iplot == 2; patch('vertices', [64 0; 65 0; 65 5.0; 64 5.0],'faces', [1, 2, 3, 4],'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); end
if iplot == 3; patch('vertices', [57 0; 58 0; 58 5.0; 57 5.0],'faces', [1, 2, 3, 4],'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); end

ylim([0 7]);
ylabel('energy dissipation rate $\mathcal{G}/G_f$','interpreter','latex','FontSize',30); 

xlabel('crack length $l$ [mm]','interpreter','latex','FontSize',30); 
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');

% leg1=legend([ p1 p0 ],{'\texttt{PF-CZM}','Peridynamics: Bobaru and Zhang (2015)'},'NumColumns',1,'Location','northwest','interpreter','latex','FontSize', 24);
% set(leg1,'Box','off');

if (iplot == 2); xlim([0 75]);end;
if (iplot == 3); xlim([0 75]);end;

%ax=gca;ax.XAxis.TickLabelFormat='%,.1f'; ax.YAxis.TickLabelFormat=['%,.1f' '%,.1f'];


%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[22.5 20];fig.PaperPosition = [0. 0. 22.5 20];%fig.PaperPositionMode = 'auto';

if iplot == 2; fileName = 'bobaru-L2b-energy-x'; end
if iplot == 3; fileName = 'bobaru-L2c-energy-x'; end
print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');



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
