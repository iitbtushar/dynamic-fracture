%
clc;clear all;close all;

dtime      = 10^6; formatSpec = '%f %f %f %f'; thickness = 1.0; %mm

fileName   = 'sigma25-random/mesh1/wu-energies.dat'; fileID     = fopen(fileName,'r'); A  = textscan(fileID, formatSpec, 'HeaderLines', 1); 
fileName   = 'sigma25-random/mesh2/wu-energies.dat'; fileID     = fopen(fileName,'r'); B  = textscan(fileID, formatSpec, 'HeaderLines', 1);

fileName   = 'sigma25/mesh2/branching-energies.dat'; fileID     = fopen(fileName,'r'); C  = textscan(fileID, formatSpec, 'HeaderLines', 1);

% convert cells to matrices
A          = cell2mat(A); B          = cell2mat(B); C          = cell2mat(C);

% % append zeros for first row
% lb = size(A); lb2 = lb(2);
% A = [zeros(1,lb2);A];

%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',1.5);set(0, 'DefaultAxesFontSize',30);

%----------------------------------------------------------------------
%---- Compare energies -------------------
%----------------------------------------------------------------------
fig = figure(1); left_color = [0 0 1]; right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

yyaxis left;
p1 = plot(dtime*A(:,1),A(:,3),'b-','LineWidth',1.5);
p2 = plot(dtime*A(:,1),A(:,2),'r-','LineWidth',1.5);

q1 = plot(dtime*B(:,1),B(:,3),'b--','LineWidth',1.5);
q2 = plot(dtime*B(:,1),B(:,2),'r--','LineWidth',1.5);

% r1 = plot(dtime*C(:,1),C(:,3),'b-.','LineWidth',1.5);
% r2 = plot(dtime*C(:,1),C(:,2),'r-.','LineWidth',1.5);

ylabel('Strain/Surface Energy [mJ]','interpreter','latex','FontSize',30);
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');

yyaxis right;
p3 = plot(dtime*A(:,1),A(:,4)/1000,'k-','LineWidth',1.5);

q3 = plot(dtime*B(:,1),B(:,4)/1000,'k--','LineWidth',1.5);

% r3 = plot(dtime*C(:,1),C(:,4)/1000,'k-.','LineWidth',1.5);

xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30);
ylabel('Kinetic Energy [J]','interpreter','latex','FontSize',30);

leg1=legend([ p1 q1 p2 q2 p3 q3],{'$\Psi_s$ (Mesh1)','$\Psi_s$ (Mesh2)','$\Psi_c$ (Mesh1)','$\Psi_c$ (Mesh2)','$\mathcal{K}$ (Mesh1)','$\mathcal{K}$ (Mesh2)'},...
     'NumColumns',3,'Location','northwest','interpreter','latex','FontSize', 24);
 set(leg1,'Box','off');

%--
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';%ax.YAxis.TickLabelFormat='%,.2f';
xlim([0 85.0]); 

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 20];fig.PaperPosition = [0. 0. 30 20];%fig.PaperPositionMode = 'auto';

fileName = 'branching-sigma2-random-energy'; 
print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');
