clc;clear all;close all;

dtime      = 10^6; formatSpec = '%f %f %f %f'; thickness = 1.0; %mm

fileName = 'loading2/wu3-energies.dat'; fileID     = fopen(fileName,'r'); A          = textscan(fileID, formatSpec, 'HeaderLines', 2);

% convert cells to matrices
A          = cell2mat(A); 

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

p1 = plot(dtime*A(:,1),A(:,3),'b-','LineWidth',1.5);
p2 = plot(dtime*A(:,1),A(:,2),'r-','LineWidth',1.5);

ylabel('Strain/Surface Energy [mJ]','interpreter','latex','FontSize',30);
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30);
ylim([0 0.35]);

% yyaxis right;
% p3 = plot(dtime*A(:,1),A(:,4),'k-','LineWidth',1.5);
% q3 = plot(dtime*B(:,1),B(:,4),'k--','LineWidth',1.5);

% % r3 = plot(dtime*C(:,1),C(:,4)/1000,'k-.','LineWidth',1.5);
% ylabel('Kinetic Energy [mJ]','interpreter','latex','FontSize',30);
% ylim([0 15]);

% %crack branching
patch('vertices', [19.5, 0; 20.3, 0; 20.3, 0.3; 19.5 0.3],'faces', [1, 2, 3, 4],'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 

leg1=legend([ p1 p2 ],{'$\Psi_s$','$\Psi_c$'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 24);
 set(leg1,'Box','off');

%--
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';%ax.YAxis.TickLabelFormat='%,.2f';
xlim([0 42.0]); 

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 20];fig.PaperPosition = [0. 0. 30 20];%fig.PaperPositionMode = 'auto';

fileName = 'bobaru-L2b-energy-wave';
%print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');
