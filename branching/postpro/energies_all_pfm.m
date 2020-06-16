clc;clear all;close all;

iplot = 2; %1=all pfm | 2 Harmind(2019)

dtime      = 10^6; formatSpec = '%f %f %f'; thickness = 1.0; %mm


    fileName   = 'sigma1b025/branching-std-energies.dat'; fileID     = fopen(fileName,'r'); A  = textscan(fileID, formatSpec, 'HeaderLines', 1);
    fileName   = 'sigma1b025/branching-marigo-energies.dat'; fileID     = fopen(fileName,'r'); B  = textscan(fileID, formatSpec, 'HeaderLines', 1); 
    fileName   = 'sigma1b025/branching-implicit-energies.dat'; fileID     = fopen(fileName,'r'); C  = textscan(fileID, formatSpec, 'HeaderLines', 1);
    
    E = load('Hirmand2019/hirmand-fracture-energy.out');
    D = load('Hirmand2019/hirmand-strain-energy.out');

% convert cells to matrices
A          = cell2mat(A); B          = cell2mat(B); C          = cell2mat(C);

% append zeros for first row
lb = size(A); lb2 = lb(2);
A = [zeros(1,lb2);A]; B = [zeros(1,lb2);B]; C = [zeros(1,lb2);C];

%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',2.0);set(0, 'DefaultAxesFontSize',30);

%----------------------------------------------------------------------
%---- Compare load deformation -------------------
%----------------------------------------------------------------------
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));
if iplot ==1
p11 = plot(dtime*A(:,1),A(:,3),'b-.');
p12 = plot(dtime*B(:,1),B(:,3),'b--');
p13 = plot(dtime*C(:,1),C(:,3),'b');

p21 = plot(dtime*A(:,1),A(:,2),'r-.');
p22 = plot(dtime*B(:,1),B(:,2),'r--');
p23 = plot(dtime*C(:,1),C(:,2),'r');
end

if iplot ==2
p13 = plot(dtime*C(:,1),C(:,3),'b');
p14 = plot(D(:,1),D(:,2),'b+');

p23 = plot(dtime*C(:,1),C(:,2),'r');
p24 = plot(E(:,1),E(:,2),'r+');
end

%crack branching
if iplot ==1
patch('vertices', [31, 0; 33, 0; 33, 0.3; 31 0.3], ...
          'faces', [1, 2, 3, 4],'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
end

xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30);
ylabel('Energy [mJ]','interpreter','latex','FontSize',30);

if iplot ==1
leg1=legend([ p13 p12 p11 p23 p22 p21 ],{'$\Psi_s$ (\texttt{PF-CZM})','$\Psi_s$ (\texttt{AT1})','$\Psi_s$ (\texttt{AT2})',...
    '$\Psi_c$ (\texttt{PF-CZM})','$\Psi_c$ (\texttt{AT1})','$\Psi_c$ (\texttt{AT2})'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 20);
end
if iplot ==2
leg1=legend([ p13 p14 p23 p24 ],{'$\Psi_s$ (\texttt{PF-CZM})','$\Psi_s$ Hirmand et. al(2019)',...
    '$\Psi_c$ (\texttt{PF-CZM})','$\Psi_c$ Hirmand et. al(2019)'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 24);
end
 set(leg1,'Box','off');
 

%--
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';ax.YAxis.TickLabelFormat='%,.2f';
xlim([0 80.0]); 

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

if iplot ==1; fileName = 'branching-energy-all'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');end;
if iplot ==2; fileName = 'branching-energy-compare'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');end;