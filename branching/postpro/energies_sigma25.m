clc;clear all;close all;

iplot = 2;%1=std|2=marigo|3=pfczm

dtime      = 10^6; thickness = 1.0; %mm

    fileName   = 'sigma25/mesh0/branching-energies.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f'; A  = textscan(fileID, formatSpec, 'HeaderLines', 1); 
    fileName   = 'sigma25/mesh1/branching-energies.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f'; B  = textscan(fileID, formatSpec, 'HeaderLines', 1);
    fileName   = 'sigma25/mesh2/branching-energies.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; C  = textscan(fileID, formatSpec, 'HeaderLines', 1);

% convert cells to matrices
A          = cell2mat(A); B          = cell2mat(B);  C          = cell2mat(C);

% append zeros for first row
lb = size(A); lb2 = lb(2);
A = [zeros(1,lb2);A]; B = [zeros(1,lb2);B]; 

%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',1.5);set(0, 'DefaultAxesFontSize',20);

%----------------------------------------------------------------------
%---- Compare energies -------------------
%----------------------------------------------------------------------
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p11 = plot(dtime*A(:,1),A(:,3),'b-','LineWidth',1.5);
p12 = plot(dtime*B(:,1),B(:,3),'b--','LineWidth',1.5);
p13 = plot(dtime*C(:,1),C(:,3),'b-.','LineWidth',1.5);

p21 = plot(dtime*A(:,1),A(:,2),'r-','LineWidth',1.5);
p22 = plot(dtime*B(:,1),B(:,2),'r--','LineWidth',1.5);
p23 = plot(dtime*C(:,1),C(:,2),'r-.','LineWidth',1.5);

xlabel('Time [$\mu$s]','interpreter','latex','FontSize',20);
ylabel('Energy [mJ]','interpreter','latex','FontSize',20);

leg1=legend([ p11 p12 p13 p21 p22 p23],{'Stored energy ($h = b/2$)','Stored energy ($h = b/5$)','Stored energy ($h = b/10$)',...
    'Surface energy ($h = b/2$)','Surface energy ($h = b/5$)','Surface energy ($h = b/10$)'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 18);
 set(leg1,'Box','off');

%--
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';ax.YAxis.TickLabelFormat='%,.2f';
%xlim([0 80.0]); 

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

%fileName = 'branching-energy-sigma25'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); 
