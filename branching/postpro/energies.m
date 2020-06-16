clc;clear all;close all;

%iplot = 2;%1=std|2=marigo|3=pfczm
for iplot = 1:3

dtime      = 10^6; formatSpec = '%f %f %f'; thickness = 1.0; %mm

if iplot == 1 
    fileName   = 'sigma1b05/branching-std-energies.dat'; fileID     = fopen(fileName,'r'); A  = textscan(fileID, formatSpec, 'HeaderLines', 1); 
    fileName   = 'sigma1b025/branching-std-energies.dat'; fileID     = fopen(fileName,'r'); B  = textscan(fileID, formatSpec, 'HeaderLines', 1);
end

if iplot == 2 
    fileName   = 'sigma1b05/branching-marigo-energies.dat'; fileID     = fopen(fileName,'r'); A  = textscan(fileID, formatSpec, 'HeaderLines', 1); 
    fileName   = 'sigma1b025/branching-marigo-energies.dat'; fileID     = fopen(fileName,'r'); B  = textscan(fileID, formatSpec, 'HeaderLines', 1);
end

if iplot == 3 
    fileName   = 'sigma1b05/branching-implicit-energies.dat'; fileID     = fopen(fileName,'r'); A  = textscan(fileID, formatSpec, 'HeaderLines', 1); 
    fileName   = 'sigma1b025/branching-implicit-energies.dat'; fileID     = fopen(fileName,'r'); B  = textscan(fileID, formatSpec, 'HeaderLines', 1);
end

% convert cells to matrices
A          = cell2mat(A); B          = cell2mat(B);

% append zeros for first row
lb = size(A); lb2 = lb(2);
A = [zeros(1,lb2);A]; B = [zeros(1,lb2);B]; 

%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',2.0);set(0, 'DefaultAxesFontSize',30);

%----------------------------------------------------------------------
%---- Compare energies -------------------
%----------------------------------------------------------------------
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

p11 = plot(dtime*A(:,1),A(:,3),'b-');
p12 = plot(dtime*B(:,1),B(:,3),'b--');

p21 = plot(dtime*A(:,1),A(:,2),'r-');
p22 = plot(dtime*B(:,1),B(:,2),'r--');

xlabel('Time [$\mu$s]','interpreter','latex','FontSize',30);
ylabel('Energy [mJ]','interpreter','latex','FontSize',30);

leg1=legend([ p11 p12 p21 p22],{'$\Psi_s$ ($b = 0.50$ mm)','$\Psi_s$ ($b = 0.25$ mm)','$\Psi_c$ ($b = 0.50$ mm)','$\Psi_c$ ($b = 0.25$ mm)'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 24);
 set(leg1,'Box','off');

%--
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';ax.YAxis.TickLabelFormat='%,.2f';
xlim([0 80.0]); 

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

if iplot == 1; fileName = 'branching-energy-std'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); end;
if iplot == 2; fileName = 'branching-energy-marigo'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); end;
if iplot == 3; fileName = 'branching-energy-wu'; print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); end;
end