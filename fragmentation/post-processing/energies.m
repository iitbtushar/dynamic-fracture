clc;clear all;close all;

thickness = 70.0; %mm
dtime      = 10^6; %convert to micro second
dEnergy      = 10^3; %convery to kJ

iplot = 31; %1=>ft = 800/100;|| 2 => mesh; || 21=>mesh-Marigo|| 22=>mesh-Std || 3=>length scale ||31=> l0-marigo ||32=> l0-std 

%--ft ----
if iplot == 1
fileName   = 'mesh1/wu800.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; A  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh1/wu.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; B  = textscan(fileID, formatSpec, 'HeaderLines', 1);

char1 = 'Hirmand-StrainEnergy'; C = xlsread('cylinder-refs.xlsx',char1);
char1 = 'Hirmand-SurfaceEnergy'; D = xlsread('cylinder-refs.xlsx',char1);
end

%-- mesh -----
if iplot == 2
fileName   = 'mesh1/wu.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; A  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh2/wu.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; B  = textscan(fileID, formatSpec, 'HeaderLines', 1);
end

if iplot == 21
fileName   = 'mesh1/marigo.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; A  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh2/marigo.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; B  = textscan(fileID, formatSpec, 'HeaderLines', 1);
end

if iplot == 22
fileName   = 'mesh1/std.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; A  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh2/std.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; B  = textscan(fileID, formatSpec, 'HeaderLines', 1);
end

%-- length scale ----
if iplot == 3
fileName   = 'mesh0/wu.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; A  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh1/wu-b15.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; B  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh1/wu.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; C  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh2/wu-l0.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; D  = textscan(fileID, formatSpec, 'HeaderLines', 1);
end

if iplot == 31
fileName   = 'mesh0/marigo.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; A  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh1/marigo.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; B  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh1/marigo-b1.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; C  = textscan(fileID, formatSpec, 'HeaderLines', 1);
fileName   = 'mesh2/marigo-l0.dat'; fileID     = fopen(fileName,'r'); formatSpec = '%f %f %f %f'; D  = textscan(fileID, formatSpec, 'HeaderLines', 1);
end

% convert cells to matrices
A          = cell2mat(A); B          = cell2mat(B);  C          = cell2mat(C); D          = cell2mat(D);



%---- figure defaults ---
set(0, 'defaultAxesTickLabelInterpreter','latex');set(0, 'defaultLegendInterpreter',       'latex');
set(0, 'defaultlinelinewidth',1.5);set(0, 'DefaultAxesFontSize',30);

%----------------------------------------------------------------------
%---- Compare load deformation -------------------
%----------------------------------------------------------------------
figure(1); clf; hold on; set(gcf, 'Position', get(0, 'Screensize'));

if iplot ==1
p11 = plot(dtime*A(:,1),A(:,3)*(thickness/1000)/dEnergy,'b-','LineWidth',1.5);
p12 = plot(dtime*B(:,1),B(:,3)*(thickness/1000)/dEnergy,'b--','LineWidth',1.5);
p13 = plot(C(:,1),C(:,2),'b-.','LineWidth',1.5);

p21 = plot(dtime*A(:,1),A(:,2)*(thickness/1000)/dEnergy,'r-','LineWidth',1.5);
p22 = plot(dtime*B(:,1),B(:,2)*(thickness/1000)/dEnergy,'r--','LineWidth',1.5);
p23 = plot(D(:,1),D(:,2),'r-.','LineWidth',1.5);

elseif (iplot==3)||(iplot==31)
p11 = plot(dtime*A(:,1),A(:,3)*(thickness/1000)/dEnergy,'b-','LineWidth',1.5);
p12 = plot(dtime*B(:,1),B(:,3)*(thickness/1000)/dEnergy,'b--','LineWidth',1.5);
p13 = plot(dtime*C(:,1),C(:,3)*(thickness/1000)/dEnergy,'b-.','LineWidth',1.5);
%p14 = plot(dtime*D(:,1),D(:,3)*(thickness/1000)/dEnergy,'bo','LineWidth',1.5);

p21 = plot(dtime*A(:,1),A(:,2)*(thickness/1000)/dEnergy,'r-','LineWidth',1.5);
p22 = plot(dtime*B(:,1),B(:,2)*(thickness/1000)/dEnergy,'r--','LineWidth',1.5);
p23 = plot(dtime*C(:,1),C(:,2)*(thickness/1000)/dEnergy,'r-.','LineWidth',1.5);
%p24 = plot(dtime*D(:,1),D(:,2)*(thickness/1000)/dEnergy,'ro','LineWidth',1.5);

else
p11 = plot(dtime*A(:,1),A(:,3)*(thickness/1000)/dEnergy,'b-','LineWidth',1.5);
p12 = plot(dtime*B(:,1),B(:,3)*(thickness/1000)/dEnergy,'b--','LineWidth',1.5);

p21 = plot(dtime*A(:,1),A(:,2)*(thickness/1000)/dEnergy,'r-','LineWidth',1.5);
p22 = plot(dtime*B(:,1),B(:,2)*(thickness/1000)/dEnergy,'r--','LineWidth',1.5);
end


%crack branching
% patch('vertices', [31, 0; 33, 0; 33, 0.3; 31 0.3], ...
%           'faces', [1, 2, 3, 4],'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);

xlabel('time [$\mu$s]','interpreter','latex','FontSize',30);
ylabel('Energy [kJ]','interpreter','latex','FontSize',30);

if iplot == 1
leg1=legend([ p11 p12 p21 p22 ],{'Stored energy ($f_t = 800$ MPa)','Stored energy ($f_t = 1000$ MPa)',...
    'Surface energy ($f_t = 800$ MPa)','Surface energy ($f_t = 1000$ MPa)'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 30);
end

if (iplot == 2) || (iplot == 21) || (iplot == 22)
leg1=legend([ p11 p12 p21 p22 ],{'stored energy (Mesh1)','stored energy (Mesh2)',...
    'surface energy (Mesh1)','surface energy (Mesh2)'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 30);
end

if (iplot == 3)||(iplot == 31)
leg1=legend([ p11 p12 p13 p21 p22 p23 ],{'$\Psi_s$ ($b = 2.0$ mm)','$\Psi_s$ ($b = 1.5$ mm)','$\Psi_s$ ($b = 1.0$ mm)',...
    '$\Psi_c$ ($b = 2.0$ mm)','$\Psi_c$ ($b = 1.5$ mm)','$\Psi_c$ ($b = 1.0$ mm)'},...
     'NumColumns',2,'Location','northwest','interpreter','latex','FontSize', 20);
end

set(leg1,'Box','off');
 

%--
set(gca,'XMinorTick','on','YMinorTick','on');grid('on');box('on');
ax=gca;ax.XAxis.TickLabelFormat='%,.0f';ax.YAxis.TickLabelFormat='%,.1f';
if (iplot == 1) || (iplot == 2) || (iplot == 3); xlim([0 100.0]); ylim([0 6.5]); end
if (iplot == 21)|| (iplot == 31); xlim([0 90.0]); ylim([0 6.5]); end
if (iplot == 22); xlim([0 123.0]); ylim([0 6]); end

%---
fig = gcf;fig.PaperUnits = 'centimeters';fig.PaperType='<custom>';fig.PaperSize=[30 15];fig.PaperPosition = [0. 0. 30 15];%fig.PaperPositionMode = 'auto';

if iplot == 1; fileName = 'cylinder-energy-ft'; end;
if iplot == 2; fileName = 'cylinder-energy-mesh'; end; 
if iplot == 21; fileName = 'cylinder-energy-mesh-marigo'; end;
if iplot == 22; fileName = 'cylinder-energy-mesh-std'; end;
if iplot == 3; fileName = 'cylinder-energy-l0'; end;
if iplot == 31; fileName = 'cylinder-energy-l0-marigo'; end;
print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000'); print(fileName,'-painters','-dpdf','-r1000');
