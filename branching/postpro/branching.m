% dynamic branching
clc; clear all; close all;

E = 3.2e4; nu = 0.2; ft = 12; Gf = 3e-3; rho = 2450*10^(-12);%dynamic-branching

mu = E/(2*(1+nu)); lambda = E*nu/((1+nu)*(1-2*nu));
Cd = sqrt((lambda+2*mu)/rho);
Cs = sqrt(mu/rho); 
Cr = (0.862 + 1.14* nu)/(1+nu)* Cs; % Rayleigh wave speed %Cr = 2119*10^3;
Cr

hmin = 0.1; dt_cfl = 0.8*hmin/Cs

hmin = 0.04; dt_cfl = 0.8*hmin/Cs

hmin = 0.019; dt_cfl = 0.8*hmin/Cs