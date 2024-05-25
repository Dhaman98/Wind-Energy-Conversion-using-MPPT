%% Microgrid Simulation
% Created: ECE 530 class of fall 2023

clc
close all   % close figure windows
clear
format compact



%% Simulation settings

simu.endTime = 60*60*24;
simu.maxStepSize = 1e-1;



%% Data

load windspeedtimeseries
load loadtimeseries
load illuminationcurrenttimeseries.mat

meanWindSpeed = 10;
meanLoad = 2500e3;

% Scaling
windspeedtimeseries.Data = meanWindSpeed*windspeedtimeseries.Data;
loadtimeseries.Data = meanLoad*loadtimeseries.Data;



%% Wind Turbine Parameters

wt.rho = 1.2; % Density of air kg/m^3

% XANT 100 kW
wt.Pgen_rated = 10000e3;
wt.B = 0.01*100e3/6^2;    % Nm/(rad/s)    B*6^2 = 0.01*100e3
wt.bladeLength = 11;
wt.bladeWeight = 1000;    % From http://windpower.sandia.gov/other/041605C.pdf
wt.J = 3*1/3*wt.bladeWeight*wt.bladeLength^2;  % Missing generator, gearbox, shaft, etc...

wt.A = pi*wt.bladeLength^2;

wt.bladeActuationTau = 0.5;    % Time constant for response of hydraulic system that positions blades to desired angle

wt.w_0 = 0;

% Cp curve modeling
% lambdaai = 1/( (1/(lambda+0.08*beta)) - 0.035/(beta^3+1) )
% cp = c1*(c2/lambdaai-c3*beta-c4)*exp(-c5/lambdaai) + c6*lambda
plot_cp = 0;
if plot_cp
    wt.c = [
        0.5176
        116
        0.4
        5
        21
        0.0068
        ];
    figure
    lambda = [0:0.1:13];
    for beta = [0:5:30];
        lambdaai = 1./(1./(lambda+0.08*beta) - 0.035./(beta.^3+1));
        cp = wt.c(1)*(wt.c(2)./lambdaai - wt.c(3)*beta - wt.c(4)) .* exp(-wt.c(5)./lambdaai) + wt.c(6)*lambda;
        hold on
        plot(lambda,cp)
        hold off
    end
    axis([0 13 0 0.5])
    xlabel('lambda (tip speed ratio)')
    ylabel('Cp')
end

% From Cp(lambda) plot
wt.lambda_opt = 8.1;
wt.Cp_max = 0.48;

% Region 2 and 3 boundary -> rated rotational speed and wind speed
wt.u_rated = (wt.Pgen_rated/(wt.Cp_max*0.5*1.2*wt.A))^(1/3);   % P_rated = Cp_max*0.5*A*bladeLength*u_rated^3
wt.w_rated = wt.lambda_opt*wt.u_rated/wt.bladeLength;

% Speed controller
% Turbine rated torque ~ 16.7 kNm
wt.speedcontroller.kp = 16.7e3/7*2*2*2;       % Maps speed error to torque reference
wt.speedcontroller.ki = (16.7e3/7)/10*2*2*2;  % Maps speed error to rate of change of torque reference
wt.speedcontroller.lowerLimit = 0;      % No fan
wt.speedcontroller.upperLimit = 5*17e3;
wt.speedcontroller.kt = 1;
wt.speedcontroller.int_0 = 0;

% Power controller
wt.powercontroller.kp = 30/100e3;       % Maps power error to blade pitch reference
wt.powercontroller.ki = 30/100e3/2.5;   % Maps power error to rate of change of blade pitch reference
wt.powercontroller.int_0 = 0;
wt.powercontroller.kt = 1;
wt.powercontroller.upperLimit = 30;    % Maximum degrees of pitch
wt.powercontroller.lowerLimit = 0;     % Minimum degrees of pitch



%% Energy Storage

% PE 95% eff
% Battery 90% round trip efficient

es.Prated = 124939e3;

es.etape = 0.95;
es.etasm = sqrt(0.9);
% es.selfDischarge = 0;

es.Erated_hrs = 24;
es.Erated_kWh = es.Prated/1000*es.Erated_hrs;
es.Erated = es.Erated_kWh*(1000*3600);

es.Ppeupper = es.Prated;
es.Ppelower = -es.Prated;

es.petau = 1;
es.Ppe_0 = 80e3;

es.SOC_0 = 0.95;

% Energy storage cost estimate
energyStorageCost_dollarsperkW = 100;
energyStorageCost_dollarsperkWh = 130;
energyStorageCost_power = energyStorageCost_dollarsperkW*es.Prated/1000;
energyStorageCost_energy = energyStorageCost_dollarsperkWh*es.Erated/(3600*1000);
energyStorageCost = energyStorageCost_power + energyStorageCost_energy;



%% Microgrid

mg.H = 5;   
mg.D = 0.001;
mg.Pbase = 4583e3;
mg.fpu_0 = 1;
mg.X = 0.05;    % PU impedance in connection to main grid

es.Kgprim = 20;
es.Kgsec = es.Kgprim/(60*10);



%% Solar

% Parameters based on SolarWorld 300 module
% 300 W module, 60 cells, 5 W per cell
% Open circuit voltage of 40 V (0.66 V per cell)
% Short circuit current of 9.83 A
% Maximum power point at 32.6 V (0.54 V per cell) and 9.31 A

pv.Is = 1e-10; % Produces about  0.66 V at 9.8 A of current
pv.Rs = 0.005; % Adjusted this to get the approximate max power point
pv.Rp = 2500;  % Based off of L. Ma et al. "The Measurement of Series and Shunt Resistances of the Silicon Solar Cell Based on LabVIEW"
pv.VT = 0.026;

pv.MPPT_sampleTime = 1;  % Not optimized

pv.Prated = 13347e3;