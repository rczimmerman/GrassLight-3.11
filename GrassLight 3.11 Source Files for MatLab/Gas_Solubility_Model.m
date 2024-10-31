function Gas_Sat = Gas_Solubility_Model(Temp, Sal, Depth, pCO2_Air)

%% Gas Solubility Model using equations from Garcia and Gordon (1992) as corrected and reported by Pilson 1998
%
%   Ver 1 written by R Zimmerman 4 Jan 2023
%   Required input parameters:
%   Temp - water temperature (°C)
%   Sal - salinity (ppt or PSS)
%   Depth - water depth (m)
%

COPYRITE = 'Copyright (c) 2024 by RC Zimmerman';
VER= 'Version 1';


%% O2 polynomial coefficients for temperature and salinity
A0_O2 = 5.80871;
A1_O2 = 3.20291;
A2_O2 = 4.17787;
A3_O2 = 5.10006;
A4_O2 = -0.00987;
A5_O2 = 3.80369;

B0_O2 = -7.02e-3;
B1_O2 = -7.70e-3;
B2_O2 = -1.14e-3;
B3_O2 = -9.52e-3;

C0_O2 = -2.76e-7;

%% N2 polynomial coefficients for temperature and salinity
A1_N2 = -173.2221;
A2_N2 = 254.6078;
A3_N2 = 146.3611;
A4_N2 = -22.0933;

B1_N2 = -0.05405;
B2_N2 = 0.027266;
B3_N2 = -0.00384;

%% Scaled temperatures
T_Kelvin = Temp + 273.15;                % Water temperature in Kelvin
Ts = log(298.15 - Temp) / (273.15+Temp);    % scaled tempreature

%% Gas solubility at 1 Atm

ln_O2_sol = A0_O2 + (A1_O2*Ts) + (A2_O2*Ts^2) + (A3_O2*Ts^3) + (A4_O2*Ts^4) + (A5_O2*Ts^5) + (Sal*(B0_O2 + (B1_O2*Ts)) + (B2_O2*Ts^2) + (B3_O2*Ts^4)) + (C0_O2*Ts^2);
O2_1atm = exp(ln_O2_sol);    % µmol/Kg

ln_N2_sol = A1_N2 + (A2_N2*100/T_Kelvin) + (A3_N2*log(T_Kelvin/100)) + (A4_N2*T_Kelvin/100) + (Sal *(B1_N2+B2_N2*(T_Kelvin/100)+(B3_N2*(T_Kelvin/100)^2)));  
N2_1atm =exp(ln_N2_sol);        % µmol/Kg

% ln_H_CO2 = (9345.17/T_Kelvin) - 167.8108 + 23.3585*log(T_Kelvin) + (0.023517 - 0.00023656*T_Kelvin + 0.00000047036*T_Kelvin^2) * sal;
ln_H_CO2 = (9345.17/T_Kelvin) - 167.8108 + 23.3585*log(T_Kelvin) + (0.023517 - 0.00023656*T_Kelvin + 0.00000047036*T_Kelvin^2) * Sal;
CO2_1atm = exp(ln_H_CO2) * pCO2_Air;    % µmol/Kg
pCO2_1atm = CO2_1atm/exp(ln_H_CO2);     % ppm

%% Gas solubility at depth (µmol Kg^-1) 

O2_Depth = O2_1atm * (Depth/10+1);             % µmol/Kg, (depth/10 + 1) accounts for the increased pressure (Atm), assuming 1 atm at 0 depth
N2_Depth = N2_1atm * (Depth/10+1);             % µmol/Kg
CO2_Depth = CO2_1atm * (Depth/10+1);           % µmol/Kg
pCO2_Depth = CO2_Depth/exp(ln_H_CO2);   % ppm

%% Combine arrays into a table
Gas_Sat = table(O2_1atm, O2_Depth, N2_1atm, N2_Depth, CO2_1atm, CO2_Depth, pCO2_1atm, pCO2_Depth);

% End of function Gas Solubility Model



