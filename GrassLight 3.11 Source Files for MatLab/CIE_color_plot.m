function CIE_Data = CIE_color_plot(EdEu_Data, Biodis_Data)
%% Import data from spreadsheet
% Script for converting Rrs spectra to CIE simulated color
%       Created by R. Zimmerman  29 October 2024
%       Adapted from CIE code provided by Heidi Dierssen
%    
%   Requires Rrs values calculated by function EdEu_water_v3_for_app.m
%


% Convert to output type
Rrs = transpose(EdEu_Data.Rrs);
Rsfc = Rrs .* pi;
R_canopy = transpose(EdEu_Data.can_refl(Biodis_Data.n_flow_reduced_canopy_layers,:));

% Clear temporary variables
%clear opts tbl

%% Create the CIE functional response spectra

%load CIEint        this data file replaced by the following functions that provide the same wavelength specific functional responses
%CIE (Commission Internationale de l'Êclairage)
%CIE 1964 supplementary standard colorimetric observer
%An ideal observer whose color matching properties correspond to the CIE color matching functions for the 10º field size

CIEint(:,1) = transpose(400:1:700);         % Establish wavelengths for CIE color matching functions

% fourier coefficients for chromaticity "a" spectrum determined from 8th order fit to the CIE table data
 a_a0 =      0.3615;
 a_a1 =      0.3563;
 a_b1 =     0.05845;
 a_a2 =      0.2942;
 a_b2 =     0.08556;
 a_a3 =   -0.009145;
 a_b3 =     0.01099;
 a_a4 =     0.02644;
 a_b4 =     0.02761;
 a_a5 =     -0.0135;
 a_b5 =    0.007144;
 a_a6 =    0.004542;
 a_b6 =   -0.002896;
 a_a7 =   -0.001271;
 a_b7 =    0.003461;
 a_a8 =   -0.002759;
 a_b8 =   -0.004852;
 a_w =     0.02133;

 % 8th order Fourier function using the "a" coefficients above

CIEint(:,2) = a_a0 + a_a1 .* cos(CIEint(:,1) .* a_w) + a_b1 .* sin(CIEint(:,1) .* a_w) + a_a2 .* cos(2 .* CIEint(:,1) .* a_w) + a_b2 .* sin(2 .* CIEint(:,1) .* a_w) +...
    a_a3 .* cos(3 .* CIEint(:,1) .* a_w) + a_b3 .* sin(3 .* CIEint(:,1) .* a_w) + a_a4 .* cos(4 .* CIEint(:,1) .* a_w) + a_b4 .* sin(4 .* CIEint(:,1) .* a_w) +...
    a_a5 .* cos(5*CIEint(:,1) .* a_w) + a_b5 .* sin(5*CIEint(:,1) .* a_w) + a_a6 .* cos(6 .* CIEint(:,1) .* a_w) + a_b6 .* sin(6*CIEint(:,1) .* a_w) +...
    a_a7 .* cos(7*CIEint(:,1) .* a_w) + a_b7 .* sin(7*CIEint(:,1) .* a_w) + a_a8 .* cos(8*CIEint(:,1) .* a_w) + a_b8 .* sin(8*CIEint(:,1) .* a_w);


% Gaussian coefficients for chromaticity "b" spectrum determined from 7th order fir to the CIE table data
 b_a1 =      -3.342;
 b_b1 =       546.6; 
 b_c1 =       19.22;
 b_a2 =       2.343;
 b_b2 =         547;
 b_c2 =       18.17;
 b_a3 =    -0.00177;
 b_b3 =       568.6;
 b_c3 =       1.242;
 b_a4 =    -0.01136;
 b_b4 =       579.4;
 b_c4 =       3.661;
 b_a5 =   -0.009187;
 b_b5 =       548.5;
 b_c5 =       8.267;
 b_a6 =      0.9239;
 b_b6 =       562.7;
 b_c6 =       60.05;
 b_a7 =       1.148;
 b_b7 =       544.2;
 b_c7 =       23.07;

 % 7th order Gaussian function using the "b" coefficients above

CIEint(:,3) = b_a1 .* exp(-((CIEint(:,1)-b_b1) ./ b_c1) .^ 2) + b_a2 .* exp(-((CIEint(:,1)-b_b2) ./ b_c2) .^ 2) +... 
              b_a3 .* exp(-((CIEint(:,1)-b_b3) ./ b_c3) .^ 2) + b_a4 .* exp(-((CIEint(:,1)-b_b4) ./ b_c4) .^ 2) +...
              b_a5 .* exp(-((CIEint(:,1)-b_b5) ./ b_c5) .^ 2) + b_a6 .* exp(-((CIEint(:,1)-b_b6) ./ b_c6) .^ 2) +... 
              b_a7 .* exp(-((CIEint(:,1)-b_b7) ./ b_c7) .^ 2);


% Gaussian coeffients for chromaticity "c" spectrum determined from 3rd order fit to the CIE table data
c_a1 = 1.095;
c_b1 = 456.3;
c_c1 = 21.86;
c_a2 = 0.5689;
c_b2 = 462.5;
c_c2 = 41.61;
c_a3 = 0.8433;
c_b3 = 433.2;
c_c3 = 15.1;

% 3rd order Gaussina function using the "c" coefficients above
CIEint(:,4) = c_a1 .* exp(-((CIEint(:,1)-c_b1) ./ c_c1) .^ 2) + c_a2 .* exp(-((CIEint(:,1)-c_b2) ./ c_c2) .^ 2) + c_a3 .* exp(-((CIEint(:,1)-c_b3) ./ c_c3).^ 2);



delwv=1;                    % wavelength interval = 1 nm


xsub_Rrs=Rrs.*CIEint(:,2).*delwv;
ysub_Rrs=Rrs.*CIEint(:,3).*delwv;
zsub_Rrs=Rrs.*CIEint(:,4).*delwv;

xsub_Rsfc=Rsfc.*CIEint(:,2).*delwv;
ysub_Rsfc=Rsfc.*CIEint(:,3).*delwv;
zsub_Rsfc=Rsfc.*CIEint(:,4).*delwv;

xsub_Rcanopy=R_canopy.*CIEint(:,2).*delwv;
ysub_Rcanopy=R_canopy.*CIEint(:,3).*delwv;
zsub_Rcanopy=R_canopy.*CIEint(:,4).*delwv;

K_Rrs=100/sum(ysub_Rrs);
X_Rrs=K_Rrs*sum(xsub_Rrs);
Y_Rrs=K_Rrs*sum(ysub_Rrs);
Z_Rrs=K_Rrs*sum(zsub_Rrs);

K_Rsfc=100/sum(ysub_Rsfc);
X_Rsfc=K_Rsfc*sum(xsub_Rsfc);
Y_Rsfc=K_Rsfc*sum(ysub_Rsfc);
Z_Rsfc=K_Rsfc*sum(zsub_Rsfc);

K_Rcanopy=100./sum(ysub_Rcanopy);
X_Rcanopy=K_Rcanopy.*sum(xsub_Rcanopy);
Y_Rcanopy=K_Rcanopy.*sum(ysub_Rcanopy);
Z_Rcanopy=K_Rcanopy.*sum(zsub_Rcanopy);


Xsum_Rrs=X_Rrs+Y_Rrs+Z_Rrs;
x_Rrs=X_Rrs./Xsum_Rrs;
y_Rrs=Y_Rrs./Xsum_Rrs;
z_Rrs=Z_Rrs./Xsum_Rrs;

Xsum_Rsfc=X_Rsfc+Y_Rsfc+Z_Rsfc;
x_Rsfc=X_Rsfc./Xsum_Rsfc;
y_Rsfc=Y_Rsfc./Xsum_Rsfc;
z_Rsfc=Z_Rsfc./Xsum_Rsfc;

Xsum_Rcanopy=X_Rcanopy+Y_Rcanopy+Z_Rcanopy;
x_Rcanopy=X_Rcanopy./Xsum_Rcanopy;
y_Rcanopy=Y_Rcanopy./Xsum_Rcanopy;
z_Rcanopy=Z_Rcanopy./Xsum_Rcanopy;

%RGB values in a particular set of primaries can be transformed to and from CIE XYZ via a 3x3 matrix transform. 
%These transforms involve tristimulus values, that is a set of three linear-light components that conform to the 
%CIE color-matching functions. CIE XYZ is a special set of tristimulus values. In XYZ, any color is represented 
%as a set of positive values.
%
%To transform from XYZ to RGB (with D65 white point), the matrix transform used is [3]:

 mat =[  3.240479 -1.537150 -0.498535 ;
       -0.969256  1.875992  0.041556 ;
        0.055648 -0.204043  1.057311 ] ;


%This matrix has color in integers from 0 to 255.  Beware that sometimes you can get
%slightly negative or greater than 255.  Either reset to limits or change the K.
RGB_Rrs = mat*[X_Rrs Y_Rrs Z_Rrs]';
RGB_Rsfc = mat*[X_Rsfc Y_Rsfc Z_Rsfc]';
RGB_Rcanopy = mat*[X_Rcanopy Y_Rcanopy Z_Rcanopy]';


CIE_Data.matRGB_Rrs=transpose(RGB_Rrs./255);
CIE_Data.matRGB_Rrs(CIE_Data.matRGB_Rrs >= 1) = 1;
CIE_Data.matRGB_Rrs(CIE_Data.matRGB_Rrs <0) = 0;

CIE_Data.matRGB_Rsfc=transpose(RGB_Rsfc./255);
CIE_Data.matRGB_Rsfc(CIE_Data.matRGB_Rsfc >= 1) = 1;
CIE_Data.matRGB_Rsfc(CIE_Data.matRGB_Rsfc <0) = 0;

CIE_Data.matRGB_Rcanopy=transpose(RGB_Rcanopy./255);
CIE_Data.matRGB_Rcanopy(CIE_Data.matRGB_Rcanopy >= 1) = 1;
CIE_Data.matRGB_Rcanopy(CIE_Data.matRGB_Rcanopy <0) = 0;

% End of function CIE_color_plot