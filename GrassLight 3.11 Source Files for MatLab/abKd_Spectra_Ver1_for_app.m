function abKd_Spectra = abKd_Spectra_Ver1_for_app(WQ_Params, IOP_Constants, IOP_Spectra, E0_sfc_PAR_mubar)
% *********************************************************************************
% Matlab function abKd_Spectra_Ver1_for app.m
% Created by Richard C. Zimmerman, based on relations presented by
% Lee, Z.-P., K. Du, and R. Arnone. 2005. A model for the diffuse attenuation coefficient of downwelling irradiance. J. Geophys. Res 110: 1-10.
%
%  Data required from Excel set-up file (some may als be changed in the main app by user):
%       WQ_Params - CDOM absorbance @ 440 nm, Chl and TSM concentrations
%       IOP_Constants - values from Lee et al. 2005
%  
%   Data required from other Matlab functions
%       none
%
%   Data passed to other functions
%       abKd_Spectra - table of absorption, scattering and attenuation coefficients, per m
%
% **********************************************************************************

COPYRITE = 'Copyright (c) 2024 by RC Zimmerman';
VER= 'Version 1';

%%   Calculate water column absorption, scattering, and attenuation components (originally peformed by SpecKd which is no longer needed)

 WL = IOP_Spectra.WL;                                                                                   % Wavelength
 THETADEG = E0_sfc_PAR_mubar.THETAZ_Direct_air*180./pi;                                                 % In-air solar angle, (deg.) as required by Lee et al. (2005
 AG = WQ_Params.AG_440*exp(-IOP_Constants.SG.*(WL-440));                                                % CDOM absorption spectrum
 ACHL = (WQ_Params.CHL_A*IOP_Constants.APHST675).*IOP_Spectra.ASTCHL;                                   % Chlorophyll absorption spectrum
 ANAP = WQ_Params.TSM*(IOP_Constants.BLTRB+IOP_Constants.SIGATRB).*exp(-IOP_Constants.STRB*(WL-440));   % NAP absorption spectrum
 AT = IOP_Spectra.AW+AG+ACHL+ANAP;                                                                      % Total absorption
 
 BP = (WQ_Params.TSM.*IOP_Constants.SIGBTRB).*(555./IOP_Spectra.WL).^IOP_Constants.ETA;                 % Particulate scattering coefficient
 BB = 0.5*IOP_Spectra.BW+IOP_Constants.BB2B.*BP;                                                        % Total backscattering (bbp + bbw)
 
 Kd_aw = (1.+0.005*THETADEG).*IOP_Spectra.AW;                                                           % Kd due to absorption by water
 Kd_ag = (1.+0.005*THETADEG).*AG;                                                                       % Kd due to absorption by CDOM
 Kd_achl = (1.+0.005*THETADEG).*ACHL;                                                                   % Kd due to absorption by phytoplankton
 Kd_anap = (1.+0.005*THETADEG).*ANAP;                                                                   % Kd due to absorption by non-algal particles
 Kd_atotal = Kd_aw + Kd_ag + Kd_achl + Kd_anap;                                                         % Kd due to absorption by all optical components
 Kd_bb = 4.18*(1.-0.52.*exp(-10.8.*AT)).*BB;                                                            % Kd due to backscattering
 Kd = Kd_atotal + Kd_bb;                                                                                % Total Kd

 % Put all calculated values into a table for passing to other functions
 abKd_Spectra = table(WL, AG, ACHL, ANAP, AT, BP, BB, Kd_aw, Kd_ag, Kd_achl, Kd_anap, Kd_atotal, Kd_bb, Kd);

 %% End function abKdSpectra_Ver_1_for_app
