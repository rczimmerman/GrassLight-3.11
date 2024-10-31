function [E0_sfc_spectra, E0_sfc_PAR_mubar] = E0_surface_Ver3_for_app(Site_Data, Atm_Const, Atm_IOPs)
%     *******************************************************************
%           Function E0_surface Ver 3.0 for GrassLignt Ver 3.0
%
%           Derived from fortran subroutine E05nm Ver 2.14  written for GrassLight 2.41
%
%         Computes from astronomic data provided by the user:
%             1) downwelling irradiance spectrum just below the sea surface
%                 at 1 nm resolution
%             2) photoperiod
%
%         Usage: [E0_sfc_spectra, E0_sfc_PAR_mubar] = E0_surface_Ver3_for_app(Site_Data, Atm_Const, Atm_IOPs)
%
%         Required Inputs:
%               Site_Data:  table of site data information read from excel file or manually entered via GUI
%               Atm_Const:  table of atmospheric constants read from excel file
%               Atm_IOPs:   table of atmospheric IOPs read from excel file
%
%         Provided Outputs:
%             E0_sfc_spectra:   table of wavelength values and surface irradiances in air and in water
%             E0_sfc_PAR_mubar: table of individual values for Daylength, sun angle, surface PAR above & below water, in water average cosine at the surface
%
%     *********************************************************************

COPYRITE = 'Copyright (c) 2024 by RC Zimmerman';
VER= 'Version 3';

%% Calculate incident irradiance spectrum (W m^-2 nm^-1) just below the water surface.  Implements Gregg and Carder (1990)
   
WL = transpose(400:1:700);  %define wavelehgth range and precision (nm)

%   First section calculates constants used to determine sun position
ECC=-0.0167;                                % Solar Eccentricity
AVOGADRO = 6.02e23;                         % Avogadro's Number - particles per mole
GAMMA=2*pi*(Site_Data.JD-1)/365.;                     % Date converted to polar coordinates
DEL=(0.006918-0.399912*cos(GAMMA)+0.070257*sin(GAMMA)-0.006758*cos(2*GAMMA)+0.000907*sin(2*GAMMA)-0.002697*cos(3*GAMMA)+0.00148*sin(3*GAMMA));
ET=(0.000075+0.001868*cos(GAMMA)-0.032077*sin(GAMMA)-0.014615*cos(2*GAMMA)-0.04089*sin(2*GAMMA))*229.18;
RPHI=pi*Site_Data.XLAT/180.;                          %Latitude: convert degrees to radians            
RLONG=pi*Site_Data.XLON/180.;                         %Longitude: convert degrees to radians 
OMEGASR=(acos(-tan(RPHI)*tan(DEL)));        % Angle of sunrise
OMEGASS=-OMEGASR;                           % Angle of sunset
ND=(180./pi)*2.*OMEGASR/15.;               % Hours of sunlight
HSR=12.-ND/2.;                              % Time of sunrise
HSS=12.+ND/2.;                              % Time of sunset
daylen=ND;                                  % Photoperiod

if Site_Data.HR<HSR || Site_Data.HR>HSS                         % Turn the lights off between sunset and sunrise
    EDS(1:1:301,1) = 0;
    EDD(1:1:301,1) = 0;
    ED(1:1:301,1) = 0;
    EDD_air(1:1:301,1) = 0;
    EDS_air(1:1:301,1) = 0;
    ED0_air(1:1:301,1) = 0;
    ED0_water(1:1:301,1) = 0;
    Sum_ED0_microMoles_air = 0;
    Sum_ED0_microMoles_water = 0;
    mubard_0 = 0;
    THETAZ_Direct_air = 0;
    % Combine Arrays Into Output Tables for use by the app other subroutines and the app
    E0_sfc_spectra = table(EDD_air, EDS_air, ED0_air, ED0_water);
    E0_sfc_PAR_mubar = table(daylen, THETAZ_Direct_air, Sum_ED0_microMoles_air, Sum_ED0_microMoles_water, mubard_0);
    return
end
    OMEGA=12.*OMEGASR/(12.-HSR)+(OMEGASR/(HSR-12.))*Site_Data.HR;                                   % Sun angle based on time of day
    THETAZ_Direct_air=(acos(sin(DEL)*sin(RPHI)+cos(DEL)*cos(RPHI)*cos(OMEGA)));                     % Direct Solar zenith angle base in radians
    THETAZ_Diffuse_air = 1.067198;                                                                  % Diffuse Solar zenith angle for assuming sky irradiance is lambertian
    MTHETA=1/(cos(THETAZ_Direct_air)+0.50572*(96.07995-180.*THETAZ_Direct_air/pi)^-1.6364);         % Slant path through the atmosphere for all gases but ozone
    MOZTHETA=1.0035/sqrt((cos(THETAZ_Direct_air))^2.+0.007);                                        % Slant path through the atmosphere for ozone
    MPTHETA=MTHETA*Atm_Const.PRESS/1013.2;                                                          % Slant path corrected for atmospheric pressure (at sea surface)

%   Ozone climatology section.  Default coefficients to N. America
%   Coefficients from Table 1 in Van Heuklon (1979, Solar Energy, p. 66)
      if Site_Data.XLAT >= 0.% Northern Hemisphere
        OZ_A=150.0;           
        OZ_BETA=1.28;
        OZ_C=40.0;
        OZ_F=-30.0;
        OZ_G=20.0;
        OZ_H=3.0;             % Ozone Scale Height, cm
        if Site_Data.XLON>=0.            
           OZ_I=20.0;       % Eastern Hemisphere
        else 
           OZ_I=0.0;         % Western Hemisphere
        end
      else                    % Southern Hemisphere
        OZ_A=100.0;          
        OZ_BETA=1.5;
        OZ_C=30.0;
        OZ_F=152.625;
        OZ_G=20.0;
        OZ_H=2.0;
        OZ_I=-75.0;
      end

      OZONE=(235+(OZ_A+OZ_C*sin(0.9856*(Site_Data.JD+OZ_F)*pi/180.)+OZ_G*sin(OZ_H*pi*(Site_Data.XLON+OZ_I)/180.))*(sin(OZ_BETA*RPHI))^2)/1000.;
%
%   Compute aerosol size distribution
	F=((2.-Atm_Const.RH/100.)/(6.*(1-Atm_Const.RH/100.)))^(1./3.);          
	A1COEF=2000*Atm_Const.AM^2;
	Y(1)=log10(A1COEF)-log10((0.1/(F*0.03))^2.)-log10(F);
	A2COEF=5.866*(Atm_Const.WM-2.2);
	Y(2)=log10(A2COEF)-log10((1/(F*0.24))^2.)-log10(F);
	A3COEF=0.01527*(Atm_Const.W-2.2)*0.05;
	Y(3)=log10(A3COEF)-log10((10/(F*2.))^2.)-log10(F);
	X(1)=log10(0.1);
	X(2)=0.;
	X(3)=log10(10.);
%    
%   Get aerosol Junge coef. and amplitude from slope and intercept
    p = polyfit(X,Y,1);
	CAMPL=exp(p(1));
    GAMAER = p(2);
	CA=3.91/Atm_Const.V;
	ALPHAE0=-(GAMAER+3.);               % Angstrom formulation for aerosol optical thickness
	TAUA550=CA*Atm_Const.HA;
	BETA=TAUA550/(0.55^-ALPHAE0);
    if ALPHAE0< 0 
	  AVCOSTHETA = 0.82;
    elseif ALPHAE0>1.2
        AVCOSTHETA = 0.65;
    else 
	    AVCOSTHETA = -0.1417*ALPHAE0+0.82;
    end
	B3COEF = log10(1-AVCOSTHETA); 
	B2COEF = B3COEF*(0.0783+B3COEF*(-0.3824-0.5874*B3COEF));
	B1COEF = B3COEF*(1.459+B3COEF*(0.1595+0.4129*B3COEF));
	FA = 1-0.5*exp((B1COEF+B2COEF*cos(THETAZ_Direct_air))*cos(THETAZ_Direct_air));
	OMEGA_A = (-0.0032*Atm_Const.AM+0.972)*exp(0.000306*Atm_Const.RH);
%    
%   Reflectance section
    if Atm_Const.W>7 
	  CD=(0.49+0.065*Atm_Const.W)*0.001;  %Drag coefficient
    else 
	  CD=(0.62+1.56*Atm_Const.W)*0.001;
    end
	RHOA=1200.;
    if Atm_Const.W>7
	  RHOF=(0.000045*RHOA*CD-0.00004)*Atm_Const.W^2.;
    elseif Atm_Const.W<4
	  RHOF=0.;
    else
	  RHOF=0.000022*RHOA*CD*Atm_Const.W^2.-0.0004;
    end
	RHODSP=0.0253*exp((-0.000714*Atm_Const.W+0.0618)*((180.*THETAZ_Direct_air/pi)-40.));
    if Atm_Const.W>4.
	  RHOSSP=0.057;
    else
	  RHOSSP=0.066;
    end
	RHOD=RHODSP+RHOF;
	RHOS=RHOSSP+RHOF;
%    
%   Calculate of flux and transmittance arrays at 1 nm intervals at the sea surface

    F0 = Atm_IOPs.H0.*(1+ECC*cos(2.*pi.*(Site_Data.JD-3)/365.)).^2;                         % Extraterrestrial solar irradiance (top of atmosphere)
	TR = exp(-(MPTHETA)./(0.0000000001156406.*WL.^4-0.000001335.*WL.^2));                   % Total Rayleign scattering coefficient
    TOZ = exp(-Atm_IOPs.AOZ.*OZONE*MOZTHETA);                                               % Ozone transmittance
    TOX = exp(-((1.41.*Atm_IOPs.AOX.*MPTHETA)./(1+118.3.*Atm_IOPs.AOX.*MPTHETA).^0.45));    % Oxygen transmittance
    TW = exp(-((0.238.*Atm_IOPs.AW_INC.*Atm_Const.WV.*MTHETA)./(1+20.07.*Atm_IOPs.AW_INC.*Atm_Const.WV.*MTHETA).*0.45));  % Water transmittance
    TAUA = BETA.*(WL./1000.).^(-ALPHAE0);                                                   % Aerosol optical thickness
    TA = exp(-TAUA.*MTHETA);                                                                % Transmittance after aerosol scattering and absorption
    TAA = exp(-(1.-OMEGA_A).*TAUA.*MTHETA);                                                 % Transmittance after aerosol absorption (no scattering)
    TAS = exp (-OMEGA_A.*TAUA.*MTHETA);                                                     % Transmittance after aerosol scattering (no absorption)
    IR = F0.*cos(THETAZ_Direct_air).*TOZ.*TOX.*TW.*TAA.*(1-TR.^0.95).*0.5;                  % Rayleigh scattering
    IA = F0.*cos(THETAZ_Direct_air).*TOZ.*TOX.*TW.*TAA.*TR.^1.5.*(1-TAS).*FA;               % Aerosol scattering
    EDD_air = F0.*cos(THETAZ_Direct_air).*TR.*TA.*TOZ.*TOX.*TW.*(1-RHOD);                   % Direct irradiance
    EDS_air =(IR + IA).*(1-RHOS);                                                           % Diffues(sky) Irradiance
    ED0_air = EDD_air +EDS_air;                                                             % Total Downwelling irradiance at the surface
    Sum_ED0_Watts_air = sum(ED0_air);                                                       % Total Downwelling irradiance summed across wavelentth

    %Convert above water light energy to umol quanta and sum across wavelength
    ED0_microMoles_air = ED0_air.*WL.*0.008351;
    Sum_ED0_microMoles_air =sum(ED0_microMoles_air);

    % Transmit the direct light across the air-water interace
    sin_THETA_Direct_water = sin(THETAZ_Direct_air)/1.341;                                      % refract the direct light across the air-sea interface; 1.341 is the refractive index of seawater
    THETA_Direct_water = asin(sin_THETA_Direct_water);
    Spec_reflectance_direct = (1/2)*((sin(THETAZ_Direct_air - THETAZ_Direct_air))^2/(2*(sin(THETAZ_Direct_air + THETAZ_Direct_air))^2) + ((tan(THETAZ_Direct_air - THETAZ_Direct_air))^2/(tan(THETAZ_Direct_air + THETAZ_Direct_air))^2));   % Snell's Law
    EDD_water = EDD_air .* (1 - Spec_reflectance_direct);
    cos_THETA_Direct_water = cos(THETA_Direct_water);                                             % cosine of the direct beam in water

    % transmit the diffuse light across the air-water interface
    EDS_water = EDS_air .* (1-0.066);                               % average reflectance loss of diffuse irradiance at the sea surface
    sin_THETA_Diffuse_water = sin(deg2rad(60))/1.341;               % refract the diffuse light across the air-sea interface; assuming a lambertian sky. 1.341 is the refractive index of seawater
    THETA_Diffuse_water = asin(sin_THETA_Diffuse_water);
    cos_THETA_Diffuse_water = cos(THETA_Diffuse_water);             % Cosine of the direct beam in water

    % sum the direct, diffuse and total fluxes, all units in W m^-2
    ED0_water = EDD_water + EDS_water;
    Sum_EDD_water = sum(EDD_water);
    Sum_EDS_water = sum(EDS_water);
    Sum_ED0_water = sum(ED0_water);

    % compute the average cosine for downwelling irradinace at the surface
    frac_EDS = Sum_EDS_water/Sum_ED0_water;
    mubard_0 = cos_THETA_Direct_water*(1-frac_EDS) + cos_THETA_Diffuse_water*(frac_EDS);

    % Convert submarine light energy to micromoles quanta and sum to get total micromoles quanta 
    ED0_microMoles_water = EDS_water.*WL.*0.008351;                     % Spectral diffuse watts to umol quanta
    Sum_ED0_microMoles_water = sum(ED0_microMoles_water); 

    %% Combine Arrays Into Output Tables for use by the app other subroutines and the app

    E0_sfc_spectra = table(WL, EDD_air, EDS_air, ED0_air, ED0_water);
    E0_sfc_PAR_mubar = table(daylen, THETAZ_Direct_air, Sum_ED0_microMoles_air, Sum_ED0_microMoles_water, mubard_0);
    
    % End of Function E0_surface_Ver3_for_app
