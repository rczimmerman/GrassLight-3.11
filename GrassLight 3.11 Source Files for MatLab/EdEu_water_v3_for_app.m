function EdEu_Data = EdEu_water_v3_for_app(Site_Data, IOP_Spectra, WQ_Params, Plant_Arch, Sed_Ref, LOPs, E0_sfc_spectra, E0_sfc_PAR_mubar, abKd_Spectra, Biodis_Data, PlantTypeDropDown, SubstrateTypeDropDown)

%% *************************************************************************************************************************************************************
% 
% Function to calculate the propagation of downwelling irradiance through the water column and submerged plant canopy of biomass distribution determined by function biodis.
%                                           written for GrassLight 3.0
%
% Adapted from Fortran biodis.for subroutine writte for GrassLight Ver 2.14 by RC Zimmerman 1 Jun 2024 
%
%
%   Usage: EdEu_Data = EdEu_water_v3_for_app(Site_Data, IOP_Spectra, WQ_Params, Plant_Arch, Sed_Ref, LOPs, E0_sfc_spectra, E0_sfc_PAR_mubar, abKd_Spectra, Biodis_Data, PlantTypeDropDown, SubstrateTypeDropDown)
%
%
% Other functions required:
%       read_GL3_input_data_for_app_302: reads Excel setup file containing data required to parameterize GrassLight
%       E0_surface_Ver3_for app:  calculates surface irradiances above and just below the water, sun angle and average cosine (in water)
%       abKd_Spectra_Ver1_for_app:  calculates diffuse attenuation spectrum for the water column; required to attenuate light with depth
%       biodis_v301_for_app: calculates vertical distribution of submerged plant biomass
%
% Input:  Tables of data read from setup file and calculated by other functions required for GL
%           Site_Data           = Site-relevant information (lat,lon, date, time, depth, plant type etc. stored in an excel file read by function read_GL3_input_data_for_app_302
%           IOP_Spectra         = IOP spectra for pure water and generic phytoplankton required to calculate attenuation coefficients read by function read_GL3_input_data_for_app_302
%           WQ_Params           = Water quality parameters and constituent concentrations read by function read_GL3_input_data_for_app_302
%           Plant_Arch          = Plant architecture values containing epiphyte load needed to attenuate light to the leaf surface read by function read_GL3_input_data_for_app_302
%           Sed_Ref             = Sediment reflectance spectra required to reflect light off the seafloor read by function read_GL3_input_data_for_app_302
%           LOPs                = Species-specific leaf optical properties - absorption and reflectance spectra read by function read_GL3_input_data_for_app_302
%           E0_Surface_spectra  = Surface irradiances calculated by function E0_surface_Ver3_for_app.
%           E0_Surface_PAR_mubar= Surface PAR (umol quanta m-2 s-1) and average cosine below the water at 0 depth
%           abKd_Spectra        = IOP and AOP spectra calcluated by abKd_Spectra_Ver1_for_app
%           Biodis_Data         = Vertical biomass distribution with depth calculated by function biodis_v301_for_app
%           PlantTypeDropDown   = Plant type information from app.
%           SubstrateTypeDropDown = Substrate type information from app.
%
% Output:   EdEu_Data           = Tables of data calculated by this function, including upwelling and downwelling irradiance & reflectance spectra,
%                                   water leaving radiance spectrum, and the light absorbed by the leaf in each layer.  
%                               
%%
%% *************************************************************************************************************************************************************


COPYRITE = 'Function EdEu_water_v3_for_app Copyright (c) 2024 by RC Zimmerman';
VER= 'Version 3';


%% Pre-allocate irradiance arrays as zeros

EdEu_Data.eow = zeros(Biodis_Data.nlayer, 301);
EdEu_Data.edw = zeros(Biodis_Data.nlayer, 301);                                  % Downwelling irradiance in water
EdEu_Data.euw = zeros(Biodis_Data.nlayer, 301);                                   % Upwelling irradiance in water
EdEu_Data.euw_down = zeros(Biodis_Data.nlayer, 301);
EdEu_Data.ed_leaf = zeros(Biodis_Data.nlayer, 301);                               % Downwelling irradiance incident on the leaf
EdEu_Data.eu_leaf = zeros(Biodis_Data.nlayer, 301);                               % Upwelling irradiance incident on the leaf
EdEu_Data.mubard = zeros(Biodis_Data.nlayer,1);                                  % Average cosine for downwelling irradiance
EdEu_Data.mubard(Biodis_Data.nlayer) = E0_sfc_PAR_mubar.mubard_0;                % mu_bar down at the surface of the water
EdEu_Data.piz = zeros(Biodis_Data.nlayer,1);                                     % Photosynthetically absorbed irradiance in each layer
EdEu_Data.piz_down = zeros(Biodis_Data.nlayer,1);

EdEu_Data.LU0_air = zeros(1, 301);                                                % Upwelling radiance, in air
EdEu_Data.Reflectance = zeros(1,301);                                             % refletance
EdEu_Data.Rrs = zeros(1,301);
EdEu_Data.can_refl = zeros(Biodis_Data.nlayer,301);

EdEu_Data.alpha = 1/10;                                  % quantum yield of photosyntheis for absorbed light  mol O2 (or C) per mol photons absorbed

%% Return with all values = 0 if at night
if E0_sfc_PAR_mubar.Sum_ED0_microMoles_air == 0 
    return
end

%%  Get Leaf Optical Properties
col = 1;
if PlantTypeDropDown.Value == "Zostera"
    labsl = LOPs.a_zostera;
    lrefl = LOPs.r_zostera;
    row = 1;
elseif PlantTypeDropDown.Value == "Thalassia"
    labsl = LOPs.a_thalassia;
    lrefl = LOPs.r_thalassia;
    row = 2;
elseif PlantTypeDropDown.Value == "Syringodium"
    labsl = LOPs.a_syringodium;
    lrefl = LOPs.r_syringodium;
    row = 3;
elseif PlantTypeDropDown.Value == "Phyllospadix"
    labsl = LOPs.a_phyllospadix;
    lrefl = LOPs.r_phyllospadix;
    row = 4;
elseif PlantTypeDropDown.Value == "Posidonia"
    labsl = LOPs.a_posidonia;
    lrefl = LOPs.r_posidonia;
    row = 5;
else 
    labsl = LOPs.a_ruppia;
    lrefl = LOPs.r_ruppia;
    row = 6;
end

%% Get Sediment Reflectance spectrum
 if SubstrateTypeDropDown.Value == "Bright Carbonate"
    sed_ref = Sed_Ref.Bright_Carbonate;
 elseif SubstrateTypeDropDown.Value == "Dark Carbonate"
    sed_ref = Sed_Ref.Dark_Carbonate;
 elseif SubstrateTypeDropDown.Value == "Silica Sand"
    sed_ref = Sed_Ref.Silica_Sand;
 else 
     sed_ref = Sed_Ref.Silica_Mud;
 end


%% Compute leaf aSTAR, photosynthetic absorptance and quantum yield of photosynthesis


aSTAR = transpose(labsl - labsl(length(labsl)));        % remove non-photosynthetic absorption from aSTAR calculation 
labsb = 1 - exp(-aSTAR);                                % photosynthetic absorption spectrum, per unit leaf


% Scale leaf total absorption and reflectance coefficients by projected leaf area within each layer
EdEu_Data.labsm = Biodis_Data.lap * transpose(labsl);       % leaf absorption in each layer, 1/m2 seafloor; increases with leaf density

lrefm = Biodis_Data.lap * transpose(lrefl);                 % leaf reflectance in each layer scaled to leaf area, this dot product generates a 1000 x 301 matrix                                                     

% Calculate light absorption coefficient of epiphytes
abs_epi = 0.3792;                                           % biomass-specific absorbance by leaf epiphytes from Bulthuis and Woelkerling 1983.  Constant for now but array here 
d_epi = Plant_Arch.epi(row,col)*abs_epi;                             % leaf-specific epiphyte absorbance
Absorptance_epi = 1-10.^(-d_epi);                           % convert absorbance to absorptance per unit leaf area
abs_coeff_epi = d_epi*2.303;                                % absorption coefficient for leaf epiphytes per unit leaf area
epi_abs_m = Biodis_Data.lap*abs_coeff_epi;                  % epiphyte absorption in each layer scaled to leaf area, generates a 1 x 100 matrix

% Calculate the waveband 
waveband = IOP_Spectra.WL(2,1) - IOP_Spectra.WL(1,1);       % waveband thickness needed to convert energy to quanta


%% Now, attenuate Ed down through the canopy
%       The constant 0.008351 converts energy (W/m^2/nm) to quanta (uE/m^2/s) per nm
%       The constant 0.5 is the average cosine for Lambertian scattering
EdEu_Data.edw(Biodis_Data.nlayer,:) = E0_sfc_spectra.ED0_water;                                                                                             % Downweling irradiance in the first layer; W/m2

for i = (Biodis_Data.nlayer):-1:2                                                                                                                           % layers increment from the bottom up; i = nlayer is the top of the canopy
    EdEu_Data.eow(i) = sum(EdEu_Data.edw(i,:)*IOP_Spectra.WL* 0.008351 * waveband);                                                                         % convert watts to quanta to get E_PAR, umol photons/m2/s
    EdEu_Data.euw_down(i,:) = EdEu_Data.edw(i,:) .* (transpose(1 - exp(-abKd_Spectra.BB * Biodis_Data.lathik)) + (lrefm(i,:)./ EdEu_Data.mubard(i,1)));     % first subtract the backscattered and reflected light, W/m2
    exponent_Kd = exp(-transpose(abKd_Spectra.Kd)*Biodis_Data.lathik);                                                                                      % attenuation by the water column in each layer, dimensionless
    exponent_leaf = exp(-EdEu_Data.labsm(i,:)/EdEu_Data.mubard(i));                                                                                         % attenuation by the leaf in each layer, dimensionless
    exponent_epiphyte = exp(-epi_abs_m(i,1)/EdEu_Data.mubard(i));                                                                                           % attenuation by epiphytes in each layer, dimensionless
    EdEu_Data.ed_leaf(i,:) = (EdEu_Data.edw(i,:) - EdEu_Data.euw(i,:)) .* exponent_epiphyte;                                                                % downwelling irradiance reaching leaf surface,  W/m2   
    EdEu_Data.piz_down(i) = sum(EdEu_Data.labsm(i) .* EdEu_Data.ed_leaf(i,:) *  IOP_Spectra.WL * 0.008351 * waveband);                                         % downwelling irradiance (quanta) absorbed by the leaf in layer i
    exponent_ed = exponent_Kd .* exponent_leaf .* exponent_epiphyte;                        % total downwelling light attenuation in each layer
    EdEu_Data.edw(i-1,:) = (EdEu_Data.edw(i,:) - EdEu_Data.euw(i,:)).*exponent_ed;                                      % downwelling irradiance in next layer i-1  
    EdEu_Data.mubard(i-1) = EdEu_Data.mubard(i) - (EdEu_Data.mubard(i) - 0.5) * Biodis_Data.lap(i);           % adjust mu_bar to become incrementally more isotropic as light propagates through the plant canopy
    EdEu_Data.Kd = abKd_Spectra.Kd.*EdEu_Data.mubard(i-1)./ EdEu_Data.mubard(i);
end

% and finish the bottom layer
EdEu_Data.eow(1) = sum(EdEu_Data.edw(1,:)*IOP_Spectra.WL* 0.008351 * waveband);                                         % convert watts to quanta to get E_PAR
EdEu_Data.euw_down(1,:) = EdEu_Data.edw(i,:) .* (transpose(1 - exp(-abKd_Spectra.BB * Biodis_Data.lathik)) + (lrefm(i,:)./ EdEu_Data.mubard(i,1)));                                            % first subtract the reflected light
exponent_Kd = exp(-transpose(abKd_Spectra.Kd)*Biodis_Data.lathik);                                               % attenuation by the water column in each layer
exponent_leaf = exp(-EdEu_Data.labsm(1,:)/EdEu_Data.mubard(1));                                                % attenuation by the leaf in each layer
exponent_epiphyte = exp(-epi_abs_m(1,1)/EdEu_Data.mubard(1));                                        % attenuation by epiphytes in each layer
EdEu_Data.ed_leaf(1,:) = (EdEu_Data.edw(1,:) - EdEu_Data.euw(1,:)) * exponent_epiphyte;           % downwelling irradiance reaching leaf surface
EdEu_Data.piz_down(1) = sum(EdEu_Data.labsm.* EdEu_Data.ed_leaf(1,:) .* EdEu_Data.labsm(1,:) *  IOP_Spectra.WL * 0.008351 * waveband);              % downwelling irradiance (quanta/m2 seafloor/s) absorbed by the leaf in layer i - not biomass normalized
exponent_ed = exponent_Kd .* exponent_leaf .* exponent_epiphyte;                        % total downwelling light attenuation in each layer
EdEu_Data.ed_bottom = (EdEu_Data.edw(1,:) - EdEu_Data.euw(1,:)).*exponent_ed;                                      % downwelling irradiance reaching the bottom  

%% Now reflect Edw off the bottom
mubaru = 0.5;                                               % average cosine for upwelling irradiance, which is assumed to be lambertian
EdEu_Data.Ku = EdEu_Data.Kd .* EdEu_Data.mubard(1,1)/mubaru;                              % Value of Ku, assuming all upwweellingg iirradiance is lambertian
EdEu_Data.euw(1,:) = EdEu_Data.euw_down(1,:) + (EdEu_Data.ed_bottom(1,:) .* transpose(sed_ref(:,1)));    % add bottom reflected light to backscattered light in the lowest layer
eu_ref = EdEu_Data.euw(1,:) .* (lrefm(1,:)./ EdEu_Data.mubard(1,1));           % Lose some Eu reflected back downward by the leaf canopy (ignore secondary reflectances)
exponent_Ku = exp(-transpose(EdEu_Data.Ku) * Biodis_Data.lathik);                                               % attenuation by the water column in each layer
exponent_leaf = exp(-EdEu_Data.labsm(1,:)/EdEu_Data.mubard(1));                                                % attenuation by the leaf in each layer
exponent_epiphyte = exp(-epi_abs_m(1,1)/EdEu_Data.mubard(1));                                        % attenuation by epiphytes in each layer
EdEu_Data.eu_leaf(1,:) = (EdEu_Data.euw(1,:) - eu_ref) * exponent_epiphyte;                      % upwelling irradiance reaching leaf surface
EdEu_Data.piz(1) = EdEu_Data.piz_down(1) + sum(EdEu_Data.labsm.* EdEu_Data.eu_leaf(1,:) .* EdEu_Data.labsm(i,:)*  IOP_Spectra.WL * 0.008351 * waveband);              % downwelling irradiance (quanta) absorbed by the leaf in layer i
exponent_eu = exponent_Ku .* exponent_leaf .* exponent_epiphyte;                        % total upwnwelling light attenuation in bottom-most layer
EdEu_Data.euw(2,:) = EdEu_Data.euw_down(2,:) + EdEu_Data.euw(1,:) .* exponent_eu;                      % attenuate euw(1,:) to make euw(2,:)
EdEu_Data.eow(1) = EdEu_Data.eow(1) + sum(EdEu_Data.euw(1,:)*IOP_Spectra.WL*0.008351* waveband);

% Now, propagate Eu upward through the canopy and determine the amount of Eu absorbed in each layer on the way up.
for i = 2:Biodis_Data.nlayer                                                                          % layers increment from the bottom up; i = 1 is the base of the canopy
    euref = EdEu_Data.euw(i,:).*lrefm(i,:)./ mubaru;                                                % downward reflectance of upwelling light; a loss term only
    exponent_Ku = exp(-transpose(EdEu_Data.Ku) .* Biodis_Data.lathik);                                               % attenuation by the water column in each layer
    exponent_leaf = exp(-EdEu_Data.labsm(i,:) ./ mubaru);                                                % attenuation by the leaf in each layer
    exponent_epiphyte = exp(-epi_abs_m(i,1) ./ mubaru);                                        % attenuation by epiphytes in each layer
    exponent_eu = exponent_Ku .* exponent_leaf .* exponent_epiphyte;                        % total downwelling light attenuation in each layer
%   euw(i,:) = (euw(i,:) + (euw(i-1,:) - euref)).*exponent_eu;                              % upwelling irradiance in layer i
    EdEu_Data.euw(i,:) = EdEu_Data.euw_down(i,:) + (EdEu_Data.euw(i-1,:) - euref).*exponent_eu; 
    EdEu_Data.eu_leaf(i,:) = (EdEu_Data.euw(i,:) + EdEu_Data.euw(i-1,:) - euref) .* exponent_Kd .* exponent_epiphyte;                                               % downwelling irradiance reaching leaf surface
    EdEu_Data.piz(i) = EdEu_Data.piz_down(i) + EdEu_Data.labsm(i) .* sum(EdEu_Data.eu_leaf(i,:) .* EdEu_Data.labsm(i,:) *IOP_Spectra.WL* 0.008351* waveband);          % downwelling irradiance (quanta) absorbed by the leaf in layer i
    EdEu_Data.eow(i) = EdEu_Data.eow(i) + sum(EdEu_Data.euw(i,:)*IOP_Spectra.WL*0.008351* waveband);
end

tau_water_to_air = 0.55;                                                                            % water to air transmittance, Mobley 1994
Q = pi;                                                                                             % Ratio of Eu/Lu at the sea surface; set to pi assuming Lambertian bottom
EdEu_Data.can_refl = EdEu_Data.euw ./ EdEu_Data.edw;                                                % canopy reflectance in each layer
EdEu_Data.EU0_air = EdEu_Data.euw(Biodis_Data.nlayer,:) .* tau_water_to_air;                        % upwelling irradiance in air
EdEu_Data.LU0_air = EdEu_Data.EU0_air/pi;                                                           % upwelling radiance in air
EdEu_Data.Reflectance = (EdEu_Data.EU0_air ./ transpose(E0_sfc_spectra.ED0_air))/Q;                 % Surface reflectance, in air
EdEu_Data.Rrs = EdEu_Data.Reflectance/pi;                                                           % Surface Rrs

% End of EdEu function





