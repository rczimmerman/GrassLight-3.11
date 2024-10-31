function Prod_Data = produce_v301_for_app(WQ_Params, Site_Data, Plant_Arch, Photo_Params, Biodis_Data, EdEu_Data, E0_sfc_PAR_mubar, PlantTypeDropDown)
%	*******************************************************************************
%                       Matlab Function produce_v3_for_app
% Originally derived from GL 2.14 FORTRAN subroutine  produce_214.for.  
% Calculates productivity of a seagrass canopy  for GrassLight Version 3, a stand-alone application written in Matlab
%
%                       Author: Richard C. Zimmerman
%                       Original Ver. 3 Release Date   1 June 2024
%
%   Data required from Excel set-up file (some may als be changed in the main app by user):
%       WQ_Params - water quality parameters, nutrient concentrations and carbonate chemistry values
%       Site_Data - location, date, time, water depth etc.  
%       Plant_Arch - plant architecture information (height, shoot:rhizome:root fraction, shoot density, sigmoid biomass parameters etc.)
%       Photo_Params - metabolic rate parameters specific for each species
%    
%   Data required from othe Matlab functions:
%       Biodis_Data - Biomass distribution data calculated by Matlab function biodis_v3_for app.m
%       EdEu_Data - Water column irradiance data calculated by Matlab function EdEu_v3_for_app.m
%       E0_sfc_PAR_mubar - surface irradiance data calculated by Matlab function E0_surface_Ver3_for_app.m
%		  
%   Data required from Main App
%       PlantTypeDropDown - plant type specified by set-up file or changed by user
%
%
%
%
%
%                             Update Log
%		Date			     Modification                       Ver	  	By
%	---------------  -----------------------------------       ------ -------
%
%
%	********************************************************************************

% tst = 'in Produce';       % debugging stmt
% disp(tst)

    COPYRITE = 'Copyright (c) 2024 by RC Zimmerman';
    VER= 'Version 3';

%% Initialize arrays and accumulators as zeros - necessry for runs with 0 biomass
    Prod_Data.Pm = 0;
    Prod_Data.PE = 0;
    Prod_Data.Mean_Biomass_normalized_P_vs_E = 0;
    Prod_Data.P_now_per_shoot = 0;
    Prod_Data.P_now_over_PE = 0;
    Prod_Data.Instantaneous_Leaf_P_to_R = 0;
    Prod_Data.Instantaneous_Shoot_P_to_R = 0;
    Prod_Data.Daily_P_to_R = 0;
    Prod_Data.piz_per_leaf = zeros(Biodis_Data.nlayer,1);
    Prod_Data.Biomass_normalized_P_vs_E = zeros(Biodis_Data.nlayer,1);
    Prod_Data.Daily_Biomass_normalized_P_vs_E = zeros(Biodis_Data.nlayer,1);
    Prod_Data.Daily_Biomass_normalized_P = zeros(Biodis_Data.nlayer,1);
    Prod_Data.Daily_P_per_shoot = 0;

    Prod_Data.Leaf_Daily_R =  0;
    Prod_Data.Root_Daily_R = 0;
    Prod_Data.Rhiz_Daily_R = 0;
    Prod_Data.Plant_Daily_R = 0;
    Prod_Data.del13C_new = Inf;

    Prod_Data.leaf_biomass_per_shoot = 0;
    Prod_Data.total_biomass_per_shoot = 0;
    Prod_Data.rhiz_biomass_per_shoot = 0;
    Prod_Data.root_biomass_per_shoot = 0;
    Prod_Data.leaf_carbon_per_shoot = 0;
    Prod_Data.rhiz_carbon_per_shoot = 0;
    Prod_Data.root_carbon_per_shoot = 0;
    Prod_Data.total_carbon_per_shoot = 0;

    Prod_Data.Pf_CO2 = 0;                               
    Prod_Data.Up_CO2 = 0;      
    Prod_Data.PE_CO2 = 0;

    Prod_Data.Pf_HCO3 = 0;                          
    Prod_Data.Up_HCO3 = 0;
    Prod_Data.PE_HCO3 = 0;

    Prod_Data.PE = 0;   
    
    Prod_Data.Leaf_Daily_R = 0;
    Prod_Data.Root_Daily_R = 0;
    Prod_Data.Rhiz_Daily_R = 0;
    Prod_Data.Plant_Daily_R = 0;
    Prod_Data.Daily_P_to_R = Inf;

    Prod_Data.Leaf_Resp_per_shoot =0;
    Prod_Data.Rhiz_Resp_per_shoot = 0;
    Prod_Data.Root_Resp_per_shoot = 0;
    Prod_Data.Shoot_Resp = 0;
    Prod_Data.Daily_Hsat_Requirement = 0;
    
%% Select Metabolic Parameters to be used based on species selected.  Row refers to the corresponding row in the Photo_Params and Plant_Arch Tables.  All column values are 1 for each parameter in the table since they have names
            col = 1;
        if PlantTypeDropDown.Value == "Zostera"
            row = 1;
        elseif PlantTypeDropDown.Value == "Thalassia"
            row = 2;
        elseif PlantTypeDropDown.Value == "Syringodium"
            row = 3;
        elseif PlantTypeDropDown.Value == "Phyllospadix"
            row = 4;
        elseif PlantTypeDropDown.Value == "Posidonia"
            row = 5;
        else 
            row = 6;
        end

%% If no plants, skip to end of the function

    if Plant_Arch.density(row,col) > 0          % Proceed with these calcs only if there are plants
    %% Set current velocity at leaf boundary
        leaf_curvel = 5;         % cm sec-1 at leaf boundary, independent of water current velocity specified in site data
    
    %% Comvert plant size estimates from m2 leaf to g DW and to g C
        Prod_Data.leaf_biomass_per_shoot = Plant_Arch.areash(row,col) .* 500;    % 500 g DW/m^2 leaf converts m^2 shoot^-1 to g DW shoot^-1
        Prod_Data.total_biomass_per_shoot = Prod_Data.leaf_biomass_per_shoot ./ Plant_Arch.shoot_frac(row,col);
        Prod_Data.rhiz_biomass_per_shoot = Prod_Data.total_biomass_per_shoot .* Plant_Arch.rhiz_frac(row,col);
        Prod_Data.root_biomass_per_shoot = Prod_Data.total_biomass_per_shoot .* Plant_Arch.root_frac(row,col);
        Prod_Data.leaf_carbon_per_shoot = Prod_Data.leaf_biomass_per_shoot .* 0.34;        % Assumes seagrass biomass is 34% carbon by wt
        Prod_Data.rhiz_carbon_per_shoot = Prod_Data.rhiz_biomass_per_shoot .* 0.34;
        Prod_Data.root_carbon_per_shoot = Prod_Data.root_biomass_per_shoot .* 0.34;
        Prod_Data.total_carbon_per_shoot = Prod_Data.total_biomass_per_shoot .* 0.34;
    
    %% Adjust photosynthetic parameters for temperature
        Q10 = 3;                                    % Q10 coefficient for plant metabolism
        Q10_effect = Q10^((WQ_Params.Temp - 21)/10);          % Relative effect on metabolic rate, scaled to 21 °C
        Pm_CO2 = Q10_effect * Photo_Params.Pm_CO2(row,col);               % Effect of temperature on Pm_CO2
        Pm_HCO3 = Q10_effect * Photo_Params.Pm_HCO3(row,col);             % Effect of temperature on Pm_HCO3
        Prod_Data.Pm = Pm_CO2 + Pm_HCO3;          % Resulting physiological Pmax
    
    %% Compute degree of CO2- and flow-limited photosynthesis according to McPherson-Zimmerman photosynthesis model ; all rates area umol O2/CO2 m^2 leaf s^-1
        % CO2 first
        Prod_Data.Pf_CO2 = Photo_Params.Pm_CO2(row,col)*WQ_Params.CO2/(Photo_Params.Ks_CO2(row,col) + WQ_Params.CO2);                               % Flow saturated photosynthesis dependent on [CO2] using Michaelis Menten function
        Prod_Data.Up_CO2 = Photo_Params.Upf_CO2(row,col) * (1-exp(-Photo_Params.beta_CO2(row,col) * leaf_curvel / Photo_Params.Upf_CO2(row,col))) + Photo_Params.Up0_CO2(row,col);      % CO2 permeability function
        Prod_Data.PE_CO2 = (1/2)*((Photo_Params.Ks_CO2(row,col) * Prod_Data.Up_CO2 + WQ_Params.CO2*Prod_Data.Up_CO2 + Prod_Data.Pf_CO2) - sqrt((Photo_Params.Ks_CO2(row,col) * Prod_Data.Up_CO2 + WQ_Params.CO2 * Prod_Data.Up_CO2 + Prod_Data.Pf_CO2)^2 - 4 * WQ_Params.CO2 * Prod_Data.Up_CO2 * Prod_Data.Pf_CO2));   % Hill-Whittingham equation for CO2
    
        % Now HCO3
        Prod_Data.Pf_HCO3 = Photo_Params.Pm_HCO3(row,col) * WQ_Params.HCO3 / (Photo_Params.Ks_HCO3(row,col) + WQ_Params.HCO3);                          % Flow saturated photosynthesis dependent on [CO2] and [HCO3] using Michaelis Menten function
        Prod_Data.Up_HCO3 = Photo_Params.Upf_HCO3(row,col) * (1-exp(-Photo_Params.beta_HCO3(row,col) * leaf_curvel / Photo_Params.Upf_HCO3(row,col))) + Photo_Params.Up0_HCO3(row,col); % HCO3 permeability function
        Prod_Data.PE_HCO3 = (1/2)*((Photo_Params.Ks_HCO3(row,col) * Prod_Data.Up_HCO3 + WQ_Params.HCO3 * Prod_Data.Up_HCO3 + Prod_Data.Pf_HCO3)-sqrt((Photo_Params.Ks_HCO3(row,col) * Prod_Data.Up_HCO3 + WQ_Params.HCO3 * Prod_Data.Up_HCO3 + Prod_Data.Pf_HCO3)^2-4*WQ_Params.HCO3 * Prod_Data.Up_HCO3 * Prod_Data.Pf_HCO3));   % Hill-Whittingham equation for HCO3
    
        % Now sum CO2 + HCO3
        Prod_Data.PE = Prod_Data.PE_CO2 + Prod_Data.PE_HCO3;                                          % Light-saturated, flow-dependent photoynthesis per m2 leaf area
    
    %% Compute light-limited photosythesis
        % right now
        Prod_Data.piz_per_leaf = EdEu_Data.piz ./ Biodis_Data.lai;                              % Light absorbed per m2 leaf
        Prod_Data.Biomass_normalized_P_vs_E = Prod_Data.PE .* (1 - exp(-EdEu_Data.alpha .* Prod_Data.piz_per_leaf ./ Prod_Data.PE));               % Light dependent photosyntheis in each layer, µmol O2/m2 leaf/s in each layer
        Prod_Data.Mean_Biomass_normalized_P_vs_E = mean(Prod_Data.Biomass_normalized_P_vs_E, "omitnan");                  % Mean light dependent photosynthesis of the whole canopy, µmol O2/m2 leaf/s, scaled to the relative biomass distribution
        Prod_Data.P_now_per_shoot = sum(Prod_Data.Biomass_normalized_P_vs_E .* Biodis_Data.lai,"omitnan")/ Plant_Arch.density(row,col);                     % Integrated shoot photosynthesis, µmol O2/shoot/s
        Prod_Data.P_now_over_PE = Prod_Data.Mean_Biomass_normalized_P_vs_E/Prod_Data.PE;                                                           % Canopy photosythesis relative to PE, dimensionless 
    
    
    %%  Now compute respiration rates per shoot per second
        Prod_Data.Leaf_Resp_per_shoot = Q10_effect * Photo_Params.Leaf_R_21_deg(row,col) .* Plant_Arch.areash(row,col);                            % Leaf R, µmol O2 (or C) per shoot per sec
        Prod_Data.Root_Resp_per_shoot = Prod_Data.Leaf_Resp_per_shoot / 2 .* Plant_Arch.root_frac(row,col);                               % Root R
        Prod_Data.Rhiz_Resp_per_shoot = Prod_Data.Root_Resp_per_shoot / 2 .* Plant_Arch.rhiz_frac(row,col);                               % Rhizome R
        Prod_Data.Shoot_Resp = Prod_Data.Leaf_Resp_per_shoot + Prod_Data.Root_Resp_per_shoot + Prod_Data.Rhiz_Resp_per_shoot;   % Total shoot R
    
    %% Instantaneous P:R
        Prod_Data.Instantaneous_Leaf_P_to_R = Prod_Data.P_now_per_shoot / Prod_Data.Leaf_Resp_per_shoot;
        Prod_Data.Instantaneous_Shoot_P_to_R = Prod_Data.P_now_per_shoot / Prod_Data.Shoot_Resp;
           
    %% Now compute d13C value for newly produced tissue
        Prod_Data.del13C_new = 28 * Prod_Data.P_now_over_PE - 28;                            % d13C value for light-saturated photosynthesis



%% Compute daily metabolic integrals assuming time = noon   
        WL = transpose(400:1:700);
    
        Prod_Data.Sum_Daily_Biomass_normalized_P = 0;       % set daily P accumulator to 0, final units are umol O2 (or C) per shoot per day
        if Site_Data.HR == 12                               % if simulation performed at noon 
            tod_inc = E0_sfc_PAR_mubar.daylen * 6;          % 10 minute increments per day for daily integrations
            time_inc = 1/tod_inc;                           % time increment, fractional days of sunlight
            
            if mod(tod_inc,E0_sfc_PAR_mubar.daylen) ~= 0    % see if there is a remainder
                tod_inc = tod_inc + 1;                      % add 1 to counter if remainder
            end
           
            time_of_day = 0;                                % set time of day counter to 0
            Prod_Data.Sum_Daily_Biomass_normalized_P = 0;   % set daily biomass P accumulator to 0
            Prod_Data.Sum_Daily_P_per_shoot = 0;            % set daily shoot P accumulator to 0
    
            for j = 1: tod_inc
                time_of_day = time_of_day + time_inc;       % increment time of day, for computation of fractional day of sunlight
                light_frac = sin(pi*time_of_day);           % fraction of noon irradiance, based on time of day
                sun_angle = asin(light_frac) - pi/2;        % sun angle derived from time of day
                light_correction = cos(sun_angle);          % cosine correction for sun angle
                Prod_Data.Daily_Biomass_normalized_P = (Prod_Data.PE .* (1 - exp(-(EdEu_Data.alpha * light_correction) .* Prod_Data.piz_per_leaf ./ Prod_Data.PE))) .* 600;    % µmol O2/m2 leaf/10 min
                Prod_Data.Daily_P_per_shoot = Prod_Data.Daily_Biomass_normalized_P .* Biodis_Data.lai / Plant_Arch.density(row,col);               % umol O2/shoot/10 min
                Prod_Data.Sum_Daily_Biomass_normalized_P = Prod_Data.Sum_Daily_Biomass_normalized_P + sum(Prod_Data.Daily_Biomass_normalized_P,"omitnan");            % µmol O2/m2 m2 leaf/day
                Prod_Data.Sum_Daily_P_per_shoot = Prod_Data.Sum_Daily_P_per_shoot + sum(Prod_Data.Daily_P_per_shoot,"omitnan");    
           end
            %Prod_Data.Daily_P_per_shoot = Plant_Arch.areash(row,col) .* Prod_Data.Sum_Daily_Biomass_normalized_P;                                               % µmol O2/shoot/day
        else
            Prod_Data.Daily_Biomass_normalized_P = inf;
            Prod_Data.Daily_P_per_shoot = inf;              % Set daily production estimate to "inf" if Site_Data.HR ~= 12
        end
  
%%  Now compute daily respiration rate, P:R per shoot and daily Hsat requirement
        if Site_Data.HR == 12
            Prod_Data.Leaf_Daily_R = Prod_Data.Leaf_Resp_per_shoot * 60 * 60 * 24;       % scaling up from seconds to days  umol O2 (or C) per shoot per day
            Prod_Data.Root_Daily_R = (Prod_Data.Root_Resp_per_shoot * 60 * 60 * E0_sfc_PAR_mubar.daylen) + (0.65 * Prod_Data.Root_Resp_per_shoot * 60 * 60 * (24 - E0_sfc_PAR_mubar.daylen));
            Prod_Data.Rhiz_Daily_R = (Prod_Data.Rhiz_Resp_per_shoot * 60 * 60 * E0_sfc_PAR_mubar.daylen) + (0.65 * Prod_Data.Rhiz_Resp_per_shoot * 60 * 60 * (24 - E0_sfc_PAR_mubar.daylen));
            Prod_Data.Plant_Daily_R = (Prod_Data.Leaf_Daily_R + Prod_Data.Root_Daily_R + Prod_Data.Rhiz_Daily_R);
            Prod_Data.Daily_P_to_R = Prod_Data.Sum_Daily_P_per_shoot / Prod_Data.Plant_Daily_R;
            Prod_Data.Daily_Hsat_Requirement = Prod_Data.Plant_Daily_R /(Prod_Data.PE * Plant_Arch.areash(row,col) * 60 * 60);
        else
           Prod_Data.Daily_P_to_R = Inf;
           Prod_Data.Daily_Hsat_Requirement = Inf;
        end

    end
end

%% End of Produce function
