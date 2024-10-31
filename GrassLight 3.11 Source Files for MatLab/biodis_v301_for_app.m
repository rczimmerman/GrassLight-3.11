function [Biodis_Data, Plant_Arch] = biodis_v301_for_app(Site_Data, Plant_Arch, PlantTypeDropDown)
% ***********************************************************************************************************************
%%              Matlab Function  biodis.m
% Calculates the fractional biomass distribution, leaf area index (LAI) and
% vertically projected leaf area (LAP) for up to 1000 evenly-spaced canopy
% layers, depending on water depth and canopy height. Adjusts for canopy bending as a result of current flow.
% If plant canopy is taller than water depth, the emergent biomass is evenly distributed into the upper 10% of the water column, creating a floating surface
% layer
%
%   Revised and updated for Matlab GrassLight Ver 3 by R. Zimmerman by incorporating all changes in FORTRAN subroutine biodis_215.for
%
%   Usage: [Biodis_Data, Plant_Arch] = biodis_v301_for_app(Site_Data, Plant_Arch, PlantTypeDropDown)
%
%       Required Inputs:    Site_Data - table of site data information read from excel file or manually entered via GUI
%                           Plant_Arch - table of plant architecture parameters read from excel file
%                           PlantTypeDropDown - plant type designation read from excel file or manually entered via GUI

%       Provided Outputs:   Biodis_Data - vertical distribution of plant biomass through the water column
%                           Plant_Arch - updated plant architecture values based on new biomass distributions
%
%                     
% *************************************************************************************************************************

COPYRITE = 'Copyright (c) 2024 by RC Zimmerman';
VER= 'Version 3';

%% Determine the thickness of the canopy layers

Biodis_Data.nlayer = 1000;                                  % # of layers fixed at 1000
layer = transpose(1:1:Biodis_Data.nlayer);                  % create a numeric array for each layer
Biodis_Data.lathik = Site_Data.depth/Biodis_Data.nlayer;    % thickness of each layer (in meters) depends on water depth
Biodis_Data.ht = layer .* Biodis_Data.lathik;               % define height above the seafloor

%% Initialize Canopy Biomass Arrays and accumulators
Biodis_Data.rel_biomas = zeros(Biodis_Data.nlayer,1);               % relative biomass abundance in each layer (0 to 1)
Biodis_Data.extended_rel_biomas = zeros(Biodis_Data.nlayer,1);      % relative biomass fully extended
Biodis_Data.abs_biomas = zeros(Biodis_Data.nlayer,1);               % apparantly not used
Biodis_Data.lai = zeros(Biodis_Data.nlayer,1);                      % Leaf area (m2/m2) in each layer
Biodis_Data.lap = zeros(Biodis_Data.nlayer,1);                      % horizontally projected leaf area (m2/m2) in each layer
Biodis_Data.betan = zeros(Biodis_Data.nlayer,1);                    % apparently not used
Biodis_Data.beta = zeros(Biodis_Data.nlayer,1);                     % leaf bending angle, in radians
Biodis_Data.sumlai = 0;                                             % sum LAI in all layers
Biodis_Data.sumbio = 0;                                             % sum relative biomass in all layers (should = 1)
Biodis_Data.sumlap = 0;                                             % sum horizontally projected leaf area in all layers (LAP <= LAI)
Biodis_Data.canopy_ht = 0;                                          % final realized canopy height affected by water depth and flow
Biodis_Data.laitot = 0;                                             % total LAI, may be the same as Biodis_Data.sumlai

%% Identify the plant species for biomass calculations to follow  
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
    row = 6;        %Ruppia
 end

%% if plants exist
if Plant_Arch.density(row,col) > 0
% Calculate the canopy bending angle
    
    if (Site_Data.curvel == 0) && (Plant_Arch.beta0(row,col) == 0)                              % if both current velocity and beta0 are 0
        Biodis_Data.beta(1:Biodis_Data.nlayer,1) = deg2rad(1);                                  % set bending angle to 1 deg from the vertical
        Plant_Arch.beta0(row,col) = 1;
    elseif (Site_Data.curvel > 0)                                                               % otherwise, if current velocity > 0 
        Biodis_Data.beta(1:Biodis_Data.nlayer,1) = pi/2*((1-exp(-0.036* Site_Data.curvel)));    % set beta (in radians) according to current velocity
        Plant_Arch.beta0(row,col) = rad2deg(Biodis_Data.beta(1,1));                             % and bending angle in degrees for display - only 1 value, not entire array
    elseif (Site_Data.curvel == 0) && (Plant_Arch.beta0(row,col) > 0)                           % otherwise if beta0 > 0
        Biodis_Data.beta(1:Biodis_Data.nlayer,1) = deg2rad(Plant_Arch.beta0(row,col));          % just use beta0 from the input file and ignore current velocity
    end

% Calculate the shape, biomass fraction at seabed and inflection point from canopy height
% Sigmoid model parameters for strappy seagrass (Zostera, Thalassia, % Posidonia etc.)
        Plant_Arch.areash(row,col) = 0.0063 .* Plant_Arch.maxht(row,col) + 0.019 .* Plant_Arch.maxht(row,col) .^2;  % Leaf Area per shoot
        Plant_Arch.asymp(row,col) = (2.51 .* Plant_Arch.maxht(row,col) .^ (-0.79124)/100);                          % asymptotic fraction of biomass at the base of the canopy, divide by 100 converts percent to fraction and multiply by lathik to adjust for layer thickness
        Plant_Arch.shape(row,col) = 4.746816;                                                                       % Not correlated to height. A constant defined by average of all populations in Zimmerman 2003 paper
        Plant_Arch.inflec(row,col) = 0.588 .* (1 - exp(-1.12 .* Plant_Arch.maxht(row,col)));
    %end
    
Biodis_Data.laitot = Plant_Arch.areash(row,col) .* Plant_Arch.density(row,col);           % Total LAI, m2 leaf/m2 seafloor
Biodis_Data.leaf_biomass_tot = Plant_Arch.areash(row,col) .* 500;                         %convert shoot leaf area to g DW/shoot
Biodis_Data.flow_reduced_ht = Plant_Arch.maxht(row,col) .* cos(Biodis_Data.beta(1,1));    % Adjust realized canopy height to account for leaf bending; maxht (meters) provided by user

%% IF THE PLANT CANOPY IS FULLY SUBMERGED (i.e., water column is taller than canopy)

% compute the fully extended biomass distribution using new vectorized code
    if Site_Data.depth >= Biodis_Data.flow_reduced_ht                                                    % if canopy is fully submerged; 
        Plant_Arch.canopy_ht(row,col) = Biodis_Data.flow_reduced_ht;                                     %   canopy height (m) is controlled by flow
        Biodis_Data.n_extended_layers = floor(Plant_Arch.maxht(row,col)/Site_Data.depth * 1000);         %   number of canopy layers based on canopy height and water depth
        if mod(Biodis_Data.n_extended_layers,1) ~= 0                                                     %   round up fractional canopy layers, if they exist
            Biodis_Data.n_extended_layers = Biodis_Data.n_extended_layers + 1;  
        end
        Biodis_Data.extended_canopy_layer = transpose(1:1:Biodis_Data.n_extended_layers);                % create a numeric array for each canopy layer (<= number of depth layers)
        Biodis_Data.extended_canopy_ht = Biodis_Data.extended_canopy_layer .* Biodis_Data.lathik;        % create a numeric array defining height of each layer above the seafloor, in meters

        if Plant_Arch.cantype(row,col) == 1                                                              % if canopy structure is sigmoidal
            Biodis_Data.extended_rel_biomas(1:Biodis_Data.n_extended_layers,1) = (Plant_Arch.asymp(row,col)./(1.0 + (Biodis_Data.extended_canopy_ht./Plant_Arch.inflec(row,col)).^Plant_Arch.shape(row,col))) .* (Biodis_Data.lathik/0.01);  %sigmoid biomass function for flat, strap-like leaves; (lthik/0.01) corrects for differences in layer thickness caused by water depth
            Biodis_Data.lai = Biodis_Data.extended_rel_biomas * Biodis_Data.laitot;             % LAI in each layer, m2 leaf/m2 seafloor in each layer
            Biodis_Data.lap(1:Biodis_Data.n_extended_layers,1) = Biodis_Data.lai(1:Biodis_Data.n_extended_layers,1) .* sin(Biodis_Data.beta(1:Biodis_Data.n_extended_layers,1));    % Horizontally projected LAI in each layer, m2 leaf/m2 seafloor
        else 
            Biodis_Data.extended_rel_biomas(1:Biodis_Data.n_extended_layers,1) = -1.0608E-07*Plant_Arch.canopy_ht(row,col).^2 +8.6671E-05*Plant_Arch.canopy_ht(row,col)+4.3512E-04;  %ellipsoidall biomass function from Mandy Stoughton's thesis
            Biodis_Data.lai = Biodis_Data.rel_biomas .* Biodis_Data.laitot*0.5; %MANDY's modified for Syringodium (*.5,sees half the light due to the cylindrical shaped leaves of Syringodium)
        end

% adjust for canopy distribution integration error
        Biodis_Data.sum_extended_biomas = sum(Biodis_Data.extended_rel_biomas);
        Biodis_Data.extended_rel_biomas = Biodis_Data.extended_rel_biomas ./ Biodis_Data.sum_extended_biomas;    % correct for integration error

% adjust for bending angle
        Biodis_Data.n_flow_reduced_canopy_layers = Biodis_Data.flow_reduced_ht/Site_Data.depth * 1000;
        if mod(Biodis_Data.n_flow_reduced_canopy_layers,1) ~= 0        
            Biodis_Data.n_flow_reduced_canopy_layers = Biodis_Data.n_flow_reduced_canopy_layers - mod(Biodis_Data.n_flow_reduced_canopy_layers,1) + 1; % round up fractional canopy layers 
        end
        Biodis_Data.canopy_ht = Biodis_Data.n_flow_reduced_canopy_layers * Biodis_Data.lathik;
        Biodis_Data.flow_reduced_canopy_layer = transpose(1:1:Biodis_Data.n_flow_reduced_canopy_layers);               % create a numeric array for each canopy layer (<= number of depth layers)
        Biodis_Data.rel_biomas(1:Biodis_Data.n_flow_reduced_canopy_layers,1) = Biodis_Data.extended_rel_biomas(1:Biodis_Data.n_flow_reduced_canopy_layers,1) ./ cos(Biodis_Data.beta(1,1));
        Biodis_Data.rel_biomas(Biodis_Data.n_flow_reduced_canopy_layers:1000,1) = 0;
        Biodis_Data.sum_flow_reduced_biomas = sum(Biodis_Data.rel_biomas);
        Biodis_Data.rel_biomas = Biodis_Data.rel_biomas ./ Biodis_Data.sum_flow_reduced_biomas;
        Biodis_Data.lai = Biodis_Data.rel_biomas .* Biodis_Data.laitot;     % LAI in each layer, m2 leaf/m2 seafloor
        Biodis_Data.lap = Biodis_Data.lai .* sin(Biodis_Data.beta);         % horizontal projection of lai, m2 leaf/m2 seafloor       
        Biodis_Data.sum_rel_bio = sum(Biodis_Data.rel_biomas);
        Biodis_Data.sumlai = sum(Biodis_Data.lai);                          % Total LAI, m2 leaf/m2 seafloor
        Biodis_Data.sumlap = sum(Biodis_Data.lap);                          % Total horizontally projected LAI, m2 leaf/m2 seafloor
    end


%% IF CANOPY IS TALLER THAN WATER COLUMN, COMPUTE THE EXTENDED AND FLOW-COMPRESSED CANOPY DISTRIBUTION AS BEFORE, THEN
%   COMPRESS ALL REMAINING LAYERS INTO THE UPPER 10% OF THE CANOPY AND ASSUME THOSE LEAVES ARE HORIZONTAL, i.e., LAI=LAP
%   THERE WILL BE SOME BIOMASS IN ALL 1000 SIMULATED LAYERS OF THE MODEL

% compute the fully extended biomass distribution using new vectorized code
if Site_Data.depth < Biodis_Data.flow_reduced_ht                                                    % if plants are taller than water column
        Plant_Arch.canopy_ht(row,col) = Site_Data.depth;
        Biodis_Data.n_extended_layers = floor(Plant_Arch.maxht(row,col)/Biodis_Data.lathik);       % get number of canopy layers if fully extended for this depth and layer thickness
        if mod(Plant_Arch.maxht(row,col)/(Biodis_Data.lathik*1000),1) > 0                          % check for fractional layer left over
            Biodis_Data.n_extended_layers = Biodis_Data.n_extended_layers + 1;                     % if so, incraesed n_extended_layers by 1
        end
        Biodis_Data.extended_canopy_layer = transpose(1:1:Biodis_Data.n_extended_layers);   % create an numeric array for each layer of the extended canopy
        Biodis_Data.extended_canopy_ht = Biodis_Data.extended_canopy_layer .* (Biodis_Data.lathik);     % calculate the height of each extended canopy layer
    
        if Plant_Arch.cantype(row,col) == 1
            Biodis_Data.extended_rel_biomas(1:Biodis_Data.n_extended_layers,1) = (Plant_Arch.asymp(row,col)./(1.0 + (Biodis_Data.extended_canopy_ht./Plant_Arch.inflec(row,col)).^Plant_Arch.shape(row,col))) .* (Biodis_Data.lathik/0.01); %sigmoid biomass function for flat, strap-like leaves
            Biodis_Data.lai = Biodis_Data.extended_rel_biomas .* Biodis_Data.laitot;

        else
            Biodis_Data.extended_rel_biomas(1:Biodis_Data.n_extended_layers,1) = -1.0608E-07*Plant_Arch.canopy_ht(row,col).^2 +8.6671E-05*Plant_Arch.canopy_ht(row,col)+4.3512E-04;  %ellipsoidall biomass function from Mandy Stoughton's thesis
            Biodis_Data.lai = Biodis_Data.rel_biomas .* Biodis_Data.laitot*0.5; %MANDY's modified for Syringodium (*.5,sees half the light due to the cylindrical shaped leaves of Syringodium)
            lap(1,1:Biodis_Data.n_extended_layers) = lai(1,1:Biodis_Datan_extended_layers) .* sin(beta(1,1:Biodis_Datan_extended_layers)); %MANDY's modified for Syringodium (*.5,sees half the light due to the cylindrical shaped leaves of Syringodium)
        end

% adjust for canopy distribution integration error
        Biodis_Data.sum_extended_biomas = sum(Biodis_Data.extended_rel_biomas);
        Biodis_Data.extended_rel_biomas = Biodis_Data.extended_rel_biomas ./ Biodis_Data.sum_extended_biomas;    % correct for integration error
        %lai = extended_rel_biomas.*laitot;       % distribute lai according to rel_biomas

% adjust canopy height and biomass distribution for bending angle
        Biodis_Data.n_flow_reduced_canopy_layers = floor(Biodis_Data.flow_reduced_ht/Biodis_Data.lathik);       % round down to the nearest integer
        if mod(Biodis_Data.n_flow_reduced_canopy_layers,1) ~= 0        
            Biodis_Data.n_flow_reduced_canopy_layers = Biodis_Data.n_flow_reduced_canopy_layers + 1; % round up fractional canopy layers 
        end
        Biodis_Data.flow_reduced_canopy_layer = transpose(1:1:Biodis_Data.n_flow_reduced_canopy_layers);               % create a numeric array for each canopy layer (<= number of depth layers)
        Biodis_Data.flow_reduced_biomas(1:Biodis_Data.n_flow_reduced_canopy_layers,1) = Biodis_Data.extended_rel_biomas(1:Biodis_Data.n_flow_reduced_canopy_layers,1) ./ cos(Biodis_Data.beta(1,1));
        Biodis_Data.flow_reduced_biomas(Biodis_Data.n_flow_reduced_canopy_layers:1000,1) = 0;
        Biodis_Data.sum_flow_reduced_biomas = sum(Biodis_Data.flow_reduced_biomas);
        Biodis_Data.flow_reduced_biomas = Biodis_Data.flow_reduced_biomas ./ Biodis_Data.sum_flow_reduced_biomas;
        if Biodis_Data.n_flow_reduced_canopy_layers <= 1000
            Biodis_Data.rel_biomas = Biodis_Data.flow_reduced_biomas;
        end

% now aggregate layers extending above the water
        if Biodis_Data.n_flow_reduced_canopy_layers > 1000
            Biodis_Data.sum_floating_biomas = sum(Biodis_Data.flow_reduced_biomas(1001:Biodis_Data.n_flow_reduced_canopy_layers,1));
            Biodis_Data.floating_biomas_increment = Biodis_Data.sum_floating_biomas / 100;
            Biodis_Data.rel_biomas = Biodis_Data.flow_reduced_biomas(1:1000,1);
            Biodis_Data.rel_biomas(901:1000,1) = Biodis_Data.rel_biomas(901:1000,1) + Biodis_Data.floating_biomas_increment;
            Biodis_Data.n_flow_reduced_canopy_layers = 1000;
        end
        Biodis_Data.sumbio = sum(Biodis_Data.rel_biomas);
        Biodis_Data.rel_biomas = Biodis_Data.rel_biomas ./ Biodis_Data.sumbio;
        Biodis_Data.lai = Biodis_Data.rel_biomas .* Biodis_Data.laitot;
        Biodis_Data.lap(1:900,1) = Biodis_Data.lai(1:900,1) .* sin(Biodis_Data.beta(1,1));       % horizontal projection of lai
        Biodis_Data.lap(901:1000,1) = Biodis_Data.lai(901:1000,1); 
        Biodis_Data.sumlai = sum(Biodis_Data.lai);
        Biodis_Data.sumlap = sum(Biodis_Data.lap);       
    end

end

% End of Biodis function

