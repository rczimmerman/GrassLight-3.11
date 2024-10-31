function [IOP_Spectra_Header, IOP_Spectra, Site_Data_Header, Site_Data, Atm_Const_Header, Atm_Const, Atm_IOPs_Header, Atm_IOPs, IOP_Constants_Header, IOP_Constants, WQ_Params_Header, WQ_Params, Plant_Arch_Header, Plant_Arch, Photo_Params_Header, Photo_Params, Sed_Ref_Header, Sed_Ref, LOPs_Header, LOPs] = read_GL3_input_data_for_app_302(inpath, infile);


COPYRITE = 'Copyright (c) 2024 by RC Zimmerman';
VER= 'Version 3';
%% *******************************Read Input Data from GL3_Setup_Default******************************************************

%               inpath = 'C:\RCZ\Projects\Canopy Model calcuations & notes\Matlab GrassLight\GL 3.0\GL3 New App\    % debugging lines for running as a stand-alone script            
%               infile = 'GL3_Setup_Default.xlsx';                                                                  % debugging lines for running as a stand-alone script 

setup_file = [inpath infile];       % Name of the input file comes from the GL GUI

%% Get Site Data from GL3_Setup_Defalut.xlsx


% Read Header Data First
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify the Header Range
opts.DataRange = "A1:G9";
opts.Sheet = "Site Data";

% Import the Site Data Header Info
Site_Data_Header = readtable(setup_file, opts, "UseExcel", false);

% Now get the numeric data required for GL
opts = spreadsheetImportOptions("NumVariables", 8);
opts.DataRange = "A12:H12";

% Specify column names and types
opts.VariableNames = ["XLAT", "XLON", "JD", "HR", "depth", "curvel","sed_type",'lop'];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double","string","string"];

% Import the data
Site_Data = readtable(setup_file, opts, "UseExcel", false);

clear opts;

%% Get IOP spectra from GL3_Setup_Defalut.xlsx

opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "IOP Spectra";
opts.DataRange = "A1:E9";

% Import the IOP Spectra Header Data
IOP_Spectra_Header = readtable(setup_file, opts, "UseExcel", false);

% Now get the IOP Data required for GL calculations
opts = spreadsheetImportOptions("NumVariables", 4);

% Specify sheet and range
opts.Sheet = "IOP Spectra";
opts.DataRange = "A13:D313";

% Specify column names and types
opts.VariableNames = ["WL", "AW", "BW", "ASTCHL"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Import the data
IOP_Spectra = readtable(setup_file, opts, "UseExcel", false);

clear opts;

%% Get the IOP Constants needed to calculate attenuation

opts = spreadsheetImportOptions("NumVariables", 10);

% Specify sheet and range
opts.Sheet = "IOP Constants";
opts.DataRange = "A1:J11";

% Import the header data
IOP_Constants_Header = readtable(setup_file, opts, "UseExcel", false);

% Now get the numeric data needed for GL
opts = spreadsheetImportOptions("NumVariables", 8);

% Specify sheet and range
opts.Sheet = "IOP Constants";
opts.DataRange = "A13:H13";

% Specify column names and types
opts.VariableNames = ["SG", "SIGATRB", "STRB", "BLTRB", "SIGBTRB", "ETA", "BB2B", "APHST675"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
IOP_Constants = readtable(setup_file, opts, "UseExcel", false);

clear opts;

%% Get the Water Quality Parameters needed to calculate attenuation and plant metabolism

% First read the Water Quality Header Data
opts = spreadsheetImportOptions("NumVariables", 6);

% Specify sheet and range
opts.Sheet = "Water Quality Parameters";
opts.DataRange = "A1:F27";

% Import the header data
WQ_Params_Header = readtable(setup_file, opts, "UseExcel", false);

% Now get the numeric water quality parameters
opts = spreadsheetImportOptions("NumVariables", 27);

% Specify sheet and range
opts.Sheet = "Water Quality Parameters";
opts.DataRange = "A31:AA31";

% Specify column names and types
opts.VariableNames = ["AG_440", "CHL_A", "TSM", "Temp", "Sal", "Alk", "O2_uM", "O2_percent", "pH", "scale", "SiO4", "PO4", "NH4", "NO3", "H2S", "pCO2_Air", "pCO2_Water", "CO2", "HCO3", "CO3", "TCO2", "K1K2", "KHSO4", "KHF", "Borate", "PAR1TYPE", "PAR2TYPE"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
WQ_Params = readtable(setup_file, opts, "UseExcel", false);

% % Convert pH scale from string to numeric value required by CO2SYS
if WQ_Params.scale == "total" || WQ_Params.scale =="Total"
    WQ_Params.pH_Scale = 1;
elseif WQ_Params.scale == "seawater"|| WQ_Params.scale =="Seawater"
    WQ_Params.pH_Scale = 2;
elseif WQ_Params.scale == "free" || WQ_Params.scale =="Free"
    WQ_Params.pH_Scale = 3;
elseif WQ_Params.scale == "nbs" || WQ_Params.scale =="NBS"
    WQ_Params.pH_Scale = 4;
end

clear opts;

%% Get Atmospheric IOP Arrays for H0, ozone, water and oxygen

% Atmospheric IOP Header Data First
opts = spreadsheetImportOptions("NumVariables", 3);

% Specify sheet and range
opts.Sheet = "Atmospheric Optical Properties";
opts.DataRange = "A1:C7";

% Import the header data
Atm_IOPs_Header = readtable(setup_file, opts, "UseExcel", false);

% Now get the optical property spectra
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Atmospheric Optical Properties";
opts.DataRange = "A12:E312";

% Specify column names and types
opts.VariableNames = ["LAMBDA1", "H0", "AOZ", "AW_INC", "AOX"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
Atm_IOPs = readtable(setup_file, opts, "UseExcel", false);

clear opts;

%% Get Atmospheric Constants from GL3_Setup_Defalut.xlsx

% First read the Header Data

opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Atmospheric Constants";
opts.DataRange = "A1:B10";

% Import the header data
Atm_Const_Header = readtable(setup_file, opts, "UseExcel", false);

% Now get the numeric data
opts = spreadsheetImportOptions("NumVariables", 8);

% Specify sheet and range
opts.Sheet = "Atmospheric Constants";
opts.DataRange = "A13:H13";

% Specify column names and types
opts.VariableNames = ["AM", "WM", "W", "RH", "PRESS", "WV", "HA", "V"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
Atm_Const = readtable(setup_file, opts, "UseExcel", false);

clear opts;


%% Get Canopy Architecture Parameters 

% Get the header data first
opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Plant Architecture";
opts.DataRange = "A1:B15";

% Import the data
Plant_Arch_Header = readtable(setup_file, opts, "UseExcel", false);

% Now get the numeric data
opts = spreadsheetImportOptions("NumVariables", 13);

% Specify sheet and range
opts.Sheet = "Plant Architecture";
opts.DataRange = "A18:M23";

% Specify column names and types
opts.VariableNames = ["lop", "cantype", "maxht", "shoot_frac", "rhiz_frac", "root_frac", "density", "areash", "asymp", "inflec", "shape", "beta0", "epi"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
Plant_Arch = readtable(setup_file, opts, "UseExcel", false);

% 
if Plant_Arch.density <0                    % if shoot denisty is entered as 0, just make it really small
    Plant_Arch.density  = 0;             % to prevent GL from blowing up 
end

%  Clear temporary variables
clear opts


%% Get Metabolic Parameters 
% Header data first
opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Metabolic Parameters";
opts.DataRange = "A1:B16";

% Import the header data
Photo_Params_Header = readtable(setup_file, opts, "UseExcel", false);

% Now get the numeric data
opts = spreadsheetImportOptions("NumVariables", 12);

% Specify sheet and range
opts.Sheet = "Metabolic Parameters";
opts.DataRange = "A21:L26";

% Specify column names and types
opts.VariableNames = ["Species", "Pm_CO2", "Ks_CO2", "Upf_CO2", "Up0_CO2", "beta_CO2", "Pm_HCO3", "Ks_HCO3", "Upf_HCO3", "Up0_HCO3", "beta_HCO3", "Leaf_R_21_deg"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
Photo_Params = readtable(setup_file, opts, "UseExcel", false);


% Clear temporary variables
clear opts

%% Get Sediment Reflectance Spectrum

% First read the header data
opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Sediment Reflectances";
opts.DataRange = "A1:B8";

% Import the data
Sed_Ref_Header = readtable(setup_file, opts, "UseExcel", false);

opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Sediment Reflectances";
opts.DataRange = "A11:E311";

% Specify column names and types
opts.VariableNames = ["WL", "Bright_Carbonate", "Dark_Carbonate", "Silica_Sand", "Silica_Mud"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
Sed_Ref = readtable(setup_file, opts, "UseExcel", false);

% Clear temporary variables
clear opts

%% Get Leaf Optical Properties

% First get the header data
opts = spreadsheetImportOptions("NumVariables", 1);

% Specify sheet and range
opts.Sheet = "Leaf Optical Properties";
opts.DataRange = "A1:A7";

% Import the data
LOPs_Header = readtable(setup_file, opts, "UseExcel", false);

% Now get the leaf optical property data
opts = spreadsheetImportOptions("NumVariables", 13);

% Specify sheet and range
opts.Sheet = "Leaf Optical Properties";
opts.DataRange = "A11:M311";

% Specify column names and types
opts.VariableNames = ["WL", "a_zostera", "r_zostera", "a_thalassia", "r_thalassia", "a_syringodium", "r_syringodium", "a_phyllospadix", "r_phyllospadix", "a_posidonia", "r_posidonia", "a_ruppia", "r_ruppia"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double",];

% Specify variable properties
opts = setvaropts(opts, ["a_phyllospadix", "r_phyllospadix", "a_posidonia", "r_posidonia", "a_ruppia", "r_ruppia"], "EmptyFieldRule", "auto");

% Import the data
LOPs = readtable(setup_file, opts, "UseExcel", false);

% Clear temporary variables
clear opts;



% ******************************************End Read Input Data from Excel % File ***************************************************