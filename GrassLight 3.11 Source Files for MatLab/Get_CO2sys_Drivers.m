function [Param1 Param2] = Get_CO2sys_drivers(WQ_Params)
%   -----------------------------------------------------------------------------------------------------------------------------------------
%  Function to assign values for Param1 and Param2 from the data file 
%         or from the CO2SYSInputParametersSDropDown callback function in the GL GUI to drive the CO2sys calculation
%
% 
%       Version         Created by          Date            Comments
%       --------------------------------------------------------------------------------------------------------
%          1.0              RCZ         20 March 2023     Intended to simplify callback code lines in GL3 App 
%
%   Usage: [Param1 Param2] = Get_CO2sys_drivers(WQ_Params)
%       Param1: Value of first parameter required to drive CO2SYS calculation of the other parameters
%       Param2: Value of second parameter required to drive CO2Sys calculation of the other parameters
%
%   Call this function just before calling CO2SYS_for_app
%
%   Data required:
%        WQ_Params: structured array of water quality data read from Water Quality Paramaters worksheet in GL3_Setup_Default.xlsx
%
%    Chemical definitions for Parameter Type values are listed in the header info in function CO2SYS_for_ m
%   ------------------------------------------------------------------------------------------------------------------------------------------

if WQ_Params.PAR1TYPE == 1 &  WQ_Params.PAR2TYPE == 2 |  WQ_Params.PAR1TYPE == 2 &  WQ_Params.PAR2TYPE == 1
    Param1 =  WQ_Params.Alk;
    Param2 =  WQ_Params.TCO2;
    return
elseif  WQ_Params.PAR1TYPE == 1 &  WQ_Params.PAR2TYPE == 3 |  WQ_Params.PAR1TYPE == 3 &  WQ_Params.PAR2TYPE == 1
    Param1 =  WQ_Params.Alk;
    Param2 =  WQ_Params.pH;
    return
elseif  WQ_Params.PAR1TYPE == 1 &  WQ_Params.PAR2TYPE == 4 |  WQ_Params.PAR1TYPE == 4 &  WQ_Params.PAR2TYPE == 1
    Param1 =  WQ_Params.Alk;
    Param2 =  WQ_Params.pCO2_Water;
    return
elseif  WQ_Params.PAR1TYPE == 2 &  WQ_Params.PAR2TYPE == 3 |  WQ_Params.PAR1TYPE == 3 &  WQ_Params.PAR2TYPE == 2
    Param1 =  WQ_Params.TCO2;
    Param2 =  WQ_Params.pH;
    return
elseif  WQ_Params.PAR1TYPE == 2 &  WQ_Params.PAR2TYPE == 4 |  WQ_Params.PAR1TYPE == 4 &  WQ_Params.PAR2TYPE == 2
     Param1 =  WQ_Params.TCO2;
     Param2 =  WQ_Params.pCO2_Water;
     return
else
      Param1 =  WQ_Params.pH;
      Param2 =  WQ_Params.pCO2_Water;
end