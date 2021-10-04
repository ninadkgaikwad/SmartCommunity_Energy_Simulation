function [u_k_hvac_bc] = HEMS_AC_Controller(BaselineController_Input)

% Author: Ninad Kiran Gaikwad
% Date: Feb/12/2021
% Description: HEMS_Baseline_Controller - Baseline Controller Logic

%% HEMS_Baseline_Controller - Baseline Controller Logic

%% Getting desired Data from the BaselineController_Input - Struct

T_AC_Base=BaselineController_Input.T_AC_Base;
T_AC_DeadBand=BaselineController_Input.T_AC_DeadBand;
u_k_hvac_prev=BaselineController_Input.u_k_hvac_prev;
T_house_now=BaselineController_Input.T_house_now;

%% Rule Based Logic Controller - HVAC

if ((u_k_hvac_prev==0) && (T_house_now>=(T_AC_Base+T_AC_DeadBand)))

     
     u_k_hvac_bc=1;

elseif ((u_k_hvac_prev==1) && ((T_house_now<=(T_AC_Base-T_AC_DeadBand))))

    
    u_k_hvac_bc=0;
    
else
    
    u_k_hvac_bc=u_k_hvac_prev;

end 
    

end

