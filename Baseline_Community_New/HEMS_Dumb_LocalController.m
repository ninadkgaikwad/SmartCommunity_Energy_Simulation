function [U_k] = HEMS_Dumb_LocalController(X_k_Plant,W_k_Plant,HEMSPlant_Params,Community_Params,Simulation_Params)

% Author: Ninad Kiran Gaikwad
% Date: Feb/12/2021
% Description: HEMS_Dumb_LocalController - Dumb Local Logic

%% HEMS_Dumb_LocalController - Dumb Local Logic

%% Getting desired Data from Input - Structs

% From X_k_Plant
T_house_now = X_k_Plant(1,7,:);
E_bat_now = X_k_Plant(1,4,:);

% From W_k_Plant
Weather_k_Plant=W_k_Plant.Weather_k_Plant;
LoadData_k_Plant=W_k_Plant.LoadData_k_Plant;

GHI=Weather_k_Plant.GHI;
T_am=Weather_k_Plant.T_am;
Ws=Weather_k_Plant.Ws;

E_Load_Desired=LoadData_k_Plant.E_Load_Desired;
E_LoadData=LoadData_k_Plant.E_LoadData;

% From HEMSPlant_Params
T_AC_Base=HEMSPlant_Params.T_AC_Base;
T_AC_DeadBand=HEMSPlant_Params.T_AC_DeadBand;

PV_TotlaModules_Num=HEMSPlant_Params.PV_TotlaModules_Num;
PV_RatedPower=HEMSPlant_Params.PV_RatedPower;
PV_TempCoeff=HEMSPlant_Params.PV_TempCoeff;
GHI_Std=HEMSPlant_Params.GHI_Std;
Temp_Std=HEMSPlant_Params.Temp_Std;
Eff_Inv=HEMSPlant_Params.Eff_Inv;
Uo=HEMSPlant_Params.Uo;
U1=HEMSPlant_Params.U1;

% From Community_Params
N_House=Community_Params.N_House;
N_PV_Bat=Community_Params.N_PV_Bat;
N_Bat=Community_Params.N_Bat;
N_PV=Community_Params.N_PV;
N_None=Community_Params.N_None;   

% From Simulation_Params
Simulation_StepSize=Simulation_Params.Simulation_StepSize;

%% Initializing U_k

U_k = zeros(1,11,N_House);

%% AC Thermostat Controller

% For Each House in this Time Step Computing PV Energy Available and Total Load
for jj=1:N_House

    % Computing U_ac 
    Baseline_ACController_Input.T_AC_Base=T_AC_Base;
    Baseline_ACController_Input.T_AC_DeadBand=T_AC_DeadBand;
    Baseline_ACController_Input.u_k_hvac_prev=X_k_Plant(1,30,jj);
    Baseline_ACController_Input.T_house_now=X_k_Plant(1,7,jj);

    [u_k_ac] = HEMS_AC_Controller(Baseline_ACController_Input);

    U_k(1,3,jj)=u_k_ac; % Updating U_ac

end  

%% Battery Logic Controller

% Creating Input for PV Energy Available Generator
PVEnergy_Generator_Input.PV_TotlaModules_Num=PV_TotlaModules_Num;
PVEnergy_Generator_Input.PV_RatedPower=PV_RatedPower;
PVEnergy_Generator_Input.PV_TempCoeff=PV_TempCoeff;

PVEnergy_Generator_Input.T_am=T_am;
PVEnergy_Generator_Input.GHI=GHI;
PVEnergy_Generator_Input.Ws=Ws;

PVEnergy_Generator_Input.Uo=Uo;
PVEnergy_Generator_Input.U1=U1;

PVEnergy_Generator_Input.Temp_Std=Temp_Std;
PVEnergy_Generator_Input.GHI_Std=GHI_Std;

PVEnergy_Generator_Input.Simulation_StepSize=Simulation_StepSize;

% Computing Total Available PV Energy
[PVEnergy_Available] = HEMS_PVEnergy_Available_Generator(PVEnergy_Generator_Input);

Total_PVEnergy_Available = (N_PV_Bat+N_PV)*PVEnergy_Available;

% Computing Total Load Desired
E_Load_Total_Desired=sum(E_Load_Desired);

% Creating Input for Battery Logic Controller
Baseline_BatteryController_Input.PVEnergy_Available=Total_PVEnergy_Available;
Baseline_BatteryController_Input.E_Tot_Desired=E_Load_Total_Desired;

[c_k,d_k] = HEMS_Battery_Controller(Baseline_BatteryController_Input);

% For Each House with Battery Storage in this Time Step
for jj=1:(N_PV_Bat+N_Bat)
    
    U_k(1,1,jj)=c_k;
    U_k(1,2,jj)=d_k;
    
end

%% Secondary Load Controller

% Basic Computation
[R,C,D]=size(E_LoadData);

% For Each House in this Time Step
for jj=1:N_House
    
    % For Each Prioritized Load Level in this Time Step
    for kk=9+1:C
    
        if (E_LoadData(1,kk,jj)>0) % Prioritized Load is Demanded
            
            U_k(1,(kk-9+3),jj)=1;
            
        else % Prioritized Load is not Demanded
            
            U_k(1,(kk-9+3),jj)=0;
            
        end
    
    end
    
end

end



