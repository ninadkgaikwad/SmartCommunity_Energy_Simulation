function [PVEnergy_Available] = HEMS_PVEnergy_Available_Generator(PVEnergy_Generator_Input)

% Author: Ninad Kiran Gaikwad
% Date: Feb/12/2021
% Description: HEMS_PVEnergyAvailable_Generator - PV Energy Model

%% HEMS_PVEnergyAvailable_Generator - PV Energy Model

%% Getting desired Data from Input - Structs

PV_TotlaModules_Num=PVEnergy_Generator_Input.PV_TotlaModules_Num;
PV_RatedPower=PVEnergy_Generator_Input.PV_RatedPower;
PV_TempCoeff=PVEnergy_Generator_Input.PV_TempCoeff;

T_am=PVEnergy_Generator_Input.T_am;
GHI=PVEnergy_Generator_Input.GHI;
Ws=PVEnergy_Generator_Input.Ws;

Uo=PVEnergy_Generator_Input.Uo;
U1=PVEnergy_Generator_Input.U1;

Temp_Std=PVEnergy_Generator_Input.Temp_Std;
GHI_Std=PVEnergy_Generator_Input.GHI_Std;

Simulation_StepSize=PVEnergy_Generator_Input.Simulation_StepSize;

%% Computing Total Available PV Energy

PV_Total_Power = PV_TotlaModules_Num * PV_RatedPower;

for ii=1:length(GHI)

    Tm=T_am(ii)+((GHI(ii))/(Uo+(U1*Ws(ii)))); % Faiman's Module Temperature Model

    PVEnergy_Available(ii)=PV_Total_Power*(1+((PV_TempCoeff/100)*(Tm-Temp_Std)))*(GHI(ii)/GHI_Std)*Simulation_StepSize/1000; % Power generated in one Module    

end

end