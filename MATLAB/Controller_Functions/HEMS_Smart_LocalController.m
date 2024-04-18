function [U_k] = HEMS_Smart_LocalController(X_k_Plant,W_k_Plant,HEMSPlant_Params,Community_Params,Simulation_Params)

% Author: Ninad Kiran Gaikwad
% Date: Mar/10/2021
% Description: HEMS_Smart_LocalController - Smart Local Controller Logic

%% HEMS_Smart_LocalController - Smart Local Controller Logic

%% Getting desired Data from Input - Structs

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
ACLoad_StartUp_Power=HEMSPlant_Params.ACLoad_StartUp_Power;

PV_TotlaModules_Num=HEMSPlant_Params.PV_TotlaModules_Num;
PV_RatedPower=HEMSPlant_Params.PV_RatedPower;
PV_TempCoeff=HEMSPlant_Params.PV_TempCoeff;
GHI_Std=HEMSPlant_Params.GHI_Std;
Temp_Std=HEMSPlant_Params.Temp_Std;
Eff_Inv=HEMSPlant_Params.Eff_Inv;
Uo=HEMSPlant_Params.Uo;
U1=HEMSPlant_Params.U1;

Eff_Charging_Battery=HEMSPlant_Params.Eff_Charging_Battery;
Eff_Discharging_Battery=HEMSPlant_Params.Eff_Discharging_Battery;

MaxRate_Charging=HEMSPlant_Params.MaxRate_Charging;
MaxRate_Discharging=HEMSPlant_Params.MaxRate_Discharging;
Battery_Energy_Max=HEMSPlant_Params.Battery_Energy_Max;
Battery_Energy_Min=HEMSPlant_Params.Battery_Energy_Min;
MaxRate_Discharging_StartUp=HEMSPlant_Params.MaxRate_Discharging_StartUp;

E_AC=HEMSPlant_Params.E_AC;

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

%% Getting Control Commands from the Local Dumb Controller

[U_k_DumbLocalController] = HEMS_Dumb_LocalController(X_k_Plant,W_k_Plant,HEMSPlant_Params,Community_Params,Simulation_Params);

%% Basic Computations

%------------------------- Load Energy at DC Side ------------------------%

E_Load_Desired_DC=E_Load_Desired/Eff_Inv;
E_LoadData_DC=[E_LoadData(1,1:9,:) E_LoadData(1,9+1:end,:)/Eff_Inv];

%------------------- Getting Total PV Energy Available -------------------%

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

%---------- Battery Energy Available for Discharging dispatch ------------%

E_Bat_Discharging_Dispatch=zeros(1,1,N_House);

for jj=1:N_PV_Bat+N_Bat % For Battery installed Houses
    
    E_Bat_Discharging_Dispatch(1,1,jj)=(U_k_DumbLocalController(1,2,jj))*...
        (min(MaxRate_Discharging*Simulation_StepSize,...
        (E_bat_now(1,1,jj)-Battery_Energy_Min)*Eff_Discharging_Battery));
    
end

%----------- Battery Energy Available for Charging dispatch --------------%

E_Bat_Charging_Dispatch=zeros(1,1,N_House);

for jj=1:N_PV_Bat+N_Bat % For Battery installed Houses
    
    E_Bat_Charging_Dispatch(1,1,jj)=(U_k_DumbLocalController(1,1,jj))*...
        (min(MaxRate_Charging*Simulation_StepSize,...
        (Battery_Energy_Max-(E_bat_now(1,1,jj)))/Eff_Charging_Battery));
    
end

%% Computing Important Indices for Gauging Status of ACs and Batteries

%----------------------- Important Battery Indices -----------------------%

% Battery Discharging Indices
BatteryDischarging_Indices=find(U_k_DumbLocalController(1,2,:)==1);
Battery_AboveMin_Indices=find(E_bat_now(1,1,:)>Battery_Energy_Min);
Battery_ActualDischarging_Indices=intersect(BatteryDischarging_Indices,Battery_AboveMin_Indices);

% Battery Charging Indices
BatteryCharging_Indices=find(U_k_DumbLocalController(1,1,:)==1);
Battery_BelowMax_Indices=find(E_bat_now(1,1,:)<Battery_Energy_Max);
Battery_ActualCharging_Indices=intersect(BatteryCharging_Indices,Battery_BelowMax_Indices);

%------------------------- Important AC Indices --------------------------%

% AC Turn-ON Conditions
AC_Indices=1:N_House;

% AC Turning On Indices
TurningOn_AC_Vector=U_k_DumbLocalController(1,3,:)-X_k_Plant(1,30,:);
TurningOn_AC_Indices=find(TurningOn_AC_Vector==1);

% AC Turned On Indices
TurnedOn_AC_Vector=U_k_DumbLocalController(1,3,:)+X_k_Plant(1,30,:);
TurnedOn_AC_Indices=find(TurnedOn_AC_Vector==2);

% AC Turned Off Indices
TurnedOff_AC_Vector=U_k_DumbLocalController(1,3,:)+X_k_Plant(1,30,:);
TurnedOff_AC_Indices=find(TurnedOff_AC_Vector==0);

% AC Turning Off Indices
TurningOff_AC_Vector=U_k_DumbLocalController(1,3,:)-X_k_Plant(1,30,:);
TurningOff_AC_Indices=find(TurningOff_AC_Vector==-1);

% AC Not Turning On Indices
Not_TurningOn_AC_Indices=AC_Indices(~ismember(AC_Indices,TurningOn_AC_Indices));

%% Computing Energy Mismatch and Number of different Status ACs

%---------------------- Computing Energy Mismatch ------------------------%

% Total energy Desired
TotalEnergy_Desired=((E_AC/Eff_Inv)*(sum(U_k_DumbLocalController(1,3,:))))+sum(E_Load_Desired_DC)+sum(E_Bat_Charging_Dispatch(1,1,Battery_ActualCharging_Indices));

% Energy Mismatch
AvailableEnergy=Total_PVEnergy_Available+sum(E_Bat_Discharging_Dispatch(1,1,Battery_ActualDischarging_Indices));
E_Mis=AvailableEnergy-TotalEnergy_Desired;

%--------------- Computing Number of different Status ACs ----------------%

Num_TurningOn_AC=length(TurningOn_AC_Indices);
Num_TurnedOn_AC=length(TurnedOn_AC_Indices);
Num_TurningOff_AC=length(TurningOff_AC_Indices);
Num_TurnedOff_AC=length(TurnedOff_AC_Indices);

%------------ Computing Total Startup Power Required for ACs -------------%

Battery_CapableofDischarge_Indices=intersect(find(U_k_DumbLocalController(1,2,:)),find(E_bat_now(1,1,:)>Battery_Energy_Min));
TotalPower_Available_Bat=length(Battery_CapableofDischarge_Indices)*(MaxRate_Discharging_StartUp);
TotalPower_Available_PV=(Total_PVEnergy_Available)/(Simulation_StepSize);
Total_StartUpPower_Available=TotalPower_Available_Bat+TotalPower_Available_PV;
Total_StartUpPower_AC=Num_TurningOn_AC*ACLoad_StartUp_Power;

%% Computing Control Commands

if (E_Mis>=0) % Loads can be serviced completely

    % Computing Battery Charging-Discharging Control Commands (Charging Batteries not shed)
    U_k(1,1:2,:)=U_k_DumbLocalController(1,1:2,:);
    
    % Computing Loads On-Off Control Commands (All Loads can be serviced)
    U_k(1,4:end,:)=U_k_DumbLocalController(1,4:end,:);

    % Computing Actual Load Serviced AC
    if(Total_StartUpPower_Available>=Total_StartUpPower_AC) % Enough Batteries/PV to Power the Turn On of ACs

        % Computing AC On-Off Control Commands (All Turning On ACs can be Turned on and all Turned On ACs can be serviced)
        U_k(1,3,:)=U_k_DumbLocalController(1,3,:);

    else

        Num_TurnOn_AC_Possible=min([Num_TurningOn_AC,floor(Total_StartUpPower_Available/ACLoad_StartUp_Power)]);

        if (isempty(TurningOn_AC_Indices)) % Check for using randsample()
            TurnOn_AC_Possible_Indices=[];
        else
            TurnOn_AC_Possible_Indices=randsample(TurningOn_AC_Indices,Num_TurnOn_AC_Possible);                       
        end
        TurnOn_AC_NotPossible_Indices=TurningOn_AC_Indices(~ismember(TurningOn_AC_Indices,TurnOn_AC_Possible_Indices));

        % Turn On AC
        U_k(1,3,TurnOn_AC_Possible_Indices)=U_k_DumbLocalController(1,3,TurnOn_AC_Possible_Indices);

        % Didnt Turn On AC
        U_k(1,3,TurnOn_AC_NotPossible_Indices)=zeros(1,1,length(TurnOn_AC_NotPossible_Indices));

        % Not Turning On AC (Rest of the ACs)
        U_k(1,3,Not_TurningOn_AC_Indices)=U_k_DumbLocalController(1,3,Not_TurningOn_AC_Indices);

    end

elseif (E_Mis<0) % Loads cannot be serviced completely

    if (abs(E_Mis)<=sum(E_Bat_Charging_Dispatch(1,1,BatteryCharging_Indices))) % E_Mis can be shedded through idleing Charging Batteries

        % Computing Battery Discharging Control Commands (Charging Batteries are shed)
        U_k(1,2,:)=U_k_DumbLocalController(1,2,:);

        % Computing Loads On-Off Control Commands (All Loads can be serviced)
        U_k(1,4:end,:)=U_k_DumbLocalController(1,4:end,:);

        % Computing Charging Batteries to be made Idle
        [~,Sorted_BatteryChargingLevel_ChargingBat_Indices]=sort(E_Bat_Charging_Dispatch(1,1,BatteryCharging_Indices));
        [Descending_BatteryChargingLevel_ChargingBat_Indices]=flip(Sorted_BatteryChargingLevel_ChargingBat_Indices);

        ChargingBat_Energy_Sum=0; % Initalization
        IdleBat_Counter=0; % Initialization
        for jj=Descending_BatteryChargingLevel_ChargingBat_Indices % For all charging batteries

            if (ChargingBat_Energy_Sum>=abs(E_Mis))

                break;

            else

                IdleBat_Counter=IdleBat_Counter+1; % Incrementing IdleBat_Counter

                ChargingBat_Energy_Sum=ChargingBat_Energy_Sum+E_Bat_Charging_Dispatch(1,1,jj);

                IdleBattery_Indices(IdleBat_Counter)=jj;                     

            end

        end
        
        BatteryCharging_Possible_Indices=setdiff(BatteryCharging_Indices,IdleBattery_Indices);
        
        % Computing Battery Charging Control Commands (Charging Batteries are shed)
        U_k(1,1,BatteryCharging_Possible_Indices)=U_k_DumbLocalController(1,1,BatteryCharging_Possible_Indices);         
        
        % Computing Battery Charging Control Commands (Charging Batteries are shed)
        U_k(1,1,IdleBattery_Indices)=0;        

        % Computing Actual Load Serviced AC
        if(Total_StartUpPower_Available>=Total_StartUpPower_AC) % Enough Batteries/PV to Power the Turn On of ACs

            % Computing AC On-Off Control Commands (All Turning On ACs can be Turned on and all Turned On ACs can be serviced)
            U_k(1,3,:)=U_k_DumbLocalController(1,3,:);

        else

            Num_TurnOn_AC_Possible=min([Num_TurningOn_AC,floor(Total_StartUpPower_Available/ACLoad_StartUp_Power)]);

            if (isempty(TurningOn_AC_Indices)) % Check for using randsample()
                TurnOn_AC_Possible_Indices=[];
            else
                TurnOn_AC_Possible_Indices=randsample(TurningOn_AC_Indices,Num_TurnOn_AC_Possible);                      
            end                
            TurnOn_AC_NotPossible_Indices=TurningOn_AC_Indices(~ismember(TurningOn_AC_Indices,TurnOn_AC_Possible_Indices));

            % Turn On AC
            U_k(1,3,TurnOn_AC_Possible_Indices)=U_k_DumbLocalController(1,3,TurnOn_AC_Possible_Indices);

            % Didnt Turn On AC
            U_k(1,3,TurnOn_AC_NotPossible_Indices)=zeros(1,1,length(TurnOn_AC_NotPossible_Indices));

            % Not Turning On AC (Rest of the ACs)
            U_k(1,3,Not_TurningOn_AC_Indices)=U_k_DumbLocalController(1,3,Not_TurningOn_AC_Indices);
        
        end 
        

elseif (abs(E_Mis)<=sum(E_Bat_Charging_Dispatch(1,1,BatteryCharging_Indices))+((E_AC/Eff_Inv)*(sum(U_k_DumbLocalController(1,3,:))))) % E_Mis can be shedded through making idle Charging batteries turning off ACs

        % Computing Battery Discharging Control Commands (All Charging Batteries are shed)
        U_k(1,2,:)=U_k_DumbLocalController(1,2,:);
        
        % Computing Battery Charging Control Commands (All Charging Batteries are shed)
        U_k(1,1,:)=zeros(1,1,N_House);           

        % Computing Loads On-Off Control Commands (All Loads can be serviced)
        U_k(1,4:end,:)=U_k_DumbLocalController(1,4:end,:);
        
        % Computing Number of ACs to be turned Off
        Number_AC_Shedded=min([ceil((abs(E_Mis)-sum(E_Bat_Charging_Dispatch(1,1,BatteryCharging_Indices)))/E_AC),N_House]);

        if ((Num_TurningOn_AC-Number_AC_Shedded)>=0) % E_mis can be shedded through TurningOn ACs

            % Recomputing Total_StartUpPower_AC
            Total_StartUpPower_AC_1=(Num_TurningOn_AC-Number_AC_Shedded)*ACLoad_StartUp_Power;        

            % Computing Actual Load Serviced AC
            if(Total_StartUpPower_Available>=Total_StartUpPower_AC_1) % Enough Batteries/PV to Power the Turn On of ACs

                Num_TurnOn_AC_Possible = Num_TurningOn_AC-Number_AC_Shedded;

                if (isempty(TurningOn_AC_Indices)) % Check for using randsample()
                    TurnOn_AC_Possible_Indices=[];
                else
                    TurnOn_AC_Possible_Indices=randsample(TurningOn_AC_Indices,Num_TurnOn_AC_Possible);                      
                end        
                TurnOn_AC_NotPossible_Indices=TurningOn_AC_Indices(~ismember(TurningOn_AC_Indices,TurnOn_AC_Possible_Indices));

                % Turn On AC
                U_k(1,3,TurnOn_AC_Possible_Indices)=U_k_DumbLocalController(1,3,TurnOn_AC_Possible_Indices);

                % Didnt Turn On AC
                U_k(1,3,TurnOn_AC_NotPossible_Indices)=zeros(1,1,length(TurnOn_AC_NotPossible_Indices));

                % Not Turning On AC (Rest of the ACs)
                U_k(1,3,Not_TurningOn_AC_Indices)=U_k_DumbLocalController(1,3,Not_TurningOn_AC_Indices);

            else

                Num_TurnOn_AC_Possible=min([floor(Total_StartUpPower_Available/(ACLoad_StartUp_Power)),Num_TurningOn_AC-Number_AC_Shedded]);

                if (isempty(TurningOn_AC_Indices)) % Check for using randsample()
                    TurnOn_AC_Possible_Indices=[];
                else
                    TurnOn_AC_Possible_Indices=randsample(TurningOn_AC_Indices,Num_TurnOn_AC_Possible);                      
                end 
                TurnOn_AC_NotPossible_Indices=TurningOn_AC_Indices(~ismember(TurningOn_AC_Indices,TurnOn_AC_Possible_Indices));

                % Turn On AC
                U_k(1,3,TurnOn_AC_Possible_Indices)=U_k_DumbLocalController(1,3,TurnOn_AC_Possible_Indices);

                % Didnt Turn On AC
                U_k(1,3,TurnOn_AC_NotPossible_Indices)=zeros(1,1,length(TurnOn_AC_NotPossible_Indices));

                % Not Turning On AC (Rest of the ACs)
                U_k(1,3,Not_TurningOn_AC_Indices)=U_k_DumbLocalController(1,3,Not_TurningOn_AC_Indices);

            end  

        else % E_Mis can be shedded through turning off turning on ACs and turned on ACs           

            Num_ExtraAC_TurnedOff=abs(Num_TurningOn_AC-Number_AC_Shedded);

            Num_TurnOn_AC_Possible=min([floor(Total_StartUpPower_Available/(ACLoad_StartUp_Power)),0]);

            if (isempty(TurningOn_AC_Indices)) % Check for using randsample()
                TurnOn_AC_Possible_Indices=[];
            else
                TurnOn_AC_Possible_Indices=randsample(TurningOn_AC_Indices,Num_TurnOn_AC_Possible);                      
            end
            TurnOn_AC_NotPossible_Indices=TurningOn_AC_Indices(~ismember(TurningOn_AC_Indices,TurnOn_AC_Possible_Indices));

            % Didnt Turn On AC (All Turning On ACs were not able to turn on)
            U_k(1,3,TurnOn_AC_NotPossible_Indices)=zeros(1,1,length(TurnOn_AC_NotPossible_Indices));

            if (isempty(TurnedOn_AC_Indices)) % Check for using randsample()
                TurnedOn_AC_Possible_Indices=[];
            else
                TurnedOn_AC_Possible_Indices=randsample(TurnedOn_AC_Indices,Num_TurnedOn_AC-Num_ExtraAC_TurnedOff);                      
            end 
            TurnedOn_AC_NotPossible_Indices=TurnedOn_AC_Indices(~ismember(TurnedOn_AC_Indices,TurnedOn_AC_Possible_Indices));

            % TurnedOn_AC_Possible_Indices AC (Rest of the ACs)
            U_k(1,3,TurnedOn_AC_Possible_Indices)=U_k_DumbLocalController(1,3,TurnedOn_AC_Possible_Indices);

            % TurnedOn_AC_NotPossible_Indices AC 
            U_k(1,3,TurnedOn_AC_NotPossible_Indices)=zeros(1,1,length(TurnedOn_AC_NotPossible_Indices));
 
            Remaing_AC_Indices=setdiff(1:N_House,union(TurningOn_AC_Indices,TurnedOn_AC_Indices));

            % Rest of the ACs
            U_k(1,3,Remaing_AC_Indices)=U_k_DumbLocalController(1,3,Remaing_AC_Indices);

        end     

        
   elseif (abs(E_Mis)<=(sum(E_Bat_Charging_Dispatch(1,1,BatteryCharging_Indices))+((E_AC/Eff_Inv)*(sum(U_k_DumbLocalController(1,3,:))))+sum(E_Load_Desired_DC))) % E_Mis can be shedded through making charging batteries idle, turning off ACs and E_Load_Desired

        % Computing Battery Discharging Control Commands (All Charging Batteries are shed)
        U_k(1,2,:)=U_k_DumbLocalController(1,2,:);
        
        % Computing Battery Charging Control Commands (All Charging Batteries are shed)
        U_k(1,1,:)=zeros(1,1,N_House);  

        % All ACs are shed
        U_k(1,3,:)=zeros(1,1,N_House);      

        % Computing Actual Load Serviced Non-AC
        [U_k_PriorityStack] = PriorityStackController_SmartCommunity(E_LoadData_DC,(abs(E_Mis)-(sum(E_Bat_Charging_Dispatch(1,1,BatteryCharging_Indices))+((E_AC/Eff_Inv)*(sum(U_k_DumbLocalController(1,3,:)))))));   
        
        % All Loads are Shed
        U_k(1,4:end,:)=U_k_PriorityStack;

   else % E_Mis cannot be shedded through both turning off AC's and E_Load_Desired

        % Computing Battery Discharging Control Commands (All Charging Batteries are shed)
        U_k(1,2,:)=U_k_DumbLocalController(1,2,:);
        
        % Computing Battery Charging Control Commands (All Charging Batteries are shed)
        U_k(1,1,:)=zeros(1,1,N_House);  

        % All ACs are Shed
        U_k(1,3,:)=zeros(1,1,N_House); 
        
        % All Loads are Shed
        U_k(1,4:end,:)=zeros(1,8,N_House);

   end

end

end



