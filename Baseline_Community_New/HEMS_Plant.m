function [X_k_Plus_Plant] = HEMS_Plant(X_k_Plant,W_k_Plant,U_k,HEMSPlant_Params,HEMSHouse_Params,Community_Params,Simulation_Params)

% Author: Ninad Kiran Gaikwad
% Date: Mar/20/2021
% Description: HEMS_Plant - Plant Dynamics

%% HEMS_Plant - Plant Dynamics

%% Getting desired Data from Input - Structs

% From W_k_Plant
Weather_k_Plant=W_k_Plant.Weather_k_Plant;
LoadData_k_Plant=W_k_Plant.LoadData_k_Plant;

GHI=Weather_k_Plant.GHI;
DNI=Weather_k_Plant.DNI;
T_am=Weather_k_Plant.T_am;
Ws=Weather_k_Plant.Ws;
DateTime_Matrix(1,1:4)=Weather_k_Plant.DateTime_Matrix;

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

%---------------------- Computing Energy Mismatch ------------------------%

% Battery Energy Available for Discharging dispatch
E_Bat_Discharging_Dispatch=zeros(1,1,N_House);

for jj=1:N_PV_Bat+N_Bat % For Battery installed Houses
    
    E_Bat_Discharging_Dispatch(1,1,jj)=(U_k(1,2,jj))*...
        (min(MaxRate_Discharging*Simulation_StepSize,...
        (X_k_Plant(1,4,jj)-Battery_Energy_Min)*Eff_Discharging_Battery));
    
end

% Battery Energy Available for Charging dispatch
E_Bat_Charging_Dispatch=zeros(1,1,N_House);

for jj=1:N_PV_Bat+N_Bat % For Battery installed Houses
    
    E_Bat_Charging_Dispatch(1,1,jj)=(U_k(1,1,jj))*...
        (min(MaxRate_Charging*Simulation_StepSize,...
        (Battery_Energy_Max-(X_k_Plant(1,4,jj)))/Eff_Charging_Battery));
    
end

% Battery Discharging Indices
BatteryDischarging_Indices=find(U_k(1,2,:)==1);
Battery_AboveMin_Indices=find(X_k_Plant(1,4,:)>Battery_Energy_Min);
Battery_ActualDischarging_Indices=intersect(BatteryDischarging_Indices,Battery_AboveMin_Indices);

% Battery Charging Indices
BatteryCharging_Indices=find(U_k(1,1,:)==1);
Battery_BelowMax_Indices=find(X_k_Plant(1,4,:)<Battery_Energy_Max);
Battery_ActualCharging_Indices=intersect(BatteryCharging_Indices,Battery_BelowMax_Indices);

% Computing Total energy desired by other loads (Not ACs)
for jj=1:N_House % For Each House
    
    OtherLoad_Energy_Desired(1,1,jj)=U_k(1,4:11,jj)*E_LoadData_DC(1,9+1:end,jj)';
    
end

Total_OtherLoad_Energy_Desired=sum(OtherLoad_Energy_Desired,3);

% Total energy Desired
TotalEnergy_Desired=((E_AC/Eff_Inv)*(sum(U_k(1,3,:))))+Total_OtherLoad_Energy_Desired+sum(E_Bat_Charging_Dispatch(1,1,Battery_ActualCharging_Indices));

% Energy Mismatch
AvailableEnergy=Total_PVEnergy_Available+sum(E_Bat_Discharging_Dispatch(1,1,Battery_ActualDischarging_Indices));
E_Mis=AvailableEnergy-TotalEnergy_Desired;

%------------ Computing Total Startup Power Required for ACs -------------%

% Computing Turning On ACs
TurningOn_AC_Vector=U_k(1,3,:)-X_k_Plant(1,30,:);
TurningOn_AC_Indices=find(TurningOn_AC_Vector==1);
Num_TurningOn_AC=length(TurningOn_AC_Indices);

% Computing Quantities for Energy and Power Constraints
Battery_CapableofDischarge_Indices=intersect(find(U_k(1,2,:)),find(X_k_Plant(1,4,:)>Battery_Energy_Min));
TotalPower_Available_Bat=length(Battery_CapableofDischarge_Indices)*(MaxRate_Discharging_StartUp);
TotalPower_Available_PV=(Total_PVEnergy_Available)/(Simulation_StepSize);
Total_StartUpPower_Available=TotalPower_Available_Bat+TotalPower_Available_PV;
Total_StartUpPower_AC=Num_TurningOn_AC*ACLoad_StartUp_Power;

%% Computing PV Energy Available for each House with PV

for jj=1:N_House

    % Computing PV Energy Available
    if((jj<=N_PV_Bat)||((jj>(N_PV_Bat+N_Bat))&&(jj<=(N_PV_Bat+N_Bat+N_PV)))) % For PV installed Houses
        X_k_Plant(1,1,jj)=PVEnergy_Available;
    end

end   

%% Setting up Weather Dtata Struct for House Thermal Dynamics Simulation

% Setting up HEMSWeatherData_Output1
HEMSWeatherData_Output1.Ws=Ws;
HEMSWeatherData_Output1.T_am=T_am;
HEMSWeatherData_Output1.GHI=GHI;
HEMSWeatherData_Output1.DNI=DNI;
HEMSWeatherData_Output1.DateTime_Matrix=DateTime_Matrix;

%% Initializing X_k_Plus_Plant

X_k_Plus_Plant=X_k_Plant;

%% Plant Physics

if (E_Mis>=0) % All Loads can be serviced

    % Computing Actual Load Serviced AC
    if(Total_StartUpPower_Available>=Total_StartUpPower_AC) % Enough Batteries/PV to Power the Turn On of ACs
        
        % Energy Consumed by AC and Other Loads 
        X_k_Plus_Plant(1,11:20,:) = [(E_AC/Eff_Inv)*U_k(1,3,:)...
                                      OtherLoad_Energy_Desired...
                                      U_k(1,4:11,:).*E_LoadData_DC(1,9+1:end,:)];  

        % AC and Other Loads Current On-Off Status 
        X_k_Plus_Plant(1,21:29,:) = U_k(1,3:11,:);

        % AC and Other Loads Previous On-Off Status (None On)
        X_k_Plus_Plant(2,30:38,:) = X_k_Plus_Plant(1,21:29,:);  
        
        % Computing Charging Load
        Total_ChargingEnergy_Available=-(sum(X_k_Plus_Plant(1,11,:))+sum(X_k_Plus_Plant(1,12,:)))+sum(X_k_Plus_Plant(1,1,:))+sum(E_Bat_Discharging_Dispatch(1,1,BatteryDischarging_Indices));

        if (Total_ChargingEnergy_Available>0) % Energy Available for charging

            Average_BatteryCharge=min([Total_ChargingEnergy_Available/length(Battery_ActualCharging_Indices),sum(E_Bat_Charging_Dispatch(1,1,BatteryCharging_Indices)/length(Battery_ActualCharging_Indices))]);

            [~,Ascending_BatteryLevel_ChargingBat_Indices_1]=sort(X_k_Plant(1,4,:));

            Counter_ChargingBat=0;

            if (~isempty(Battery_ActualCharging_Indices))

                for kk=1:length(Ascending_BatteryLevel_ChargingBat_Indices_1)    

                    for mm=1:length(Battery_ActualCharging_Indices)  

                        if (Battery_ActualCharging_Indices(mm)==Ascending_BatteryLevel_ChargingBat_Indices_1(kk))

                            Counter_ChargingBat=Counter_ChargingBat+1;

                            Ascending_BatteryLevel_ChargingBat_Indices(Counter_ChargingBat)=Battery_ActualCharging_Indices(mm);

                        end

                    end

                end  

            else

                Ascending_BatteryLevel_ChargingBat_Indices=[];

            end             

            BatteryCharge_Sum=0; % Initialization

            BatteryCharge_Counter=0; % Initialization

            Total_BatteryCharge_RequiredNow=Total_ChargingEnergy_Available; % Initialization

            for jj=Ascending_BatteryLevel_ChargingBat_Indices % For each Charging Battery with Ascending battery level

                if (E_Bat_Charging_Dispatch(1,1,jj)>Average_BatteryCharge)

                    X_k_Plus_Plant(1,5,jj)=min([Average_BatteryCharge+((E_Bat_Charging_Dispatch(1,1,jj)-Average_BatteryCharge)*((Battery_Energy_Max-X_k_Plus_Plant(1,4,jj))/(Battery_Energy_Max-Battery_Energy_Min))),...
                        Total_BatteryCharge_RequiredNow, (Battery_Energy_Max-X_k_Plus_Plant(1,4,jj))/(Eff_Charging_Battery),MaxRate_Charging*Simulation_StepSize]);

                else

                    X_k_Plus_Plant(1,5,jj)=min([Average_BatteryCharge-((Average_BatteryCharge-E_Bat_Charging_Dispatch(1,1,jj))*(1-((Battery_Energy_Max-X_k_Plus_Plant(1,4,jj))/(Battery_Energy_Max-Battery_Energy_Min)))),...
                        Total_BatteryCharge_RequiredNow, (Battery_Energy_Max-X_k_Plus_Plant(1,4,jj))/(Eff_Charging_Battery),MaxRate_Charging*Simulation_StepSize]);

                end

                if (X_k_Plus_Plant(1,5,jj)>0)

                    % Incrementing BatteryDischarge_Counter
                    BatteryCharge_Counter=BatteryCharge_Counter+1;  

                end

                % Incrementing BatteryDischarge_Sum
                %BatteryCharge_Sum=BatteryCharge_Sum+House_DataMatrix(ii,5,jj); 

                % Computing new Total_BatteryCharge_RequiredNow
                %Total_BatteryCharge_RequiredNow=Total_BatteryCharge_RequiredNow-BatteryCharge_Sum;

                Total_BatteryCharge_RequiredNow=Total_BatteryCharge_RequiredNow-X_k_Plus_Plant(1,5,jj);

                if (Total_BatteryCharge_RequiredNow<0) % Total_BatteryCharge_RequiredNow cannot be less than 0

                    Total_BatteryCharge_RequiredNow=0;

                end

                % Computing new Average_BatteryDischarge
                %Average_BatteryDischarge=(Total_BatteryDischarge_Required-BatteryDischarge_Sum)/(length(BatteryDischarging_Indices)-BatteryDischarge_Counter);

            end        

        else % Energy not Available for charging

            X_k_Plus_Plant(1,5,:)=zeros(1,1,N_House);        

        end         

        % Computing Battery Discharge Required
        Total_BatteryDischarge_Required=(sum(X_k_Plus_Plant(1,11,:))+sum(X_k_Plus_Plant(1,12,:)))+sum(X_k_Plus_Plant(1,5,BatteryCharging_Indices))-sum(X_k_Plus_Plant(1,1,:));

        if (Total_BatteryDischarge_Required>0) % Battery discharge required to service loads

            Average_BatteryDischarge=min([Total_BatteryDischarge_Required/length(Battery_ActualDischarging_Indices),sum(E_Bat_Discharging_Dispatch(1,1,BatteryDischarging_Indices))/length(Battery_ActualDischarging_Indices)]);

            [~,Sorted_BatteryLevel_DischargingBat_Indices]=sort(X_k_Plant(1,4,:));
            [Descending_BatteryLevel_DischargingBat_Indices_1]=flip(Sorted_BatteryLevel_DischargingBat_Indices);

            Counter_DichargingBat=0;

            if (~isempty(Battery_ActualDischarging_Indices))

                for kk=1:length(Sorted_BatteryLevel_DischargingBat_Indices)    

                    for mm=1:length(Battery_ActualDischarging_Indices)

                        if (Battery_ActualDischarging_Indices(mm)==Descending_BatteryLevel_DischargingBat_Indices_1(kk))

                            Counter_DichargingBat=Counter_DichargingBat+1;

                            Descending_BatteryLevel_DischargingBat_Indices(Counter_DichargingBat)=Battery_ActualDischarging_Indices(mm);

                        end

                    end

                end  

            else

                Descending_BatteryLevel_DischargingBat_Indices=[];

            end                

            BatteryDischarge_Sum=0; % Initialization

            BatteryDischarge_Counter=0; % Initialization

            Total_BatteryDischarge_RequiredNow=Total_BatteryDischarge_Required; % Initialization

            for jj=Descending_BatteryLevel_DischargingBat_Indices % For each Discharging Battery with descending battery level

                if (E_Bat_Discharging_Dispatch(1,1,jj)>Average_BatteryDischarge)

                    X_k_Plus_Plant(1,6,jj)=min([Average_BatteryDischarge+((E_Bat_Discharging_Dispatch(1,1,jj)-Average_BatteryDischarge)*(1-((Battery_Energy_Max-X_k_Plus_Plant(1,4,jj))/(Battery_Energy_Max-Battery_Energy_Min)))),...
                        Total_BatteryDischarge_RequiredNow, (X_k_Plus_Plant(1,4,jj)-Battery_Energy_Min)*(Eff_Discharging_Battery),MaxRate_Discharging*Simulation_StepSize]);

                else

                    X_k_Plus_Plant(1,6,jj)=min([Average_BatteryDischarge-((Average_BatteryDischarge-E_Bat_Discharging_Dispatch(1,1,jj))*((Battery_Energy_Max-X_k_Plus_Plant(1,4,jj))/(Battery_Energy_Max-Battery_Energy_Min))),...
                        Total_BatteryDischarge_RequiredNow, (X_k_Plus_Plant(1,4,jj)-Battery_Energy_Min)*(Eff_Discharging_Battery),MaxRate_Discharging*Simulation_StepSize]);

                end

                if (X_k_Plus_Plant(1,6,jj)>0)

                    % Incrementing BatteryDischarge_Counter
                    BatteryDischarge_Counter=BatteryDischarge_Counter+1;  

                end

                % Incrementing BatteryDischarge_Sum
                %BatteryDischarge_Sum=BatteryDischarge_Sum+House_DataMatrix(ii,6,jj); 

                % Computing new Total_BatteryDischarge_RequiredNow
                %Total_BatteryDischarge_RequiredNow=Total_BatteryDischarge_RequiredNow-BatteryDischarge_Sum;

                Total_BatteryDischarge_RequiredNow=Total_BatteryDischarge_RequiredNow-X_k_Plus_Plant(1,6,jj);

                if (Total_BatteryDischarge_RequiredNow<0) % Total_BatteryDischarge_RequiredNow cannot be less than 0

                    Total_BatteryDischarge_RequiredNow=0;

                end                    

                % Computing new Average_BatteryDischarge
                %Average_BatteryDischarge=(Total_BatteryDischarge_Required-BatteryDischarge_Sum)/(length(BatteryDischarging_Indices)-BatteryDischarge_Counter);

            end

        else % Battery Discharge not required to service loads

            X_k_Plus_Plant(1,6,:)=zeros(1,1,N_House);

        end

        % Computing New Battery State
        X_k_Plus_Plant(2,4,:)=X_k_Plus_Plant(1,4,:)+(X_k_Plus_Plant(1,5,:)*Eff_Charging_Battery)-(X_k_Plus_Plant(1,6,:)/Eff_Discharging_Battery);

        % Computing PV Used 
        Total_PVEnergy_Used_Possible = (sum(X_k_Plus_Plant(1,11,:))+sum(X_k_Plus_Plant(1,12,:))) + sum(X_k_Plus_Plant(1,5,:)) - sum(X_k_Plus_Plant(1,6,:));

        if (Total_PVEnergy_Used_Possible>0)

            Average_PvEnergy_Used = Total_PVEnergy_Used_Possible/(N_PV_Bat+N_PV);

            X_k_Plus_Plant(1,2,(1:N_PV_Bat))= Average_PvEnergy_Used*ones(1,1,length((1:N_PV_Bat))); 
            X_k_Plus_Plant(1,2,(N_PV_Bat+N_Bat+1:N_PV_Bat+N_Bat+N_PV))= Average_PvEnergy_Used*ones(1,1,length((N_PV_Bat+N_Bat+1:N_PV_Bat+N_Bat+N_PV)));

            % Computing PV Unused Energy
            X_k_Plus_Plant(1,3,:)=X_k_Plus_Plant(1,1,:)-X_k_Plus_Plant(1,2,:);  

            %PV_Unused_Negative=length(find(X_k_Plus_Plant(1,3,:)<0)); % Debugger

            %PV_Unused_Negative1=PV_Unused_Negative;

        else
            % Nothing
        end  
        
        % House Thermal Dynamics Simulation

        for jj=1:N_House % For each House in Smart Comuunity       

            % Creating 
            HEMSHouse_States=[]; % Initialization

            % Setting up House States initial conditions for the first time

            HEMSHouse_States.T_wall1=X_k_Plant(1,8,jj);
            HEMSHouse_States.T_ave1=X_k_Plant(1,7,jj);
            HEMSHouse_States.T_attic1=X_k_Plant(1,9,jj);
            HEMSHouse_States.T_im1=X_k_Plant(1,10,jj);        

            % Getting Correct Actual AC ON-OFF for this time Step for the House
            HEMSHouse_States.u_k_hvac=X_k_Plus_Plant(1,21,jj);

            % Calling external function for House Thermal Dynamics
            [HEMSHouseRCModel_Output1] = HEMS_HouseRCModel(HEMSHouse_Params,HEMSWeatherData_Output1,Simulation_Params,HEMSPlant_Params,HEMSHouse_States);

            % Setting up House States
            HEMSHouse_States.T_wall1=HEMSHouseRCModel_Output1.T_wall(end);
            HEMSHouse_States.T_ave1=HEMSHouseRCModel_Output1.T_ave(end);
            HEMSHouse_States.T_attic1=HEMSHouseRCModel_Output1.T_attic(end);
            HEMSHouse_States.T_im1=HEMSHouseRCModel_Output1.T_im(end);

            % Updating House_DataMatrix with House Thermal States
            X_k_Plus_Plant(2,7:10,jj)=[HEMSHouseRCModel_Output1.T_ave(end),HEMSHouseRCModel_Output1.T_wall(end),HEMSHouseRCModel_Output1.T_attic(end),HEMSHouseRCModel_Output1.T_im(end)];

        end            

    else % Not Enough Batteries/PV to Power the Turn On of ACs - Nothing will work

        % PV Energy Used and PV Energy Unused
        X_k_Plus_Plant(1,2,:) = zeros(1,1,N_House); % PV Energy Used
        X_k_Plus_Plant(1,3,:) = X_k_Plant(1,1,:); % PV Energy Unused

        % Battery Energy State/Battery Charging Energy/Battery Discharging Energy
        X_k_Plus_Plant(2,4,:) = X_k_Plant(1,4,:); % Battery Energy State
        X_k_Plus_Plant(1,5,:) = zeros(1,1,N_House); % Battery Charging Energy 
        X_k_Plus_Plant(1,6,:) = zeros(1,1,N_House); % Battery Discharging Energy 

        % Energy Consumed by AC and Other Loads (None Consumed)
        X_k_Plus_Plant(1,11:20,:) = zeros(1,10,N_House);

        % AC and Other Loads Current On-Off Status (None On)
        X_k_Plus_Plant(1,21:29,:) = zeros(1,9,N_House);

        % AC and Other Loads Previous On-Off Status (None On)
        X_k_Plus_Plant(2,30:38,:) = X_k_Plus_Plant(1,21:29,:);    

        % House Thermal Dynamics Simulation

        for jj=1:N_House % For each House in Smart Comuunity       

            % Creating 
            HEMSHouse_States=[]; % Initialization

            % Setting up House States initial conditions for the first time

            HEMSHouse_States.T_wall1=X_k_Plant(1,8,jj);
            HEMSHouse_States.T_ave1=X_k_Plant(1,7,jj);
            HEMSHouse_States.T_attic1=X_k_Plant(1,9,jj);
            HEMSHouse_States.T_im1=X_k_Plant(1,10,jj);        

            % Getting Correct Actual AC ON-OFF for this time Step for the House
            HEMSHouse_States.u_k_hvac=X_k_Plus_Plant(1,21,jj);

            % Calling external function for House Thermal Dynamics
            [HEMSHouseRCModel_Output1] = HEMS_HouseRCModel(HEMSHouse_Params,HEMSWeatherData_Output1,Simulation_Params,HEMSPlant_Params,HEMSHouse_States);

            % Setting up House States
            HEMSHouse_States.T_wall1=HEMSHouseRCModel_Output1.T_wall(end);
            HEMSHouse_States.T_ave1=HEMSHouseRCModel_Output1.T_ave(end);
            HEMSHouse_States.T_attic1=HEMSHouseRCModel_Output1.T_attic(end);
            HEMSHouse_States.T_im1=HEMSHouseRCModel_Output1.T_im(end);

            % Updating House_DataMatrix with House Thermal States
            X_k_Plus_Plant(2,7:10,jj)=[HEMSHouseRCModel_Output1.T_ave(end),HEMSHouseRCModel_Output1.T_wall(end),HEMSHouseRCModel_Output1.T_attic(end),HEMSHouseRCModel_Output1.T_im(end)];

        end    

    end

elseif (E_Mis<0) % All loads cannot be serviced

    % PV Energy Used and PV Energy Unused
    X_k_Plus_Plant(1,2,:) = zeros(1,1,N_House); % PV Energy Used
    X_k_Plus_Plant(1,3,:) = X_k_Plant(1,1,:); % PV Energy Unused
    
    % Battery Energy State/Battery Charging Energy/Battery Discharging Energy
    X_k_Plus_Plant(2,4,:) = X_k_Plant(1,4,:); % Battery Energy State
    X_k_Plus_Plant(1,5,:) = zeros(1,1,N_House); % Battery Charging Energy 
    X_k_Plus_Plant(1,6,:) = zeros(1,1,N_House); % Battery Discharging Energy 
    
    % Energy Consumed by AC and Other Loads (None Consumed)
    X_k_Plus_Plant(1,11:20,:) = zeros(1,10,N_House);
    
    % AC and Other Loads Current On-Off Status (None On)
    X_k_Plus_Plant(1,21:29,:) = zeros(1,9,N_House);
    
    % AC and Other Loads Previous On-Off Status (None On)
    X_k_Plus_Plant(2,30:38,:) = X_k_Plus_Plant(1,21:29,:);    
   
    % House Thermal Dynamics Simulation

    for jj=1:N_House % For each House in Smart Comuunity       
              
        % Creating 
        HEMSHouse_States=[]; % Initialization

        % Setting up House States initial conditions for the first time

        HEMSHouse_States.T_wall1=X_k_Plant(1,8,jj);
        HEMSHouse_States.T_ave1=X_k_Plant(1,7,jj);
        HEMSHouse_States.T_attic1=X_k_Plant(1,9,jj);
        HEMSHouse_States.T_im1=X_k_Plant(1,10,jj);        
        
        % Getting Correct Actual AC ON-OFF for this time Step for the House
        HEMSHouse_States.u_k_hvac=X_k_Plus_Plant(1,21,jj);
        
        % Calling external function for House Thermal Dynamics
        [HEMSHouseRCModel_Output1] = HEMS_HouseRCModel(HEMSHouse_Params,HEMSWeatherData_Output1,Simulation_Params,HEMSPlant_Params,HEMSHouse_States);

        % Setting up House States
        HEMSHouse_States.T_wall1=HEMSHouseRCModel_Output1.T_wall(end);
        HEMSHouse_States.T_ave1=HEMSHouseRCModel_Output1.T_ave(end);
        HEMSHouse_States.T_attic1=HEMSHouseRCModel_Output1.T_attic(end);
        HEMSHouse_States.T_im1=HEMSHouseRCModel_Output1.T_im(end);

        % Updating House_DataMatrix with House Thermal States
        X_k_Plus_Plant(2,7:10,jj)=[HEMSHouseRCModel_Output1.T_ave(end),HEMSHouseRCModel_Output1.T_wall(end),HEMSHouseRCModel_Output1.T_attic(end),HEMSHouseRCModel_Output1.T_im(end)];

    end    

end

end



