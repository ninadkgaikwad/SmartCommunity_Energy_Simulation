function [Plant_Performance]=HEMS_Plant_Performance_Computer(HEMS_PerformanceComputation)

% Author: Ninad Kiran Gaikwad
% Date: Mar/15/2021
% Description: HEMS_Plant_Performance_Computer - Computing Performance Measure of the Simulation

%% HEMS_Plant_Performance_Computer - Computing Performance Measure of the Simulation

%% Getting the desired data from the HEMS_Plant_Baseline_FigurePlotter_Input - Struct

%---------------------HEMS_PerformanceComputation-------------------------%
X_k_Plant_History=HEMS_PerformanceComputation.X_k_Plant_History;
U_k_History=HEMS_PerformanceComputation.U_k_History; 
E_LoadData=HEMS_PerformanceComputation.E_LoadData;
E_Load_Desired=HEMS_PerformanceComputation.E_Load_Desired;
HEMSWeatherData_Output=HEMS_PerformanceComputation.HEMSWeatherData_Output;
HEMSPlant_Params=HEMS_PerformanceComputation.HEMSPlant_Params;
Community_Params=HEMS_PerformanceComputation.Community_Params;

%-----------------------HEMSWeatherData_Output----------------------------%
Ws=HEMSWeatherData_Output.Ws;
T_am=HEMSWeatherData_Output.T_am;
GHI=HEMSWeatherData_Output.GHI;
DNI=HEMSWeatherData_Output.DNI;
DateTimeVector=HEMSWeatherData_Output.DateTimeVector;
DateTime_Matrix=HEMSWeatherData_Output.DateTime_Matrix;

%---------------------------HEMSPlant_Params------------------------------%

E_AC=HEMSPlant_Params.E_AC;
T_AC_max=HEMSPlant_Params.T_AC_max;
T_AC_min=HEMSPlant_Params.T_AC_min;
ACLoad_StartUp_Power=HEMSPlant_Params.ACLoad_StartUp_Power;
Eff_Inv=HEMSPlant_Params.Eff_Inv;

Battery_Energy_Max=HEMSPlant_Params.Battery_Energy_Max;
Battery_Energy_Min=HEMSPlant_Params.Battery_Energy_Min;
MaxRate_Charging=HEMSPlant_Params.MaxRate_Charging;
MaxRate_Discharging=HEMSPlant_Params.MaxRate_Discharging;
Eff_Charging_Battery=HEMSPlant_Params.Eff_Charging_Battery;
Eff_Discharging_Battery=HEMSPlant_Params.Eff_Discharging_Battery;
MaxRate_Discharging_StartUp=HEMSPlant_Params.MaxRate_Discharging_StartUp;

%---------------------------Community_Params.-----------------------------%

N_House=Community_Params.N_House;
N_PV_Bat=Community_Params.N_PV_Bat;
N_Bat=Community_Params.N_Bat;
N_PV=Community_Params.N_PV;
N_None=Community_Params.N_None; 

%% Basic Computation

% House Numbers
N1=N_PV_Bat;
N2=N_Bat;
N3=N_PV;
N4=N_None;

% Truncating Plant History
X_k_Plant_History=X_k_Plant_History(1:end-1,:,:);

% Computing Number of Days of Simulation
Start_DateTime=DateTimeVector(1);
End_DateTime=DateTimeVector(end);

TotalNum_Days_Simulation=daysact(Start_DateTime,End_DateTime);

%% Performance Metric Computation

% For All Houses AC
for jj=1:N_House
    
    % For One House AC
    AC_Death_TimeTotal=0;
    for ii=1:length(X_k_Plant_History(:,7,jj))
        if(X_k_Plant_History(ii,7,jj)>(T_AC_max+2))
            AC_Death_TimeTotal=AC_Death_TimeTotal+(10/60); % Hours
        end
    end
    AC_Death_AvgPerDay(jj)=AC_Death_TimeTotal/(TotalNum_Days_Simulation);

end

% For All Houses All Other Loads
All_Served=X_k_Plant_History(:,12,:);
All_Desired=E_Load_Desired(:,1,:);

for jj=1:N_House
    
    % For One House All Other Loads
    All_Desired_Points=0;
    All_Served_Points=0;
    for ii=1:length(All_Served)
        if(All_Desired(ii,1,jj)>0)
            All_Desired_Points=All_Desired_Points+1;
        end
        if((All_Served(ii,1,jj)<All_Desired(ii,1,jj))&&(~(All_Served(ii,1,jj)<0)))
            % NC Load Not Served
        elseif ((All_Served(ii,1,jj)==All_Desired(ii,1,jj))&&(~(All_Served(ii,1,jj)<=0)))
            All_Served_Points=All_Served_Points+1;
        end
    end

    Percentage_All_Served(jj) = ((All_Served_Points)/(All_Desired_Points)*(100));

end

% For All Houses Critical Loads
C_Served=X_k_Plant_History(:,13,:);
C_Desired=E_LoadData(:,1,:);

for jj=1:N_House
    
    % For One House Critical Loads
    C_Desired_Points=0;
    C_Served_Points=0;
    for ii=1:length(C_Served)
        if(C_Desired(ii,1,jj)>0)
            C_Desired_Points=C_Desired_Points+1;
        end
        if((C_Served(ii,1,jj)<C_Desired(ii,1,jj))&&(~(C_Served(ii,1,jj)<0)))
            % NC Load Not Served
        elseif ((C_Served(ii,1,jj)==C_Desired(ii,1,jj))&&(~(C_Served(ii,1,jj)<=0)))
            C_Served_Points=C_Served_Points+1;
        end
    end

    Percentage_C_Served(jj) = ((C_Served_Points)/(C_Desired_Points)*(100));

end

% For All Houses
for jj=1:N_House
    
    % For One House PRM and SRM    
    PRM(jj)=1-(AC_Death_AvgPerDay(jj)/24);

    SRM_C(jj)=Percentage_C_Served(jj)/100;

    SRM_All(jj)=Percentage_All_Served(jj)/100;

end

%RP= (h*PRM) + ((1-h)*SRM_C);

% Computing Performance Measure for entire Community
AC_Death_AvgPerDay_Community=mean(AC_Death_AvgPerDay);
Percentage_All_Served_Community=mean(Percentage_All_Served);
Percentage_C_Served_Community=mean(Percentage_C_Served);
PRM_Community=mean(PRM);
SRM_C_Community=mean(SRM_C);
SRM_All_Community=mean(SRM_All);

%% Priniting Results

% Priniting Quantitative Results for each House
for ii=1:N_House
    
    fprintf('\n')
    fprintf('Average Fridge Death for House:%d (PLP) = %f Hours/Day',ii,AC_Death_AvgPerDay(ii));
    fprintf('\n')
    fprintf('Percentage of Non-Critical Loads served for House:%d (SLP) = %f %/Day',ii,Percentage_C_Served(ii));
    fprintf('\n')
    fprintf('Percentage of All Loads served for House:%d (SLP) = %f %/Day',ii,Percentage_All_Served(ii));
    fprintf('\n')
    fprintf(' House:%d (PRM) = %f Hours/Day',ii,PRM(ii));
    fprintf('\n')
    fprintf(' House:%d (SRM_C) = %f %/Day',ii,SRM_C(ii));
    fprintf('\n')
    fprintf(' House:%d (SRM_All) = %f %/Day',ii,SRM_All(ii));
    fprintf('\n')    
    fprintf('\n')

end


% Priniting Quantitative Results for Community
fprintf('Average Fridge Death for Community (PLP) = %f Hours/Day',AC_Death_AvgPerDay_Community);
fprintf('\n')
fprintf('Percentage of Non-Critical Loads served for Community (SLP) = %f %/Day',Percentage_C_Served_Community);
fprintf('\n')
fprintf('Percentage of All Loads served for Community (SLP) = %f %/Day',Percentage_All_Served_Community);
fprintf('\n')
fprintf(' Community (PRM) = %f Hours/Day',PRM_Community);
fprintf('\n')
fprintf(' Community (SRM_C) = %f %/Day',SRM_C_Community);
fprintf('\n')
fprintf(' Community (SRM_All) = %f %/Day',SRM_All_Community);
fprintf('\n')    
fprintf('\n')

%% Creating Plant_Performance

Plant_Performance = [];

Plant_Performance.AC_Death_AvgPerDay=AC_Death_AvgPerDay;
Plant_Performance.Percentage_All_Served=Percentage_All_Served;
Plant_Performance.Percentage_C_Served=Percentage_C_Served;
Plant_Performance.PRM=PRM_Community;
Plant_Performance.SRM_C=SRM_C_Community;
Plant_Performance.SRM_All=SRM_All;

Plant_Performance.AC_Death_AvgPerDay_Community=AC_Death_AvgPerDay_Community;
Plant_Performance.Percentage_All_Served_Community=Percentage_All_Served_Community;
Plant_Performance.Percentage_C_Served_Community=Percentage_C_Served_Community;
Plant_Performance.PRM_Community=PRM_Community;
Plant_Performance.SRM_C_Community=SRM_C_Community;
Plant_Performance.SRM_All_Community=SRM_All_Community;
    
end

