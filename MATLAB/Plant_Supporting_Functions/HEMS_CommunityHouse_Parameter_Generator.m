function [HEMSPlant_Params,HEMSHouse_Params] = HEMS_CommunityHouse_Parameter_Generator(Community_Params,Simulation_Params)

% Author: Ninad Kiran Gaikwad
% Date: Feb/10/2021
% Description: HEMS_CommunityHouse_Parameter_Generator - Parameter Generation

%% Getting Variables from Input Structs

% From Community Params
N_House=Community_Params.N_House;
N_PV_Bat=Community_Params.N_PV_Bat;
N_Bat=Community_Params.N_Bat;
N_PV=Community_Params.N_PV;
N_None=Community_Params.N_None;  

% From Simulation Params
FileRes=Simulation_Params.FileRes;
StepSize=Simulation_Params.StepSize;
Simulation_StepSize=Simulation_Params.Simulation_StepSize;

%% HEMS_CommunityHouse_Parameter_Generator - Parameter Generation 

%-------------------------- Location Information -------------------------%

Location_Name='Gainesville_Fl_USA';
hem=-1;        % 1 - Eastern Hemisphere ; -1 - Western Hemisphere
Lat=18.17;     % Decimal Deg
Long= 66.74;   % Decimal Deg
TimeZone=4;    % Hours
Ltm=60;        % Deg Decimal - Time Meridian (Atlantic Standard Time)

%----------------------- House Thermal Specifications --------------------%

% Thermal Resistances (K/W)
R_w=0.0134;             % External Wall
R_attic=0.0235;         % Attic
R_roof=0.00156;         % Roof
R_im=0.00171;           % Internal Mass
R_win=0.021;            % Windows

% Thermal Capacitances (J/K)
C_w=10383364;           % External Wall
C_attic=704168;         % Attic
C_im=23396403;          % Internal Mass
C_in=8665588;           % Intdoor Air
C1=0.691;
C2=0.784;
C3=0.1;

% Internal Heat Load [W]
Human_Num=4;
Human_Heat=100;         % 2000 kCal/Day = 2000*10^3*4.18 J/Day = 96.91 W
Appliance_Heat=60+65+100;

Q_ihl=Human_Num*Human_Heat+Appliance_Heat; % [W] Includes Humans and all internal appliances which contribute to the Heat gain/loss of the house

Q_ac=0;                 % [W] AC Load

% Ventilation and Infiltration Heat Load
Cp=1;                   % @ 300K, 1atm [kJ/kg K] Specific Heat of Air
V=0;                    % Mass Flow Rate (m^3/s)
Den_Air=1.18;           % @ 300K, 1atm [kg/m^3] Air Density
C_oew=0.01;             % Coefficient used to scale the wind speed to the infiltration rate

% Q_solar
SHGC=0.8;               % From ASHRAE - Fundamentals

% Windows Roofs and Walls Orientation
Alpha_w=0.9;            % Radiation absorption Coefficient
Alpha_r=0.9;            % Radiation absorption Coefficient
% Wall
Area_w=[33,33,33,33];   % m^(2)
Tilt_w=[90,90,90,90];   % Deg (Angle Between wall Surface and Ground Surface measured from Ground Surface)
Azi_w=[0,90,-90,180];   % Deg (Angle of the wall Surface in the E-W plane wrt N-S axis [S-E-N -- 0-90-180 ; S-W-N -- 0-(-90)-(-180) ] )
% Roof
Area_r=[100];           % m^(2)
Tilt_r=[0];             % Deg (Angle Between roof Surface and Ground Surface)
Azi_r=[0];              % Deg (Angle of the roof Surface in the E-W plane wrt N-S axis [S-E-N -- 0-90-180 ; S-W-N -- 0-(-90)-(-180) ] )
% Window
Area_win=[10,10,10,10]; % m^(2)
Tilt_win=[90,90,90,90]; % Deg (Angle Between window Surface and Ground Surface)
Azi_win=[0,90,-90,180]; % Deg (Angle of the window Surface in the E-W plane wrt N-S axis [S-E-N -- 0-90-180 ; S-W-N -- 0-(-90)-(-180) ] )

%-------------------------- AC Specifications --------------------------%

% Controller Load
ACLoad_Num=1;
ACLoad_Power=3000; % Tesla Site Value

AC_VolatageDip=0.3; % 30% Volage Dip at Startup causes 30% Current Dip
AC_StartUp_LRA_Factor=6; % LRA - Locked Rotor Current usually 3-8 times rated current
ACLoad_StartUp_Power=(1-AC_VolatageDip)^(2)*AC_StartUp_LRA_Factor*ACLoad_Power/1000; % In kW

AC_COP=3.33;
T_AC_Base=24;
T_AC_DeadBand = 1;

T_AC_max=T_AC_Base+T_AC_DeadBand;
T_AC_min=T_AC_Base-T_AC_DeadBand;

T_AC= T_AC_Base;

AC_Indicator=1;
AC_Indicator1=1;

E_AC=ACLoad_Num*ACLoad_Power*Simulation_StepSize*(1/1000);

%-------------------------- PV Specifications --------------------------%

% Tesla Module SC325

PV_TotlaModules_Num=31;  % User Input - 10kW system -- 31 Modules
PV_RatedPower=325;      % User Input Watts
PV_TempCoeff=-0.31;     % User Input %/DegC
GHI_Std = 1000;         % kW/m^2
Temp_Std = 25;          % Deg C
Eff_Inv=0.9;            % Inverter Efficiency
   
% Faiman Model for Module Temperature Computation
Uo=25;                  % Values of U_{0} varied from 23.5 to 26.5 with a combined fit = 25 W/m^{2}K
U1=6.84 ;               % Values of U_{1} varied from 6.25 to 7.68 with a combined fit = 6.84 W/m^{2}K

PV_Total_Power = PV_TotlaModules_Num * PV_RatedPower;

%-------------------------- Battery Specifications --------------------------%

           % User Input Volts
DOD=0;                        % Depth of Discharge (For logetivity)
Eff_Charging_Battery=0.9;       % Battery Charging Efficiency
Eff_Discharging_Battery=0.9;    % Battery Discharging Efficiency
N1=1;                           % User Input - Battery Max Changing Factor

MaxRate_Charging=7.6;

P_AC=(ACLoad_Num*ACLoad_Power)*(1/1000);

% MaxRate_Discharging=Bat_TotalStrings_Num*(Bat_Capacity2/C_Hours2)*(System_Voltage)*(1/1000); % kWatts
MaxRate_Discharging=5;
MaxRate_Discharging_StartUp=7;

Battery_Energy_Max = 13.5*N1;
Battery_Energy_Min=Battery_Energy_Max*DOD;
 
HEMSPlant_Params=[]; % Empty Struct

HEMSPlant_Params.Location_Name=Location_Name;
HEMSPlant_Params.hem=hem;
HEMSPlant_Params.Lat=Lat;
HEMSPlant_Params.Long=Long;
HEMSPlant_Params.TimeZone=TimeZone;
HEMSPlant_Params.Ltm=Ltm;

HEMSPlant_Params.ACLoad_Num=ACLoad_Num;
HEMSPlant_Params.ACLoad_Power=ACLoad_Power;
HEMSPlant_Params.AC_COP=AC_COP;
HEMSPlant_Params.T_AC_Base=T_AC_Base;
HEMSPlant_Params.T_AC_DeadBand=T_AC_DeadBand;
HEMSPlant_Params.AC_Indicator=AC_Indicator;
HEMSPlant_Params.AC_Indicator1=AC_Indicator1;
HEMSPlant_Params.E_AC=E_AC;
HEMSPlant_Params.T_AC_max=T_AC_max;
HEMSPlant_Params.T_AC_min=T_AC_min;
HEMSPlant_Params.ACLoad_StartUp_Power=ACLoad_StartUp_Power;

HEMSPlant_Params.PV_TotlaModules_Num=PV_TotlaModules_Num;
HEMSPlant_Params.PV_RatedPower=PV_RatedPower;
HEMSPlant_Params.PV_TempCoeff=PV_TempCoeff;
HEMSPlant_Params.GHI_Std=GHI_Std;
HEMSPlant_Params.Temp_Std=Temp_Std;
HEMSPlant_Params.Eff_Inv=Eff_Inv;
HEMSPlant_Params.Uo=Uo;
HEMSPlant_Params.U1=U1;

HEMSPlant_Params.DOD=DOD;
HEMSPlant_Params.Eff_Charging_Battery=Eff_Charging_Battery;
HEMSPlant_Params.Eff_Discharging_Battery=Eff_Discharging_Battery;
HEMSPlant_Params.N1=N1;
HEMSPlant_Params.MaxRate_Charging=MaxRate_Charging;
HEMSPlant_Params.MaxRate_Discharging=MaxRate_Discharging;
HEMSPlant_Params.Battery_Energy_Max=Battery_Energy_Max;
HEMSPlant_Params.Battery_Energy_Min=Battery_Energy_Min;
HEMSPlant_Params.MaxRate_Discharging_StartUp=MaxRate_Discharging_StartUp;

%-------------------------- House Specifications -------------------------%

HEMSHouse_Params=[]; % Empty Struct

HEMSHouse_Params.R_w=R_w;
HEMSHouse_Params.R_attic=R_attic;
HEMSHouse_Params.R_roof=R_roof;
HEMSHouse_Params.R_im=R_im;
HEMSHouse_Params.R_win=R_win;
HEMSHouse_Params.C_w=C_w;
HEMSHouse_Params.C_attic=C_attic;
HEMSHouse_Params.C_im=C_im;
HEMSHouse_Params.C_in=C_in;
HEMSHouse_Params.C1=C1;
HEMSHouse_Params.C2=C2;
HEMSHouse_Params.C3=C3;
HEMSHouse_Params.Human_Num=Human_Num;
HEMSHouse_Params.Human_Heat=Human_Heat;
HEMSHouse_Params.Appliance_Heat=Appliance_Heat;
HEMSHouse_Params.Q_ac=AC_COP*ACLoad_Power;
HEMSHouse_Params.Cp=Cp;
HEMSHouse_Params.V=V;
HEMSHouse_Params.Den_Air=Den_Air;
HEMSHouse_Params.C_oew=C_oew;
HEMSHouse_Params.SHGC=SHGC;
HEMSHouse_Params.Alpha_w=Alpha_w;
HEMSHouse_Params.Alpha_r=Alpha_r;
HEMSHouse_Params.Area_w=Area_w;
HEMSHouse_Params.Tilt_w=Tilt_w;
HEMSHouse_Params.Azi_w=Azi_w;
HEMSHouse_Params.Area_r=Area_r;
HEMSHouse_Params.Tilt_r=Tilt_r;
HEMSHouse_Params.Azi_r=Azi_r;
HEMSHouse_Params.Area_win=Area_win;
HEMSHouse_Params.Tilt_win=Tilt_win;
HEMSHouse_Params.Azi_win=Azi_win;

HEMSHouse_Params.ACLoad_Num=ACLoad_Num;
HEMSHouse_Params.ACLoad_Power=ACLoad_Power;
HEMSHouse_Params.AC_COP=AC_COP;
HEMSHouse_Params.T_AC_Base=T_AC_Base;
HEMSHouse_Params.T_AC_DeadBand=T_AC_DeadBand;
HEMSHouse_Params.AC_Indicator=AC_Indicator;
HEMSHouse_Params.AC_Indicator1=AC_Indicator1;


end

