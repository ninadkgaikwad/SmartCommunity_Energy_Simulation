%% Main File - Baseline Controller for Community of Houses [HEMS]

% Author: Ninad Kiran Gaikwad
% Date: Feb/1/2021
% Description: Main Baseline Controller for Community of Houses

clear all;
clc;
close all;

%-------------------------------- Add Paths ------------------------------%

addpath(pwd,'../Controller_Functions');
addpath(pwd,'../NSRDB_DataProcessingFunctions');
addpath(pwd,'../PecanStreet_DataProcessing_Functions');
addpath(pwd,'../Plant_Supporting_Functions');
addpath(pwd,'../PlantFunctions');
addpath(pwd,'../CodeFromSWEEFA');

%% Simulation - User Inputs

%----------------------------- OS Information ----------------------------%

OS=2; % 1 - Linux ; 2 - Windows-Laptop ; 3 - Windows-PC

%-------------------------- Simulation Step Sizes ------------------------%

FileRes=10;                         % in Minutes                                                                                                                                                                
Simulation_StepSize = FileRes/60;   % in Hours
StepSize = FileRes*60;              % in Seconds
SmartCommunity_ControllerType=1;    % 1 = Smart Local Controller ; 2 = Dumb Local Controller

%-------------------------- Simulation Parameters ------------------------%

SimulationType = 0; % Important for Single Large House Simulation [1,2,3,4,5,6,7,8] == [N_PV_Bat_EV, N_PV_Bat, N_PV_EV, N_Bat_EV, N_PV, N_Bat, N_EV, N_None]

LoadDataType=2; % 1 - File Generated from Preprocessed Pecan Street Data files ; 2 - .mat File already exists

WeatherDataType=2; % 1 - File Generated from Preprocessed NSRDB File ; 2 - .mat File already exists

Single_House_Plotting_Index=1; % House Index for Single House Plotting

%------------------------- Community Specification -----------------------%

N_PV_Bat=1;     % Houses with both PV and Battery
N_PV=1;         % Houses with just PV
N_Bat=1;        % Houses with just Battery
N_None=1;       % Houses with niether PV and Battery

% Computing Total Number of Houses
N_House=N_PV_Bat+N_PV+N_Bat+N_None;

N_House_Vector=[N_PV_Bat,N_Bat,N_PV,N_None];

%----------------- Plant Initial Condition Specification -----------------%

% House Temperature Intial Condition
T_AC_Base=24;
T_House_Variance=0.5;

% Battery Initial Condition
N1=1;                           % User Input - Battery Max Changing Factor
Battery_Energy_Max = 13.5*N1;

%--------------------- Simulation Period Specification -------------------%

% Load Computation Start Date
StartYear=2017;     % User Defined
StartMonth=9;      % User Defined
StartDay=11;         % User Defined
StartTime=0;        % User Defined

% Load Computation End Date
EndYear=2017;       % User Defined
EndMonth=9;        % User Defined
EndDay=18;          % User Defined
EndTime=24-(FileRes/60);          %24-(FileRes/60);  

%----------------------- Folder Paths Specification ----------------------%

ImageFolder_Name='Gainesville_BaseLine_7DayTest_SC_PVBat1_Bat1_PV1_None1_SCL1_';

SimulationData_FileName='FigurePlotterData_Gainesville_BaseLine_7DayTest_SC_PVBat1_Bat1_PV1_None1_SLC1';

SimulationPerformanceData_FileName='PerformanceData_Gainesville_BaseLine_7DayTest_SC_PVBat1_Bat1_PV1_None1_SLC1';

LoadData_FileName='PecanStreet_LoadData_SC_PVBat1_Bat1_PV1_None1';

WeatherData_FileName='Gainesville_Irma';

%-------------------- Weather Data Location and Period -------------------%

% Getting to Weather Data Folder in the Correct OS Folder
if (OS==1) % Linux
    
    WeatherDataFile_Path='/home/ninadgaikwad/Dropbox (UFL)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/19_Resiliency/codes/Matlab_Scripts_New/CCTA_2020/Laptop_Final_Improved/DwellTime_CNCL_WithoutL1/Data/Gainesville_2017_To_2017_WeatherData_NSRDB_30minTo10minRes.csv';
    
    LoadDataFolder_Path='/home/ninadgaikwad/Dropbox (UFL)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/20_Gaikwad_SmartCommunity/data/PreProcessedFiles/10minute_data_austin_HouseWise/';
    
elseif (OS==2) % Windows-Laptop
    
     WeatherDataFile_Path='C:\Users\ninad\Dropbox (UFL)\NinadGaikwad_PhD\Gaikwad_Research\Gaikwad_Research_Work\19_Resiliency\codes\Matlab_Scripts_New\CCTA_2020\Laptop_Final_Improved\DwellTime_CNCL_WithoutL1\Data\Gainesville_2017_To_2017_WeatherData_NSRDB_30minTo10minRes.csv';
    
     LoadDataFolder_Path='C:\Users\ninad\Dropbox (UFL)\NinadGaikwad_PhD\Gaikwad_Research\Gaikwad_Research_Work\20_Gaikwad_SmartCommunity\data\PreProcessedFiles\10minute_data_austin_HouseWise\';
          
elseif (OS==3) % Windows-PC
    
     WeatherDataFile_Path='C:\Users\Me!\Dropbox (UFL)\NinadGaikwad_PhD\Gaikwad_Research\Gaikwad_Research_Work\19_Resiliency\codes\Matlab_Scripts_New\CCTA_2020\Laptop_Final_Improved\DwellTime_CNCL_WithoutL1\Data\Gainesville_2017_To_2017_WeatherData_NSRDB_30minTo10minRes.csv';
    
     LoadDataFolder_Path='C:\Users\Me!\Dropbox (UFL)\NinadGaikwad_PhD\Gaikwad_Research\Gaikwad_Research_Work\20_Gaikwad_SmartCommunity\data\PreProcessedFiles\10minute_data_austin_HouseWise\';
          
end



%% Weather Data Extraction

% Creating Simulation_Params Struct

Simulation_Params=[]; % Empty Struct

Simulation_Params.FileRes=FileRes;
Simulation_Params.Simulation_StepSize=Simulation_StepSize;
Simulation_Params.StepSize=StepSize;
Simulation_Params.SmartCommunity_ControllerType=SmartCommunity_ControllerType;

% Creating HEMSWeatherData_Input Struct

HEMSWeatherData_Input=[]; % Empty Struct

HEMSWeatherData_Input.WeatherDataFile_Path=WeatherDataFile_Path;
HEMSWeatherData_Input.StartYear=StartYear;
HEMSWeatherData_Input.StartMonth=StartMonth;
HEMSWeatherData_Input.StartDay=StartDay;
HEMSWeatherData_Input.StartTime=StartTime;
HEMSWeatherData_Input.EndYear=EndYear;
HEMSWeatherData_Input.EndMonth=EndMonth;
HEMSWeatherData_Input.EndDay=EndDay;
HEMSWeatherData_Input.EndTime=EndTime;

if (WeatherDataType==1) % We do not have Weather Data File

    [HEMSWeatherData_Output] = WeatherData_Extractor(HEMSWeatherData_Input,Simulation_Params,WeatherData_FileName);

elseif (WeatherDataType==2) % We have Weather Data File
    
    load(strcat(WeatherData_FileName,'.mat'))
    
end

%% Load Data Extraction

%Type=1; % Type of Load Data Extraction
if (LoadDataType==1) % We do not have Load Data File
    
    [PecanStreet_Data_Output] = PecanStreet_Data_Extractor(HEMSWeatherData_Input,Simulation_Params,LoadDataFolder_Path,N_House_Vector,SimulationType,LoadData_FileName);

elseif (LoadDataType==2) % We already have Load Data File
    
    load(strcat(LoadData_FileName,'.mat'));
    
end

%% Basic Computation

%-------------------- Creating Community_Params Struct -------------------%

Community_Params.N_House=N_House;
Community_Params.N_PV_Bat=N_PV_Bat;
Community_Params.N_Bat=N_Bat;
Community_Params.N_PV=N_PV;
Community_Params.N_None=N_None;

%------------------- From Extracted Weather Data -------------------------%

Ws=HEMSWeatherData_Output.Ws;
T_am=HEMSWeatherData_Output.T_am;
GHI=HEMSWeatherData_Output.GHI;
DNI=HEMSWeatherData_Output.DNI;
DateTimeVector=HEMSWeatherData_Output.DateTimeVector;
DateTime_Matrix=HEMSWeatherData_Output.DateTime_Matrix;

Simulation_Steps_Total=length(DateTimeVector);

%------------------------ From Extracted Load Data -----------------------%

% Getting Renewable Source Data
SolarGen_Data=PecanStreet_Data_Output(:,4+1:6,:);

Battery_ChargerDischarge_Data=PecanStreet_Data_Output(:,7,:);

EVCharging_Data=PecanStreet_Data_Output(:,7+1:9,:);

E_LoadData=PecanStreet_Data_Output(:,:,:);

% Making Negatives (-) = 0 in LoadData
E_LoadData(E_LoadData(:,9+1:end,:)<0)=0;

% Creating 8 Level Priority Load Data
E_Load_P1=sum(E_LoadData(:,9+1:9+12,:),2); % Priority Level1 Sum
E_Load_P2=sum(E_LoadData(:,21+1:21+5,:),2); % Priority Level2 Sum
E_Load_P3=sum(E_LoadData(:,26+1:26+3,:),2); % Priority Level3 Sum
E_Load_P4=sum(E_LoadData(:,29+1:29+7,:),2); % Priority Level4 Sum
E_Load_P5=sum(E_LoadData(:,36+1:36+6,:),2); % Priority Level5 Sum
E_Load_P6=sum(E_LoadData(:,42+1:42+6,:),2); % Priority Level6 Sum
E_Load_P7=sum(E_LoadData(:,48+1:48+3,:),2); % Priority Level7 Sum
E_Load_P8=sum(E_LoadData(:,51+1:51+6,:),2); % Priority Level8 Sum

E_LoadData=[E_LoadData(:,1:9,:)...
            E_Load_P1 E_Load_P2 E_Load_P3 E_Load_P4...
            E_Load_P5 E_Load_P6 E_Load_P7 E_Load_P8];

% Creating E_Load_Desired from LoadData (Summing all loads for a given house)
E_Load_Desired_Array=sum(E_LoadData(:,9+1:end,:),2);

% Creating E_Load_Desired from E_Load_Array
for ii=1:N_House
    E_Load_Desired(:,ii)=E_Load_Desired_Array(:,:,ii);
end

%% Community-House Parameter Generation

[HEMSPlant_Params,HEMSHouse_Params] = HEMS_CommunityHouse_Parameter_Generator(Community_Params,Simulation_Params);

%% Initial States,Disturbance and Control Struct Creation

%--------------------------- Initial Conditions --------------------------%

rng(1); % Setting Randomness Seed for repeatability

T_House_Initial=T_AC_Base*ones(1,N_House)+(T_House_Variance*randn(1,N_House));

E_Bat_Initial=[Battery_Energy_Max*ones(1,N_PV_Bat+N_Bat),zeros(1,length(1:N_House-(N_PV_Bat+N_Bat)))];

%------------------------ Current State X_k_Plant ------------------------%

X_k_Plant=zeros(2,38,N_House); % Initializing Size

for ii=1:N_House
    
    % House Temperatures
    X_k_Plant(1,7:10,ii)=[T_House_Initial(ii), T_House_Initial(ii), T_House_Initial(ii), T_House_Initial(ii)];
    
    % Battery Initial Conditions
    if(ii<=N_PV_Bat+N_Bat)
        X_k_Plant(1,4,ii)=[E_Bat_Initial(ii)];        
    end
    
    % AC on-off Status Previous
    X_k_Plant(1,30,ii)=[1];
    
    % Prioritized Loads on-off Status Previous
    X_k_Plant(1,31:end,ii)=ones(1,8,1);    
    
end

%---------------------- Current Disturbance W_k_Plant --------------------%

% Weather 
Weather_k_Plant=[];

Weather_k_Plant.Ws=Ws(1);
Weather_k_Plant.T_am=T_am(1);
Weather_k_Plant.GHI=GHI(1);
Weather_k_Plant.DNI=DNI(1);
Weather_k_Plant.DateTime_Matrix=DateTime_Matrix(1,1:4);

% Load Data
LoadData_k_Plant=[];

LoadData_k_Plant.E_Load_Desired=E_Load_Desired(1,:);
LoadData_k_Plant.E_LoadData=E_LoadData(1,:,:);

% W_k_Plant
W_k_Plant=[];

W_k_Plant.Weather_k_Plant=Weather_k_Plant;
W_k_Plant.LoadData_k_Plant=LoadData_k_Plant;

%------------------ Initialializing Plant State History ------------------%
X_k_Plant_History=X_k_Plant;

%------------------ Initialializing Controller History ------------------%
U_k_History=zeros(1,11,N_House);

%% HEMS Plant Simulation 

T_Start = tic; % Measuring Time of Simulation

for ii=1:Simulation_Steps_Total % For each Simulation Time Step
    
    ii % For knowing what iteration is going on
    
    % Step 1: Compute Control Command
    
    tic; % Measuring Time for Controller
    
    if (SmartCommunity_ControllerType==1) % Smart Local Controller
        
        [U_k] = HEMS_Smart_LocalController(X_k_Plant,W_k_Plant,HEMSPlant_Params,Community_Params,Simulation_Params);
        
    elseif (SmartCommunity_ControllerType==2) % Dumb Local Controller
        
        [U_k] = HEMS_Dumb_LocalController(X_k_Plant,W_k_Plant,HEMSPlant_Params,Community_Params,Simulation_Params);
        
    end
    
    TimePer_Controller(ii)=toc; % Measuring Time for Controller
    
    % Step 2: Compute Next Plant State
    
    [X_k_Plus_Plant] = HEMS_Plant(X_k_Plant,W_k_Plant,U_k,HEMSPlant_Params,HEMSHouse_Params,Community_Params,Simulation_Params);
    
    % Step 3: Update Plant History
    
    X_k_Plant_History(ii:ii+1,:,:)=X_k_Plus_Plant;
    
    % Step 4: Update Current Plant State
    
    X_k_Plant(1,:,:)=X_k_Plus_Plant(2,:,:);
    
    % Step 5: Update Current Disturbance
    
    if (ii<Simulation_Steps_Total)
    
        % Weather 
        Weather_k_Plant=[];

        Weather_k_Plant.Ws=Ws(ii+1);
        Weather_k_Plant.T_am=T_am(ii+1);
        Weather_k_Plant.GHI=GHI(ii+1);
        Weather_k_Plant.DNI=DNI(ii+1);
        Weather_k_Plant.DateTime_Matrix=DateTime_Matrix(ii+1,:);

        % Load Data
        LoadData_k_Plant=[];

        LoadData_k_Plant.E_Load_Desired=E_Load_Desired(ii+1,:);
        LoadData_k_Plant.E_LoadData=E_LoadData(ii+1,:,:);

        % W_k_Plant
        W_k_Plant=[];

        W_k_Plant.Weather_k_Plant=Weather_k_Plant;
        W_k_Plant.LoadData_k_Plant=LoadData_k_Plant; 
        
    else
        
        % Do Nothing as no update required
        
    end    
    
    % Step 6: Update Controller History
    U_k_History(ii,:,:)=U_k;    
    
end

% Computing Time for completion of Simulation
TimeCompletion_Simulation=toc(T_Start); % Measuring Time for Simulation 
Avg_TimePer_Controller=mean(TimePer_Controller);

fprintf('\n Time to complete simulation = %f',TimeCompletion_Simulation)
fprintf('\n Average Time per Controller = %f',Avg_TimePer_Controller)

%% Creating Plots

% Folder Paths for saving Variables
Baseline_Output_Images_Path=strcat(pwd,'/Baseline_Output_Images','/',ImageFolder_Name);

% Creating the HEMS_Plant_Baseline_FigurePlotter_Input - Struct
HEMS_Plant_FigurePlotter_Input = []; % Empty Structs

HEMS_Plant_FigurePlotter_Input.X_k_Plant_History=X_k_Plant_History;
HEMS_Plant_FigurePlotter_Input.U_k_History=U_k_History; 
HEMS_Plant_FigurePlotter_Input.E_LoadData=E_LoadData(:,9+1:end,:)/HEMSPlant_Params.Eff_Inv;
HEMS_Plant_FigurePlotter_Input.E_Load_Desired=E_Load_Desired_Array/HEMSPlant_Params.Eff_Inv;
HEMS_Plant_FigurePlotter_Input.HEMSWeatherData_Output=HEMSWeatherData_Output;
HEMS_Plant_FigurePlotter_Input.HEMSPlant_Params=HEMSPlant_Params;
HEMS_Plant_FigurePlotter_Input.Community_Params=Community_Params;
HEMS_Plant_FigurePlotter_Input.Baseline_Output_Images_Path=Baseline_Output_Images_Path;
HEMS_Plant_FigurePlotter_Input.Single_House_Plotting_Index=Single_House_Plotting_Index;
HEMS_Plant_FigurePlotter_Input.Simulation_Params=Simulation_Params;

HEMS_Plant_FigurePlotter_Input.TimeCompletion_Simulation=TimeCompletion_Simulation;
HEMS_Plant_FigurePlotter_Input.Avg_TimePer_Controller=Avg_TimePer_Controller;
HEMS_Plant_FigurePlotter_Input.TimePer_Controller=TimePer_Controller;

% Saving FigurePlotterData
save(SimulationData_FileName,'HEMS_Plant_FigurePlotter_Input');

% Plotting using External Function
HEMS_Plant_FigurePlotter(HEMS_Plant_FigurePlotter_Input);

%% Performance Computation

% Creating the HEMS_Plant_Baseline_FigurePlotter_Input - Struct
HEMS_PerformanceComputation = []; % Empty Structs

HEMS_PerformanceComputation.X_k_Plant_History=X_k_Plant_History;
HEMS_PerformanceComputation.U_k_History=U_k_History; 
HEMS_PerformanceComputation.E_LoadData=E_LoadData(:,9+1:end,:)/HEMSPlant_Params.Eff_Inv;
HEMS_PerformanceComputation.E_Load_Desired=E_Load_Desired_Array/HEMSPlant_Params.Eff_Inv;
HEMS_PerformanceComputation.HEMSWeatherData_Output=HEMSWeatherData_Output;
HEMS_PerformanceComputation.HEMSPlant_Params=HEMSPlant_Params;
HEMS_PerformanceComputation.Community_Params=Community_Params;

% using External Function
[Plant_Performance]=HEMS_Plant_Performance_Computer(HEMS_PerformanceComputation);

save(SimulationPerformanceData_FileName,'Plant_Performance');