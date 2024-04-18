function [HEMSWeatherData_Output] = WeatherData_Extractor(HEMSWeatherData_Input,Simulation_Params,WeatherData_FileName)

% Author: Ninad Kiran Gaikwad
% Date: Jun/10/2019
% Description: WEATHERDATA_EXTRACTOR - Extracts Required Weather Data from Pre-Processed Weather CSV File

%% WEATHERDATA_EXTRACTOR - Extracts Required Weather Data from Pre-Processed Weather CSV File

%% Getting Required data from HEMSWeatherData_Input - Struct

%---------------------- HEMSWeatherData_Input ----------------------------%
WeatherDataFile_Path=HEMSWeatherData_Input.WeatherDataFile_Path;
StartYear=HEMSWeatherData_Input.StartYear;
StartMonth=HEMSWeatherData_Input.StartMonth;
StartDay=HEMSWeatherData_Input.StartDay;
StartTime=HEMSWeatherData_Input.StartTime;
EndYear=HEMSWeatherData_Input.EndYear;
EndMonth=HEMSWeatherData_Input.EndMonth;
EndDay=HEMSWeatherData_Input.EndDay;
EndTime=HEMSWeatherData_Input.EndTime;

%------------------------- Simulation_Params------------------------------%
FileRes=Simulation_Params.FileRes;
Simulation_StepSize=Simulation_Params.Simulation_StepSize;
StepSize=Simulation_Params.StepSize;

%% Getting DateTime and Weather Data

FullData=csvread(WeatherDataFile_Path);

% Getting the Date-Time Matrix
DateTime_Matrix=FullData(:,1:4);    % Day Month Year DecimalTime
[Row,Col]=size(DateTime_Matrix);

%% Computing Correct Indices for Weather Data Extraction

% Computing End Time
% EndTime=24-(FileRes/60);

% Creating a DateTimeMatrix_ForSlicer
DateTimeMatrixAggregate_ForSlicer=horzcat(DateTime_Matrix,zeros(Row,1));

% Computing required Start and End Indices
[ ~,StartIndex_Aggregate,EndIndex_Aggregate ] = DateTimeSeriesSlicer(DateTimeMatrixAggregate_ForSlicer,1,FileRes,StartYear,EndYear,StartMonth,EndMonth,StartDay,EndDay,StartTime,EndTime);

%% Getting required Weather Data

% Getting Wind Speed Data
Ws=FullData(StartIndex_Aggregate:EndIndex_Aggregate,16);         % m/s

% Getting Temperature Data
T_am=FullData(StartIndex_Aggregate:EndIndex_Aggregate,20);       % Deg C

% Getting Irradiance Data
GHI=FullData(StartIndex_Aggregate:EndIndex_Aggregate,7);         % W/m^(2)
DNI=FullData(StartIndex_Aggregate:EndIndex_Aggregate,6);         % W/m^2

%% Creating DateTime Vector
Counter=0; % Initializing Counter

for ii=StartIndex_Aggregate:EndIndex_Aggregate
    
    % Incrementing Counter
    Counter=Counter+1;
    
   Day=DateTime_Matrix(ii,1);
   Month=DateTime_Matrix(ii,2);
   Year=DateTime_Matrix(ii,3);
   Time=DateTime_Matrix(ii,4);

   % Getting Hrs Mins Secs from Decimal Time
   [ hr,minn,sec ] = DeciToHM( Time );

    % Creating DateTime Vector
    DateTimeVector(Counter,1)=datetime(Year,Month,Day,hr,minn,sec);   

end

% Getting the DateTime Matrix
DateTime_Matrix=FullData(StartIndex_Aggregate:EndIndex_Aggregate,1:4); % Day Month Year DecimalTime

%% Creating HEMSWeatherData_Output - Struct

HEMSWeatherData_Output=[];

HEMSWeatherData_Output.Ws=Ws;
HEMSWeatherData_Output.T_am=T_am;
HEMSWeatherData_Output.GHI=GHI;
HEMSWeatherData_Output.DNI=DNI;
HEMSWeatherData_Output.DateTimeVector=DateTimeVector;
HEMSWeatherData_Output.DateTime_Matrix=DateTime_Matrix;

save(WeatherData_FileName,'HEMSWeatherData_Output');
%% Extra Code for debugging - Creating One Day Uniform Data

% Getting Full Weather Dataset for the desired Period
Full_Data=FullData(StartIndex_Aggregate:EndIndex_Aggregate,:);

% Creating Uniform single day Data Weatherfile

for ii=1:360
    if (ii==1)
        SingleDay_WeatherFile_1=Full_Data;
    end
    SingleDay_WeatherFile_1=vertcat(SingleDay_WeatherFile_1,Full_Data);
end

SingleDay_WeatherFile=SingleDay_WeatherFile_1;
