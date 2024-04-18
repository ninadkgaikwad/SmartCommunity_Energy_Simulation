%% Weather Data Analysis

clear all;
clc;

%% User Input - Aggregation Information

FileRes=5; % User Input in Minutes

EndTime=24-(FileRes/60);

SolarIrradiance_Columns=[1 2 3 4 5 6]; % User Input [DHI,DNI,GHI,Clear Sky DHI,Clear Sky DNI,Clear Sky GHI]

Other_WeatherData_Columns=[16 12  13 ]; % User Input [Temperature, Windspeed, Precipitable Water]

All_WeatherData_Columns=4+[SolarIrradiance_Columns, Other_WeatherData_Columns ];

Aggregation_Days= 7;% User Input

Aggregation_Hours=(Aggregation_Days-1)*24; % Number of Hours

Aggregation_Duration=duration(Aggregation_Hours,0,0); % Light CeilingFan Refrigerator Total

% Aggregation Start Datetime

StartYear=2017; % User Defined
StartMonth=1; % User Defined
StartDay=1; % User Defined

% Aggregation End Datetime
EndYear=2017; % User Defined
EndDay=30; % User Defined
EndMonth=12; % User Defined

%% Getting the Pre-Processed WeatherFile File

Weather_FileName='GarzasAdjuntas_2017_To_2017_WeatherData_NSRDB_30minTo5minRes.csv'; % Enter File Name

WeatherFile_FolderName='Data'; % Constant Folder

WeatherFile_FullPath=strcat(WeatherFile_FolderName,'/',Weather_FileName);

WeatherFileLoad_Full=csvread(WeatherFile_FullPath); % WeatherFileLoad_Full

[R,C] = size(WeatherFileLoad_Full);

% Getting Approprite Columns from the Weather Data file

WeatherFile_ReqCols=WeatherFileLoad_Full(:,[1:4, All_WeatherData_Columns]);

SolarData=WeatherFile_ReqCols(:,5:4+length(SolarIrradiance_Columns));

OtherWeatherData=WeatherFile_ReqCols(:,5+length(SolarIrradiance_Columns):4+length(SolarIrradiance_Columns)+length(Other_WeatherData_Columns));

[ ~,AggregateStartIndex,AggregateEndIndex ] = DateTimeSeriesSlicer(WeatherFile_ReqCols,1,FileRes,StartYear,EndYear,StartMonth,EndMonth,StartDay,EndDay,0,EndTime);

%% Aggregating the Weather Variables

WeatherFile_CurrentStarting_RowNum=AggregateStartIndex; % Initializing

Counter=0; % Initializing Counter

WeatherFile_Aggregate=zeros(1,4+4+length(SolarIrradiance_Columns)+length(Other_WeatherData_Columns)); % StartDateTime EndDateTime Light CeilingFan Refrigerator Total 

EndDateTime_HouseLoad_Full=datetime(WeatherFile_ReqCols(AggregateEndIndex,3),WeatherFile_ReqCols(AggregateEndIndex,2),WeatherFile_ReqCols(AggregateEndIndex,1));

while (~isnan(WeatherFile_CurrentStarting_RowNum))
    
    if (WeatherFile_CurrentStarting_RowNum>R)
       
        break;
        
    end
    
    % Getting the Current Date 
    Day_Current=WeatherFile_ReqCols(WeatherFile_CurrentStarting_RowNum,1);
    Month_Current=WeatherFile_ReqCols(WeatherFile_CurrentStarting_RowNum,2);
    Year_Current=WeatherFile_ReqCols(WeatherFile_CurrentStarting_RowNum,3);
    
    % Creating a Datetime object    
    CurrentDate=datetime(Year_Current,Month_Current,Day_Current,0,0,0);
    
    % Computing Next Date
    NextDate=CurrentDate+Aggregation_Duration;
    
    if (NextDate<=EndDateTime_HouseLoad_Full)
        
        % Incrementing Counter
        Counter=Counter+1 % Debugger
    
        % Getting Next Date
        Day_Next=NextDate.Day;
        Month_Next=NextDate.Month;
        Year_Next=NextDate.Year;

        % Getting the Appropriate Indices     
        [ ~,StartIndex,EndIndex ] = DateTimeSeriesSlicer(WeatherFile_ReqCols,1,FileRes,Year_Current,Year_Next,Month_Current,Month_Next,Day_Current,Day_Next,0,EndTime);

        % Computing Insolation Values (kWh)
        Insolation_RowMatrix=zeros(1,length(SolarIrradiance_Columns)); % Initialization
        
        Current_SolarData_ColumnMatrix=SolarData(StartIndex:EndIndex,:);
        
        for jj=1:length(SolarIrradiance_Columns) % For Each Solar Irradiance Column
            
            Current_SolarData_SingleColumn=Current_SolarData_ColumnMatrix(:,jj);
            
            Current_SolarData_SingleColumn_NonZeroIndices=find(Current_SolarData_SingleColumn);
            
            Current_SolarData_SingleColumn_NonZero=Current_SolarData_SingleColumn(Current_SolarData_SingleColumn_NonZeroIndices);
            
            Current_Insolation=sum(Current_SolarData_SingleColumn_NonZero)*(FileRes/60)*(1/1000);
            
            Insolation_RowMatrix(1,jj)=Current_Insolation;
            
        end
        
        % Computing the Average of Other Weather Variables        
        Current_OtherWeather_ColumnMatrix=OtherWeatherData(StartIndex:EndIndex,:);
        
        [R_OW,C_OW]=size(Current_OtherWeather_ColumnMatrix);
        
        Current_OtherWeatherRowMatrix=sum(Current_OtherWeather_ColumnMatrix)/R_OW;        
                
        % Appropriately Filling HouseLoad_Aggregate
        
        WeatherFile_Aggregate(Counter,1:end)=[Day_Current,Month_Current,Year_Current,0,Day_Next,Month_Next,Year_Next ,EndTime, Insolation_RowMatrix ,Current_OtherWeatherRowMatrix ];
        
        % Reinitializing HouseLoad_CurrentStarting_RowNum
        
        WeatherFile_CurrentStarting_RowNum=EndIndex+1;
        
    else
        
        % While Loop Break Condition
        
        WeatherFile_CurrentStarting_RowNum=nan;
        
    end

    
end

%% Writing the Weather Data Matrix in CSV File

Weather_FileName='Adjuntas_2017-2017_WeatherAggregate_7Days.csv'; % Enter File Name

WeatherFile_FolderName='WeatherData_AggregateCSVFiles'; % Constant Folder

WeatherFile_FullPath=strcat(WeatherFile_FolderName,'/',Weather_FileName);

csvwrite(WeatherFile_FullPath,WeatherFile_Aggregate);
