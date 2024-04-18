%% Author: Ninad Kiran Gaikwad
%% Date: Apr-01-2019
%% Weather Data Plotting

clear all;
clc;
close all;

%% Load File

[File,Path] = uigetfile('*.csv','Select Weather File');
FullPath = strcat(Path,File);
FullData=csvread(FullPath);

% Getting the Date-Time Matrix
DateTime_Matrix=FullData(:,1:4); % Day Month Year DecimalTime
[Row,Col]=size(DateTime_Matrix);

%% User INPUT - Date Time for Weather Data Plotting

FileRes=5; % Min - User Input

% Start Datetime
StartYear=2017; % User Defined
StartMonth=1; % User Defined
StartDay=1; % User Defined

% End Datetime
EndYear=2017; % User Defined
EndMonth=12; % User Defined
EndDay=1; % User Defined

%% Getting Appropriate Weather Data

EndTime=24-(FileRes/60);

% Creating a DateTimeMatrix_ForSlicer

DateTimeMatrixAggregate_ForSlicer=horzcat(DateTime_Matrix,zeros(Row,1));

[ ~,StartIndex_Aggregate,EndIndex_Aggregate ] = DateTimeSeriesSlicer(DateTimeMatrixAggregate_ForSlicer,1,FileRes,StartYear,EndYear,StartMonth,EndMonth,StartDay,EndDay,0,EndTime);

WeatherData=FullData(StartIndex_Aggregate:EndIndex_Aggregate,:);

DHI=WeatherData(:,5); % W/m^2
DNI=WeatherData(:,6); % W/m^2
GHI=WeatherData(:,7); % W/m^2

Ws=WeatherData(:,16); % m/s
WindDirection=WeatherData(:,18); % Degrees

T_am=WeatherData(:,20); % DegC
Pressure=WeatherData(:,21); % mbar

%% Creating DateTime Vector

for ii=1:length(GHI)
    
    Day=WeatherData(ii,1);
    Month=WeatherData(ii,2);
    Year=WeatherData(ii,3);
    Time=WeatherData(ii,4);

    % Getting Hrs Mins Secs from Decimal Time
    [ hr,minn,sec ] = DeciToHM( Time );

    % Creating DateTime Vector
    DateTimeVector(ii,1)=datetime(Year,Month,Day,hr,minn,sec);    
    
end

%% Plots

figure;
hold on
plot(DateTimeVector,GHI,'-b','LineWidth',2,'DatetimeTickFormat','HH:mm');
title('GHI');
xlabel('Time (HH:mm)');
ylabel('Power (W/m^2)');
%legend('GHI');
hold off;

figure;
hold on
plot(DateTimeVector,T_am,'-b','LineWidth',2,'DatetimeTickFormat','HH:mm');
title('Temperature');
xlabel('Time (HH:mm)');
ylabel('Temperature (DegC)');
%legend('Temperature');
hold off;

figure;
hold on
plot(DateTimeVector,Ws,'-b','LineWidth',2,'DatetimeTickFormat','HH:mm');
title('Wind Speed');
xlabel('Time (HH:mm)');
ylabel('Wind Speed (m/s)');
%legend('Wind Speed');
hold off;
