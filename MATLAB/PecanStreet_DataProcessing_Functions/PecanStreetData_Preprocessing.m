%% Pecan Street Data Preprocessing : Converting Many Files to 1

% Author: Ninad Kiran Gaikwad
% Date: Mar/15/2021
% Description: Pecan Street Data Preprocessing 


%%

clear all;
clc;

%% Step 1: Reading CSV File

DataFolderPath='C:\Users\ninad\Dropbox (UFL)\NinadGaikwad_PhD\Gaikwad_Research\Gaikwad_Research_Work\20_Gaikwad_SmartCommunity\data\RawFiles\15minute_data_newyork\';

FileName='15minute_data_newyork.csv';

OriginalResolution = 15; % Mins

NewResolution = 10; % Mins

ProcessingType = 1; % 1 - Full File to be Processed ; 2 - Part of the File to be Processed User Input is Required

AveragingPoints = 2; %

DataFullPath=strcat(DataFolderPath,FileName);

ActualFile= readtable(DataFullPath);

% Getting relevant columns of data in Array format
ID_Array=table2array(ActualFile(:,1));
%ID_Array=table2array(ActualFile(1:10000,1)); % Debugging

DateTime_String_Array=table2array(ActualFile(:,2));
%DateTime_String_Array=table2array(ActualFile(1:10000,2)); % Debugging

Data_Array=table2array(ActualFile(:,3:end));

[DataFrame_Row,DataFrame_Column]=size(Data_Array);

% Getting DataCols
DataCols=DataFrame_Column-2;
           
% Deleting ActualFile variable
clear ActualFile;

%% Step 2: Changing Date-Time Stamp Columns for Utility

DateTimeStamp=zeros(1,4); % Initialization

for i=1:length(DateTime_String_Array) % For each row in ActualFile
    
   DateTime_String= DateTime_String_Array(i);
    
   DateTime_String1 = split(datestr(DateTime_String,'yyyy-mm-dd HH:MM'),' ');
   
   DateTime_String_Date = split(DateTime_String1(1),'-');
   
   DateTime_String_Time = split(DateTime_String1(2),':');
   
   Year = str2num(DateTime_String_Date{1});
   Month = str2num(DateTime_String_Date{2});
   Day = str2num(DateTime_String_Date{3});
   
   Hour = str2num(DateTime_String_Time{1});
   Min = str2num(DateTime_String_Time{2});
   
   [ TimeDeci ] = HMToDeci( Hour,Min,0 );
   
   DateTimeStamp(i,1:4)=[Day,Month,Year,TimeDeci];
    
   i % Debugger
   
end

%% Step 3: Grouping Data according to different houses

Unique_Houses=unique(ID_Array);

File_Num=length(Unique_Houses) % Debugger

FileNum_Current=0; % Debugger

for i=1:length(Unique_Houses) % For each House
    
    % Incrementing FileNum_Current
    FileNum_Current=FileNum_Current+1
    
    % Finding Indices for the current House
    CurrentHouse_Indices=find(ID_Array==Unique_Houses(i));
    
    % Creating current House Dataframe
    CurrentHouse_Dataframe=[DateTimeStamp(CurrentHouse_Indices,:),Data_Array(CurrentHouse_Indices,:)];

    % Converting NaNs to Zeros
    CurrentHouse_Dataframe(isnan(CurrentHouse_Dataframe))=0;
    
    % Clean Data for any irregularities of date-time ordering and missing data (Negatives are not converted to 0s)  
    [ ProcessedDataFrame ] = SolarPVWeatherDataCleaner_ModifiedForPecanStreet( OriginalResolution,DataCols,AveragingPoints,CurrentHouse_Dataframe );
 
    % Changing to required Resolution
    [CurrentHouse_Dataframe_NewResFile] = PecanStreet_Low2HighRes(OriginalResolution,NewResolution,ProcessingType,ProcessedDataFrame);
    
    % Current House Name
    CurrentHouse_Name=strcat('House_',num2str(i),'_',num2str(OriginalResolution),...
        'minTo_',num2str(NewResolution),'min','.csv');
    
    % Saving House Data in a CSV File
    csvwrite(strcat(DataFolderPath,CurrentHouse_Name),CurrentHouse_Dataframe_NewResFile);
    
end


