%% Data Preprocessing : Converting Many Files to 1

% Weather Data Source   : NSRDB
% Location              : Garzas, Adjuntas, Puerto Rico
% Latitude              : 18.17 Deg
% Longitude             : -66.74 Deg
% Time Zone             : 4 Hours

%%

clear all;
clc;

%% Step 1: Reading All separate CSV Files as One

DataFolderName='Data';

DataFolderStructure=dir(DataFolderName);

FileNames_Cell=extractfield(DataFolderStructure,'name');

FolderPath=DataFolderStructure.folder;

Counter=0; % Initialization

for i=1:length(FileNames_Cell) % For each file in the Folder
    
    FileName=FileNames_Cell{i};
    
    startIndex = regexp(FileName,'\<[^a-zA-Z_0-9]\w*'); % To check for temporary files
    
   if (length(startIndex)==0) % If actual File
       
       % Incrementing Counter
       
       Counter=Counter+1
       
       FullPathName=strcat(FolderPath,'/',FileName); % Creating Full Path for the Actual File
       
       if (Counter==1)
           
           ActualFile= csvread(FullPathName,3);
           
       else
       
           ActualFile_1= csvread(FullPathName,3);

           ActualFile= vertcat(ActualFile,ActualFile_1);
       
       end
            
       
   end
    
end

%% Step 2: Changing Date-Time Stamp Columns for Utility

DateTimeStamp=zeros(1,4); % Initialization

for i=1:length(ActualFile) % For each row in ActualFile
    
   Hour=ActualFile(i,4);
   
   Min=ActualFile(i,5);
   
   [ TimeDeci ] = HMToDeci( Hour,Min,0 );
   
   DateTimeStamp(i,1:4)=[ActualFile(i,3),ActualFile(i,2),ActualFile(i,1),TimeDeci];
    
   i % Debugger
   
end

ActualFile=ActualFile(:,6:end); % Removing older Date-Time Stamp

ActualFile=horzcat(DateTimeStamp,ActualFile); % Adding New Date-Time Stamp

%% Step 3: Write to a CSV File

CombinedFileName='GarzasAdjuntas_1998_To_2017_WeatherData_NSRDB.csv';

csvwrite(strcat(DataFolderName,'/',CombinedFileName),ActualFile);


