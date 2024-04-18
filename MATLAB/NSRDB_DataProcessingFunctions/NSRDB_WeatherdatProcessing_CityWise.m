%% NSRDB - Weather Data Processing Script - City Wise

% Author: Ninad Kiran Gaikwad
% Date: Jan/31/2020
% Description: NSRDB - Weather Data Processing Script - City Wise

%% Step 1 - Getting the Raw Weather Data Folder

MainFolder_Path='/home/ninadgaikwad/Dropbox (UFL)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/NSRDB_DATA_PROCESSING/NSRDB_USA_WeatherData_ByCity_New';

DestinationFolder_Path='/home/ninadgaikwad/Dropbox (UFL)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/NSRDB_DATA_PROCESSING/NSRDB_USA_Database_1998_2018_ByCity_DiffRes';

 
%% Step 2 - Getting List of all Sub-Folders

DataFolderStructure1=dir(MainFolder_Path);

FileNames_Cell1=extractfield(DataFolderStructure1,'name');

FolderPath1=DataFolderStructure1.folder;

% Initialization
FolderCounter=0;

%% Step 3 - Starting the Main Folder Loop for Weather File Processing

for i=1:length(FileNames_Cell1) % For each folder in the Folder
    
    FileName1=FileNames_Cell1{i};
    
    startIndex1 = regexp(FileName1,'\<[^a-zA-Z_0-9]\w*'); % To check for temporary files
    
   if (length(startIndex1)==0) % If actual File
       
        % Incrementing FileCounter
        FolderCounter=FolderCounter+1         
       
       % Getting CityStateName
       CityStateName=FileName1;
       
       % Subfolder Path
       FolderPath2=[FolderPath1,'/',FileName1];
       
       % Getting Files inside the CityState Folder
        DataFolderStructure2=dir(FolderPath2);

        FileNames_Cell2=extractfield(DataFolderStructure2,'name');

        FolderPath2=DataFolderStructure2.folder; 
        
        % Initialization
        FileCounter=0;
        
        %% Step 4 - Starting the Main Subfolder Loop for Weather File Processing
        
        for i=1:length(FileNames_Cell2) % For each file in the Subfolder
            
            FileName3=FileNames_Cell2{i};

            startIndex2 = regexp(FileName3,'\<[^a-zA-Z_0-9]\w*'); % To check for temporary files

           if (length(startIndex2)==0) % If actual File
               
                % Incrementing FileCounter
                FileCounter=FileCounter+1               

               % Getting Real File Name
               FileName4=FileName3;
               
               % Getting Information from the File Name
               FileName4_Contents=split(FileName4,'_');
               Lat=FileName4_Contents{2};
               Long=FileName4_Contents{3};
               Year=FileName4_Contents{4};

               % Subfolder Path
               WeatherFilePath=[FolderPath2,'/',FileName4];               
               
               % Reading the Weather File
               ActualFile= csvread(WeatherFilePath,3);
               
                %% Step 2: Changing Date-Time Stamp Columns for Utility

                DateTimeStamp=zeros(1,4); % Initialization

                for i=1:length(ActualFile) % For each row in ActualFile

                   Hour=ActualFile(i,4);

                   Min=ActualFile(i,5);

                   [ TimeDeci ] = HMToDeci( Hour,Min,0 );

                   DateTimeStamp(i,1:4)=[ActualFile(i,3),ActualFile(i,2),ActualFile(i,1),TimeDeci];

                end

                ActualFile=ActualFile(:,6:end); % Removing older Date-Time Stamp

                ActualFile=horzcat(DateTimeStamp,ActualFile); % Adding New Date-Time Stamp               
               
               %% Step 6 - Preprocessing Time Resolution of the Files
               
               % Res:30-30
               [NewResFile_Res30] = NSRDB_Low2HighRes(30,30,1,ActualFile);
               
               % Res:30-15
               
               [NewResFile_Res15] = NSRDB_Low2HighRes(30,15,1,ActualFile);
               
               % Res:30-10
               [NewResFile_Res10] = NSRDB_Low2HighRes(30,10,1,ActualFile);
               
               % Res:30-5
               [NewResFile_Res5] = NSRDB_Low2HighRes(30,5,1,ActualFile);
               
               %% Step 7 - Creating New Directories and Sub-Directories
               CurrentPath=pwd;
               
               cd(DestinationFolder_Path);
               
               mkdir(CityStateName);
               
               cd([DestinationFolder_Path,'/',CityStateName]);
               
               mkdir('Res_30');
               
               mkdir('Res_15');
               
               mkdir('Res_10');
               
               mkdir('Res_5');
               
               cd(CurrentPath);
               
               %% Step 8 - Saving Processed Weather Data Files 
               % Getting the Correct File Names              
                NewResWeatherDataFileName_Res30=[CityStateName,'_Lat-',num2str(Lat),'_Long-','_',num2str(Long),num2str(Year),'_To_',num2str(Year),'_WeatherData_NSRDB_','30minTo30' ,'minRes.csv'];
                NewResWeatherDataFileName_Res15=[CityStateName,'_Lat-',num2str(Lat),'_Long-','_',num2str(Long),num2str(Year),'_To_',num2str(Year),'_WeatherData_NSRDB_','30minTo15' ,'minRes.csv'];
                NewResWeatherDataFileName_Res10=[CityStateName,'_Lat-',num2str(Lat),'_Long-','_',num2str(Long),num2str(Year),'_To_',num2str(Year),'_WeatherData_NSRDB_','30minTo10' ,'minRes.csv'];
                NewResWeatherDataFileName_Res5=[CityStateName,'_Lat-',num2str(Lat),'_Long-','_',num2str(Long),num2str(Year),'_To_',num2str(Year),'_WeatherData_NSRDB_','30minTo5' ,'minRes.csv'];

                % Getting Correct Folder Paths
                DestinationFolder_Path_Res30=['/home/ninadgaikwad/Dropbox (UFL)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/NSRDB_DATA_PROCESSING/NSRDB_USA_Database_1998_2018_ByCity_DiffRes/',CityStateName,'/Res_30'];

                DestinationFolder_Path_Res15=['/home/ninadgaikwad/Dropbox (UFL)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/NSRDB_DATA_PROCESSING/NSRDB_USA_Database_1998_2018_ByCity_DiffRes/',CityStateName,'/Res_15'];

                DestinationFolder_Path_Res10=['/home/ninadgaikwad/Dropbox (UFL)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/NSRDB_DATA_PROCESSING/NSRDB_USA_Database_1998_2018_ByCity_DiffRes/',CityStateName,'/Res_10'];

                DestinationFolder_Path_Res5=['/home/ninadgaikwad/Dropbox (UFL)/NinadGaikwad_PhD/Gaikwad_Research/Gaikwad_Research_Work/NSRDB_DATA_PROCESSING/NSRDB_USA_Database_1998_2018_ByCity_DiffRes/',CityStateName,'/Res_5'];

                
                % Writing the Processed Weather Data to appropriate Files
                csvwrite(strcat(DestinationFolder_Path_Res30,'/',NewResWeatherDataFileName_Res30),NewResFile_Res30);
                csvwrite(strcat(DestinationFolder_Path_Res15,'/',NewResWeatherDataFileName_Res15),NewResFile_Res15);
                csvwrite(strcat(DestinationFolder_Path_Res10,'/',NewResWeatherDataFileName_Res10),NewResFile_Res10);
                csvwrite(strcat(DestinationFolder_Path_Res5,'/',NewResWeatherDataFileName_Res5),NewResFile_Res5);

           end
           
        end  
        
   end
    
end

