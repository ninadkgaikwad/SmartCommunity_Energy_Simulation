function [PecanStreet_Data_Output] = PecanStreet_Data_Extractor(HEMSWeatherData_Input,Simulation_Params,PecanStreet_Data_FolderPath,N_House_Vector,Type,Data_MatFile_Name)

%% Pecan Street Data Extraction : Converting Many Files to 1

% Author: Ninad Kiran Gaikwad
% Date: Mar/15/2021
% Description: Pecan Street Data Extractor

%% Getting Required data from HEMSWeatherData_Input - Struct

%---------------------- HEMSWeatherData_Input ----------------------------%
StartMonth=HEMSWeatherData_Input.StartMonth;
StartDay=HEMSWeatherData_Input.StartDay;
StartTime=HEMSWeatherData_Input.StartTime;
EndMonth=HEMSWeatherData_Input.EndMonth;
EndDay=HEMSWeatherData_Input.EndDay;
EndTime=HEMSWeatherData_Input.EndTime;

%---------------------- HEMSWeatherData_Output ---------------------------%
%DateTime_Matrix=HEMSWeatherData_MPC_Output.DateTime_Matrix;

%------------------------- Simulation_Params------------------------------%
FileRes=Simulation_Params.FileRes;

%% Getting All Files Data

DataFolderName=PecanStreet_Data_FolderPath;

DataFolderStructure=dir(DataFolderName);

FileNames_Cell=extractfield(DataFolderStructure,'name');

%FolderPath=DataFolderStructure.folder;

File_Counter=0; % Initialization

for ii=1:length(FileNames_Cell) % For each file in the Folder
    
    FileName=FileNames_Cell{ii};
    
    startIndex = regexp(FileName,'\<[^a-zA-Z_0-9]\w*'); % To check for temporary files
    
   if (isempty(startIndex)) % If actual File
             
       File_Counter=File_Counter+1; % Incrementing File_Counter
       
       FullPathName=strcat(DataFolderName,FileName); % Creating Full Path for the Actual File
           
       AllFiles_Data_Cell(1,File_Counter)= {csvread(FullPathName)};            
       
   end
    
end

%% Getting All files which have the required DateTime Range

DateTime_NoError_Counter=0; % Initialization

%AllFiles_CorrectDates_Data_Matrix=NaN; % Initialization

for ii=1:length(AllFiles_Data_Cell) % For each file AllFiles_Data_Cell
    
    Single_DataFile_Matrix=AllFiles_Data_Cell{1,ii};
    
    % Getting StartYear and EndYear
    StartYear=Single_DataFile_Matrix(1,3);
    EndYear=StartYear;
    
    % Creating a DateTimeMatrix_ForSlicer
    DateTimeMatrixAggregate_ForSlicer=horzcat(Single_DataFile_Matrix(:,1:4),zeros(length(Single_DataFile_Matrix(:,1)),1));

    % Computing required Start and End Indices
    [ OriginalSeries,StartIndex_Aggregate,EndIndex_Aggregate ] = DateTimeSeriesSlicer_PecanStreetData(DateTimeMatrixAggregate_ForSlicer,1,FileRes,StartYear,EndYear,StartMonth,EndMonth,StartDay,EndDay,StartTime,EndTime);

    % Checking if desired Dates were presentin the current file
    if (isstring(OriginalSeries)) % Desired Dates not present in the current file
        
        % Nothing Happens - We skip over the current file
        
    else % Desired Dates present in the current file
        
        % Getting the Data for the desired dates
        Data=Single_DataFile_Matrix(StartIndex_Aggregate:EndIndex_Aggregate,:);
        
        % Incrementing DateTime_NoError_Counter
        DateTime_NoError_Counter=DateTime_NoError_Counter+1;
        
        % Storing the Files with correct dates along with DateTime
        AllFiles_CorrectDates_Data_Matrix(:,:,DateTime_NoError_Counter)=Data;
        
    end     
    
end

%% Arranging Data Columns in Descending order of SolarPV, Battery, EV, Priority for Loads

if (isnan(AllFiles_CorrectDates_Data_Matrix))
    
    warning('Desired Dates not found in files')
    
    PecanStreet_Data_Output='None';
    
else
    
    AllFiles_CorrectDates_Priority_Data_Matrix=[AllFiles_CorrectDates_Data_Matrix(:,1:4,:)... % Date-Time
                                                AllFiles_CorrectDates_Data_Matrix(:,66+4:67+4,:) AllFiles_CorrectDates_Data_Matrix(:,13+4,:) AllFiles_CorrectDates_Data_Matrix(:,14+4:15+4,:)... % Solar PV, Battery, EV
                                                AllFiles_CorrectDates_Data_Matrix(:,61+4:62+4,:) AllFiles_CorrectDates_Data_Matrix(:,25+4,:) AllFiles_CorrectDates_Data_Matrix(:,39+4,:) AllFiles_CorrectDates_Data_Matrix(:,37+4,:) AllFiles_CorrectDates_Data_Matrix(:,40+4,:) AllFiles_CorrectDates_Data_Matrix(:,38+4,:) AllFiles_CorrectDates_Data_Matrix(:,60+4,:) AllFiles_CorrectDates_Data_Matrix(:,71+4,:) AllFiles_CorrectDates_Data_Matrix(:,49+4,:) AllFiles_CorrectDates_Data_Matrix(:,53+4:54+4,:) ... Level 1 Priority - Fridge, Freezer, Kitchen
                                                AllFiles_CorrectDates_Data_Matrix(:,8+4:12+4,:)... % Level 2 Priority - Bedrooms
                                                AllFiles_CorrectDates_Data_Matrix(:,47+4:48+4,:) AllFiles_CorrectDates_Data_Matrix(:,50+4,:) ... % Level 3 Priority - Living Rooms, Office Room
                                                AllFiles_CorrectDates_Data_Matrix(:,17+4,:) AllFiles_CorrectDates_Data_Matrix(:,23+4,:) AllFiles_CorrectDates_Data_Matrix(:,64+4,:) AllFiles_CorrectDates_Data_Matrix(:,22+4,:) AllFiles_CorrectDates_Data_Matrix(:,59+4,:) AllFiles_CorrectDates_Data_Matrix(:,69+4,:) AllFiles_CorrectDates_Data_Matrix(:,74+4,:) ... % Level 4 Priority Clothes, Garbage Disposal, P{umps
                                                AllFiles_CorrectDates_Data_Matrix(:,63+4,:) AllFiles_CorrectDates_Data_Matrix(:,6+4:7+4,:) AllFiles_CorrectDates_Data_Matrix(:,19+4:21+4,:) ... % Level 5 Priority - Security, Bathrooms, Dinning Room, Dishwasher
                                                AllFiles_CorrectDates_Data_Matrix(:,28+4:29+4,:) AllFiles_CorrectDates_Data_Matrix(:,70+4,:) AllFiles_CorrectDates_Data_Matrix(:,65+4,:) AllFiles_CorrectDates_Data_Matrix(:,51+4:52+4,:) ... % Level 6 Priority - Remaining Rooms, Outside Lights
                                                AllFiles_CorrectDates_Data_Matrix(:,5+4,:) AllFiles_CorrectDates_Data_Matrix(:,68+4,:) AllFiles_CorrectDates_Data_Matrix(:,75+4,:) ... % Level 7 Priority - Aquarium, Lawn Sprinklers, Wine Cooler
                                                AllFiles_CorrectDates_Data_Matrix(:,35+4,:) AllFiles_CorrectDates_Data_Matrix(:,58+4,:) AllFiles_CorrectDates_Data_Matrix(:,57+4,:) AllFiles_CorrectDates_Data_Matrix(:,36+4,:) AllFiles_CorrectDates_Data_Matrix(:,72+4:73+4,:)]; % Level 8 Priority - Pool, Jacuzzi, Water Heater

    %% Arranging AllFiles_CorrectDates_Priority_Data_Matrix as per [N_PV_Bat_EV, N_PV_Bat, N_PV_EV, N_Bat_EV, N_PV, N_Bat, N_EV, N_None]

    % Getting number of Files in AllFiles_CorrectDates_Data_Matrix
    [Rows,Columns,File_Num]=size(AllFiles_CorrectDates_Data_Matrix);

    % Computing Info for PV File Indices
    House_PVSum_Array=sum(AllFiles_CorrectDates_Priority_Data_Matrix(:,5:6,:));
    House_PV_File_Indices=find(House_PVSum_Array(1,1,:));

    % Computing Info for Bat File Indices
    House_BatSum_Array=sum(AllFiles_CorrectDates_Priority_Data_Matrix(:,13,:));
    House_Bat_File_Indices=find(House_BatSum_Array(1,1,:));

    % Computing Info for EV File Indices
    House_EVSum_Array=sum(AllFiles_CorrectDates_Priority_Data_Matrix(:,14:15,:));
    House_EV_File_Indices=find(House_EVSum_Array(1,1,:));

    % PV_Bat_EV File Indices Intersection
    House_PV_Bat_EV_File_Indices_Intersection=intersect(intersect(House_PV_File_Indices,House_Bat_File_Indices),House_EV_File_Indices);

    % PV_Bat_EV File Indices Union
    House_PV_Bat_EV_File_Indices_Union=union(union(House_PV_File_Indices,House_Bat_File_Indices),House_EV_File_Indices);

    % Only None Indices
    House_OnlyNone_File_Indices=setdiff(1:File_Num,House_PV_Bat_EV_File_Indices_Union);

    % Only PV Indices
    House_OnlyPV_File_Indices=setdiff(setdiff(setdiff(1:File_Num,House_OnlyNone_File_Indices),House_Bat_File_Indices),House_EV_File_Indices);

    % Only Bat Indices
    House_OnlyBat_File_Indices=setdiff(setdiff(setdiff(1:File_Num,House_OnlyNone_File_Indices),House_PV_File_Indices),House_EV_File_Indices);

    % Only EV Indices
    House_OnlyEV_File_Indices=setdiff(setdiff(setdiff(1:File_Num,House_OnlyNone_File_Indices),House_PV_File_Indices),House_Bat_File_Indices);

    % Only PV_Bat Indices
    House_Only_PV_Bat_File_Indices=intersect(House_OnlyPV_File_Indices,House_OnlyBat_File_Indices);

    % Only PV_EV Indices
    House_Only_PV_EV_File_Indices=intersect(House_OnlyPV_File_Indices,House_OnlyEV_File_Indices);

    % Only Bat_EV Indices
    House_Only_Bat_EV_File_Indices=intersect(House_OnlyBat_File_Indices,House_OnlyEV_File_Indices);

    %% Getting required House Data - Based on Type of Community [PV,Bat,EV,None]

    % Getting Length of N_House_Vector - Gives information about what Major devices are present
    Length_HouseVector=length(N_House_Vector);

    % Creating File Data Matrix based on Single Large House (Type=[1,2,3,4,5,6,7,8]) and Community of Houses (Type=[0])
    if (Length_HouseVector==1) % Single Large House [N_PV_Bat_EV, N_PV_Bat, N_PV_EV, N_Bat_EV, N_PV, N_Bat, N_EV, N_None]

        switch Type
            case 1 % N_PV_Bat_EV
                if (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_PV_Bat_EV_File_Indices_Intersection(1));
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    warning('Bat_EV File instead of PV_Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_Bat_EV_File_Indices(1));
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    warning('PV_EV File instead of PV_Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_EV_File_Indices(1));                
                elseif (~isempty(House_OnlyEV_File_Indices))
                    warning('EV File instead of PV_Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyEV_File_Indices(1));                    
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    warning('PV_Bat File instead of PV_Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_Bat_File_Indices(1));                        
                elseif (~isempty(House_OnlyBat_File_Indices))
                    warning('Bat File instead of PV_Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyBat_File_Indices(1));                            
                elseif (~isempty(House_OnlyPV_File_Indices))
                    warning('PV File instead of PV_Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyPV_File_Indices(1));                                
                elseif (~isempty(House_OnlyNone_File_Indices))
                    warning('None File instead of PV_Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyNone_File_Indices(1));                                    
                else 
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output=zeros(Rows,Columns);                
                end
            case 2 % N_PV_Bat
                if (~isempty(House_Only_PV_Bat_File_Indices))
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_Bat_File_Indices(1));
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    warning('PV_Bat_EV File instead of PV_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_PV_Bat_EV_File_Indices_Intersection(1));
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    warning('Bat_EV File instead of PV_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_Bat_EV_File_Indices(1));                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    warning('PV_EV File instead of PV_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_EV_File_Indices(1));                        
                elseif (~isempty(House_OnlyBat_File_Indices))
                    warning('Bat File instead of PV_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyBat_File_Indices(1));                            
                elseif (~isempty(House_OnlyPV_File_Indices))
                    warning('PV File instead of PV_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyPV_File_Indices(1)); 
                elseif (~isempty(House_OnlyEV_File_Indices))
                    warning('EV File instead of PV_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyEV_File_Indices(1));                     
                elseif (~isempty(House_OnlyNone_File_Indices))
                    warning('None File instead of PV_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyNone_File_Indices(1));                                    
                else 
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output=zeros(Rows,Columns);                
                end
            case 3 % N_PV_EV
                if (~isempty(House_Only_PV_EV_File_Indices))
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_EV_File_Indices(1));
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    warning('PV_Bat_EV File instead of PV_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_PV_Bat_EV_File_Indices_Intersection(1));
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    warning('Bat_EV File instead of PV_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_Bat_EV_File_Indices(1));                    
                elseif (~isempty(House_OnlyEV_File_Indices))
                    warning('EV File instead of PV_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyEV_File_Indices(1));                            
                elseif (~isempty(House_OnlyPV_File_Indices))
                    warning('PV File instead of PV_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyPV_File_Indices(1));  
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    warning('PV_Bat File instead of PV_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_Bat_File_Indices(1));                       
                elseif (~isempty(House_OnlyBat_File_Indices))
                    warning('Bat File instead of PV_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyBat_File_Indices(1));                     
                elseif (~isempty(House_OnlyNone_File_Indices))
                    warning('None File instead of PV_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyNone_File_Indices(1));                                    
                else 
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output=zeros(Rows,Columns);                
                end
            case 4 % N_Bat_EV
                if (~isempty(House_Only_Bat_EV_File_Indices))
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_Bat_EV_File_Indices(1));
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    warning('PV_Bat_EV File instead of Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_PV_Bat_EV_File_Indices_Intersection(1));
                elseif (~isempty(House_OnlyEV_File_Indices))
                    warning('EV File instead of Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyEV_File_Indices(1));  
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    warning('PV_EV File instead of Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_EV_File_Indices(1));                  
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    warning('PV_Bat File instead of Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_Bat_File_Indices(1));                       
                elseif (~isempty(House_OnlyBat_File_Indices))
                    warning('Bat File instead of Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyBat_File_Indices(1));  
                elseif (~isempty(House_OnlyPV_File_Indices))
                    warning('PV File instead of Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyPV_File_Indices(1));                  
                elseif (~isempty(House_OnlyNone_File_Indices))
                    warning('None File instead of Bat_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyNone_File_Indices(1));                                    
                else 
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output=zeros(Rows,Columns);                
                end                
            case 5 % N_PV
                if (~isempty(House_OnlyPV_File_Indices))
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyPV_File_Indices(1));
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    warning('PV_Bat_EV File instead of N_PV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_PV_Bat_EV_File_Indices_Intersection(1));                  
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    warning('PV_Bat File instead of N_PV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_Bat_File_Indices(1)); 
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    warning('PV_EV File instead of N_PV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_EV_File_Indices(1));  
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    warning('Bat_EV File instead of N_PV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_Bat_EV_File_Indices(1));                      
                elseif (~isempty(House_OnlyBat_File_Indices))
                    warning('Bat File instead of N_PV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyBat_File_Indices(1));  
                elseif (~isempty(House_OnlyEV_File_Indices))
                    warning('EV File instead of N_PV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyEV_File_Indices(1));                  
                elseif (~isempty(House_OnlyNone_File_Indices))
                    warning('None File instead of N_PV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyNone_File_Indices(1));                                    
                else 
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output=zeros(Rows,Columns);               
                end 
            case 6 % N_Bat                      
                if (~isempty(House_OnlyBat_File_Indices))
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyBat_File_Indices(1)); 
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    warning('PV_Bat_EV File instead of N_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_PV_Bat_EV_File_Indices_Intersection(1));                  
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    warning('PV_Bat File instead of N_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_Bat_File_Indices(1));   
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    warning('Bat_EV File instead of N_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_Bat_EV_File_Indices(1)); 
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    warning('PV_EV File instead of N_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_EV_File_Indices(1));
                elseif (~isempty(House_OnlyPV_File_Indices))
                    warning('PV File instead of N_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyPV_File_Indices(1));
                elseif (~isempty(House_OnlyEV_File_Indices))
                    warning('EV File instead of N_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyEV_File_Indices(1));                  
                elseif (~isempty(House_OnlyNone_File_Indices))
                    warning('None File instead of N_Bat');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyNone_File_Indices(1));                                    
                else 
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output=zeros(Rows,Columns);                
                end 
            case 7 % N_EV
                if (~isempty(House_OnlyEV_File_Indices))
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyEV_File_Indices(1)); 
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    warning('PV_Bat_EV File instead of N_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_PV_Bat_EV_File_Indices_Intersection(1));  
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    warning('PV_EV File instead of N_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_EV_File_Indices(1));   
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    warning('Bat_EV File instead of N_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_Bat_EV_File_Indices(1));                 
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    warning('PV_Bat File instead of N_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_Bat_File_Indices(1));
                elseif (~isempty(House_OnlyPV_File_Indices))
                    warning('PV File instead of N_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyPV_File_Indices(1));
                elseif (~isempty(House_OnlyBat_File_Indices))
                    warning('Bat File instead of N_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyBat_File_Indices(1));                  
                elseif (~isempty(House_OnlyNone_File_Indices))
                    warning('None File instead of N_EV');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyNone_File_Indices(1));                                    
                else 
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output=zeros(Rows,Columns);                
                end 
            case 8 % N_None
                if (~isempty(House_OnlyNone_File_Indices))
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyNone_File_Indices(1)); 
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    warning('PV_Bat_EV File instead of None');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_PV_Bat_EV_File_Indices_Intersection(1));                   
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    warning('PV_Bat File instead of None');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_Bat_File_Indices(1));
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    warning('PV_EV File instead of None');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_PV_EV_File_Indices(1));   
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    warning('Bat_EV File instead of None');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_Only_Bat_EV_File_Indices(1));
                elseif (~isempty(House_OnlyPV_File_Indices))
                    warning('PV File instead of None');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyPV_File_Indices(1));
                elseif (~isempty(House_OnlyBat_File_Indices))
                    warning('Bat File instead of None');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyBat_File_Indices(1));                  
                elseif (~isempty(House_OnlyEV_File_Indices))
                    warning('EV File instead of None');
                    PecanStreet_Data_Output=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,House_OnlyEV_File_Indices(1));                                    
                else 
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output=zeros(Rows,Columns);                
                end 
        end    

    elseif (Length_HouseVector==4) % Smart Community - [N_PV_Bat, N_Bat, N_PV, N_None]

        % Getting Total number of Houses
        N_House=sum(N_House_Vector);
        
        % Initializing Constants
        N1=N_House_Vector(1);
        N2=N1+N_House_Vector(2);
        N3=N2+N_House_Vector(3);
        N4=N_House;
        
        Counter_SmartCommunity_1=0; % Initializing Constants
        
        % Creating PecanStreet_Data_Output
        for ii=1:N_House % For each house in N_House
            
            if ((ii>=1) && (ii<=N1)) % N_PV_Bat
                    
                if (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);                        
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1);                            
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1); 
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                     
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end  
                
            elseif ((ii>=N1+1) && (ii<=N2)) % N_Bat
                
                if (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);                       
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);                             
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1); 
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                     
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end                  
                
            elseif ((ii>=N2+1) && (ii<=N3)) % N_PV
                    
                if (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);                       
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                             
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1); 
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                     
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end                     
                    
            elseif ((ii>=N3+1) && (ii<=N4)) % N_None
                
                if (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);                       
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                     
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1);                              
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1); 
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                                   
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end                 
                
            end
            
        end

    elseif (Length_HouseVector==8) % Smart Community - [N_PV_Bat_EV, N_PV_Bat, N_PV_EV, N_Bat_EV, N_PV, N_Bat, N_EV, N_None]

        % Getting Total number of Houses
        N_House=sum(N_House_Vector);
        
        % Initializing Constants
        N1=N_House_Vector(1);
        N2=N1+N_House_Vector(2);
        N3=N2+N_House_Vector(3);
        N4=N3+N_House_Vector(4);
        N5=N4+N_House_Vector(5);
        N6=N5+N_House_Vector(6);
        N7=N6+N_House_Vector(7);
        N8=N_House;
        
        Counter_SmartCommunity_1=0; % Initializing Constants
        
        % Creating PecanStreet_Data_Output
        for ii=1:N_House % For each house in N_House
            
            if ((ii>=1) && (ii<=N1)) % N_PV_Bat_EV                   
    
                if (~isempty(House_Only_PV_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_EV_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of N_PV_Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_File_Indices_Intersection))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of N_PV_Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of N_PV_Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);  
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of N_PV_Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                            
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of N_PV_Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1);                        
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of N_PV_Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1);                    
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of N_PV_Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end              
            
            elseif ((ii>=N1+1) && (ii<=N2)) % N_PV_Bat
                    
                if (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);                        
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1);                            
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1); 
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                     
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of PV_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end  
                
            elseif ((ii>=N2+1) && (ii<=N3)) % N_PV_EV
                    
                if (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of PV_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of PV_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);  
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of PV_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of PV_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);                            
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of PV_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1);                         
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of PV_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1);                   
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of PV_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end 
                
            elseif ((ii>=N3+1) && (ii<=N4)) % N_Bat_EV
                    
                if (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);  
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);                           
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1);                           
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1);                  
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of Bat_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end   
                
            elseif ((ii>=N4+1) && (ii<=N5)) % N_PV
                    
                if (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);                       
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                             
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1); 
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                     
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of N_PV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end                 
                
            elseif ((ii>=N5+1) && (ii<=N6)) % N_Bat
                
                if (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);                       
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);                             
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1); 
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                     
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of N_Bat');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end 
                
            elseif ((ii>=N5+1) && (ii<=N6)) % N_EV
                
                if (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of N_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);                     
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of N_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);     
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of N_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                 
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of N_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);                             
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of N_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1); 
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of N_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1);                     
                elseif (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('None File instead of N_EV');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);                                    
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end                 
                
            elseif ((ii>=N7+1) && (ii<=N8)) % N_None
                
                if (~isempty(House_OnlyNone_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyNone_File_Indices))+1);
                elseif (~isempty(House_PV_Bat_EV_File_Indices_Intersection))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat_EV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_PV_Bat_EV_File_Indices_Intersection))+1);                       
                elseif (~isempty(House_Only_PV_Bat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_Bat File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_Bat_File_Indices))+1);                    
                elseif (~isempty(House_Only_PV_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV_EV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_PV_EV_File_Indices))+1);
                elseif (~isempty(House_Only_Bat_EV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat_EV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_Only_Bat_EV_File_Indices))+1);                     
                elseif (~isempty(House_OnlyPV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('PV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyPV_File_Indices))+1);                              
                elseif (~isempty(House_OnlyBat_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('Bat File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyBat_File_Indices))+1); 
                elseif (~isempty(House_OnlyEV_File_Indices))
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('EV File instead of N_None');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=AllFiles_CorrectDates_Priority_Data_Matrix(:,:,rem(ii,length(House_OnlyEV_File_Indices))+1);                                   
                else 
                    Counter_SmartCommunity_1=Counter_SmartCommunity_1+1; % Incrementing Counter_SmartCommunity_1
                    warning('No File found for the desired dates and required type');
                    PecanStreet_Data_Output(:,:,Counter_SmartCommunity_1)=zeros(Rows,Columns);                
                end                 
                
            end
            
        end

    end                                            
                                            
end

%% Saving PecanStreet_Data_Output as .mat File

save(Data_MatFile_Name,'PecanStreet_Data_Output');

end

