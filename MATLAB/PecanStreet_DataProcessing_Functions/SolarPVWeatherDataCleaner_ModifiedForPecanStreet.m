function [ ProcessedData ] = SolarPVWeatherDataCleaner_ModifiedForPecanStreet( Res,DataCols,N,DataFile )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Headers = If Data File Contains Headers, then 1, otherwise 0
% Res = Time Step of the Data File in Minutes
% DataCols = Number of Columns in Data File Which Represents Data, other than date and time Columns

%% Computing Data for getting EXACT ROWS and COLUMNS of the PROCESSED DATA FILE

% Finding Start and End Days of the Data Set

% Start day Information
StartMonth=DataFile(1,2);

StartDay=DataFile(1,1);

StartYear=DataFile(1,3);

% End Day Information

[rnum,colnum]=size(DataFile);

EndMonth=DataFile(rnum,2);

EndDay=DataFile(rnum,1);

EndYear=DataFile(rnum,3);

% Computing Rows And Columns for the Processed Data File using Pre-defined Function

[ Rows,Cols,TotDays ] = RowsColsToComputeDataCleaning( StartYear,StartMonth,StartDay,EndYear,EndMonth,EndDay,Res,DataCols,4 );

% Initializing Processed Data File to zeros

ProcessedData=nan(Rows,Cols);

% Initializing Data Captur Matrix to Zeros

DataCapture=zeros(1,DataCols);

%% Putting Data into CORRECT ROWS & COLUMNS from Raw Data File to the Pre-Initialized Processed Data file

% Creating Date Time (Decimal Solar Time) Matrix for the given number of Days using Pre-Defined Function

[ DateTimeMatrix,~,TimeT ] = StartEndCalender( StartYear,StartMonth,StartDay,TotDays,Res,DataCols );

TimeT=TimeT'; % Converting Column Vector to Row Vector
len=length(TimeT); 

% Copying the DateTimeMatrix to the ProcessedData Matrix

ProcessedData=DateTimeMatrix;
[rrnum,ccnum]=size(ProcessedData);
len1=rrnum;
Up=1; % This variable will be updated whithin the following for loop to make UPDATEING FOR LOOP Dynamic to Compute FASTER

First=0; % Debugger

% Updating ProcessedData Data Columns in ProcessedData matrix for each row in Original Data File

for i=1:rnum
    
    First=First+1; % Debugger   
                 
    % Reading Date Time Signature of Current Data Row

    Month=DataFile(i,2);

    Day=DataFile(i,1);

    Year=DataFile(i,3);

    % Current Instant Time Information

    [ TimeDeci ] = DataFile(i,4);

    % Reading Data from the Current Row of DataFile into DataCapture Vector
    for k=4+1:(DataCols+4)

        DataCapture(1,k-4)=DataFile(i,k);

    end

    % Finding Corrected Time value for Time Deci as per the Time Signature of ProcessedData Matrix            
    for j=1:len

        diffrence(1,j)=abs(TimeDeci-TimeT(1,j));

    end

    [M,I]=min(diffrence);

    T=TimeT(1,I); % Corrected Time Value  

    % Computing CORRECT INDEX (Using: Day, Month, Year and T) of the ProcessesdData Matrix where the Current Data should be Stored

    for l=Up:len1

        if (Day==ProcessedData(l,1))&&(Month==ProcessedData(l,2))&&(Year==ProcessedData(l,3))&&(T==ProcessedData(l,4))

          break; 

        end

    end

    %Up=Up+1; % Updating Loop Starting Point for Faster Computation           

    % Storing Data from DataCapture Matrix to the ProcessedData Matrix at the Correct Location (Given by INDEX 'l')

    for m=1:DataCols

        ProcessedData(l,m+4)=DataCapture(1,m);

    end

end

%% N Point Average Method for Filling missing Data in TEMPERATURE, WIND, RELATIVE HUMIDITY

RaN=zeros(N,DataCols); % Initializing the Running Value Storage used in NaN Value Filling 

NPointAverageN=sum(RaN)/N; % Calculating the Running Average

Second=0; % Debugger

% FOR LOOP for Point-Wise Filling of NaN and Zero Values (Top to Bottom)

for i=1:Rows
    
    Second=Second+1; % Debugger
    
    RaCounter=rem(i,N); % For using Running Value Storage Vector to cyclically update its N Values with the next value
    
    if RaCounter==0
        
        RaCounter=N;
        
    end
    
    % FOR LOOP for Each Data Column
    
    for k=1:DataCols
        
        if (isnan(ProcessedData(i,k+4))) % ||(ProcessedData(i,k+4)==0
            
            NPointAverageN=sum(RaN(:,k))/N; % Calculating the Running Average
            
            ProcessedData(i,k+4)=NPointAverageN; % Updating NaN Value with Running Average Value
        
            RaN(RaCounter,k)=ProcessedData(i,k+4);
            
        else
            
            RaN(RaCounter,k)=ProcessedData(i,k+4);
            
        end
        
           

        
        
    end
    
    
    
end

% FOR LOOP for Point-Wise Filling of NaN and Zero Values (Bottom to Top)

Third=0; % Debugger

for i=Rows:-1:1
    
    Third=Third+1; % Debugger
    
    RaCounter=rem(i,N); % For using Running Value Storage Vector to cyclically update its N Values with the next value
    
    if RaCounter==0
        
        RaCounter=N;
        
    end
    
    % FOR LOOP for Each Data Column
    
    for k=1:DataCols
        
        if (isnan(ProcessedData(i,k+4))) % ||(ProcessedData(i,k+4)==0
            
            NPointAverageN=sum(RaN(:,k))/N; % Calculating the Running Average
            
            ProcessedData(i,k+4)=NPointAverageN; % Updating NaN Value with Running Average Value
        
            RaN(RaCounter,k)=ProcessedData(i,k+4);
            
        else
            
            RaN(RaCounter,k)=ProcessedData(i,k+4);
            
        end
        
           

        
        
    end
    
    
    
end



end



