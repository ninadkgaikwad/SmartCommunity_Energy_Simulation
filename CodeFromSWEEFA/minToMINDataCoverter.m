function [ ProcessedData1 ] = minToMINDataCoverter( DataCols,ResOriginal,ResNew,AvgOrAdd,Headers )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Note: ResOriginal Should Always Be Smaller Than ResNew; And 60mins Should be integral multiples of both

% Note: Same Function can be used to convert files to Half-Hourly or Hourly Files

%% Getting the Raw Data File in the MATLAB WorkSpace

% File Selection
[Filename,Pathname]=uigetfile({'*.*'},'Raw Data File Selector');

Fullpathname=strcat(Pathname,Filename);

[ProcessedData]=xlsread(Fullpathname,1);

%% Computing Size of ProcessedData Matrix

% [Row1,Col1]=size(ProcessedData);
% 
% ProcessedData(Row1+1,:)=ProcessedData(Row1,:); % For avoiding the problem with starting with 0 time field
% 
% ProcessedData(Row1+1,4)=0; % For avoiding the problem with starting with 0 time field

[Row,Col]=size(ProcessedData);

RowNew=Row*(ResOriginal/ResNew); % [Embedded Formula] Computing Number of Rows for the PRocessedData1 Matrix

RowNew=ceil(RowNew);

% Initializing The New ProcessedData Matrix

ProcessedData1=zeros(RowNew,(4+DataCols));

% Computing Number off ROWS to be AVERAGED or ADDED

NumRows=ResNew/ResOriginal;

% [Modification] Recoding for Fractional NumRows Problem

NumRows_Frac = rem(ResNew, ResOriginal); % If NumRows_Frac==0 , There is no Fractional NumRows Problem (Use Earlier Algorithm) ; NumRows_Frac~=0 , There is Fractional NumRows Problem (Use New Algorithm)

if (Headers==1)
    
    % Getting Headers
    
    [~ ,~,DataFile]=xlsread(Fullpathname,1);
    
    % Getting the Header Text
    
    Header1 = DataFile(1,5:(DataCols+4));    
    
    % Clearing the not needed DataFile Variable
    
    clearvars DataFile
    
    Header = {'Day', 'Month', 'Year', 'Time'};
    
    % Concatenating Headers derived frome the Original File
    
    Header = [Header, Header1];    
    
elseif (Headers==0)
    
    Header = {'Day', 'Month', 'Year', 'Time'};
    
end

if (NumRows_Frac == 0) % Use Earlier Algorithm

    % Initializing Index for ProcessedData1 Matrix

    Index1=1;

    %% FOR LOOP for Averging and Adding to get Desired ResNew according to AvgOrAdd

    % FOR LOOP for each ROW of ProcessedData Skipping by NumRows

    for i=2:NumRows:Row

        % Correcting for the 0 Hour field on the First Row of the ProcessedData Matrix

        ProcessedData1(1,:)=ProcessedData(1,:); % Copying Correct DateTime Stamp

        % Incrementing Index1 for placing Data in Correct Rows of ProcessedData1 Matrix

        Index1 =Index1+1;

        % FOR LOOP for each DataCol

        for k=1:DataCols

            Indicator=AvgOrAdd(1,k); % For indication Values Should be Averaged or Added

            Add=0; % Initializing Add Variable to Zero

            Avg=zeros(1,NumRows); % Initializing Avg Vector to Zeros

            % FOR LOOP for Averaging & Adding as Per RES Values

            for j=1:NumRows

                RowIndex=i+(j-1); % Computing RowIndex

                if Indicator==1 % ADDITION

                    Add=Add+ProcessedData(RowIndex,(k+4));

                elseif Indicator==0 % AVERAGE

                    Avg(1,j)=ProcessedData(RowIndex,(k+4));

                end

            end        

            ProcessedData1(Index1,1:4)=ProcessedData(RowIndex,1:4); % Copying Correct DateTime Stamp

            if Indicator==1 % ADDITION

                 ProcessedData1(Index1,(k+4))=Add; % Copying Correct Data Value

            elseif Indicator==0 % AVERAGE

                 ProcessedData1(Index1,(k+4))=sum(Avg)/NumRows; % Copying Correct Data Value

            end

        end


    end
    
elseif (NumRows_Frac ~= 0) % [Modification] Use New Algorithm
    
    % Getting StartYear,StartMonth,StartDay,EndYear,EndMonth,EndDay from the Input File
    
    StartYear = ProcessedData(1,3);   
    StartMonth = ProcessedData(1,2);
    StartDay = ProcessedData(1,1);
    
    EndYear = ProcessedData(Row,3);
    EndMonth = ProcessedData(Row,2);
    EndDay = ProcessedData(Row,1);
    
    % Computing Rows And Columns for the Processed Data File using Pre-defined Function

    [ Rows1,Cols1,TotDays ] = RowsColsToComputeDataCleaning( StartYear,StartMonth,StartDay,EndYear,EndMonth,EndDay,ResNew,DataCols,4 );
    
    % Initializing Processed Data File to zeros

    ProcessedData1=zeros(Rows1,Cols1);
    
    % Creating Date Time (Decimal Time) Matrix for the given number of Days using Pre-Defined Function

    [ DateTimeMatrix,~,TimeT] = StartEndCalender( StartYear,StartMonth,StartDay,TotDays,ResNew,DataCols );
    
    % Filling the ProcessedData1 Matrix with Correct Time-Stamps
    
    ProcessedData1(:,1:4)=DateTimeMatrix(:,1:4);
    
    %% FOR LOOP for Averging and Adding to get Desired ResNew according to AvgOrAdd
    
    Index1=0; % Initializing The Row Counter for the New ProcessedData1 Matrix

    % FOR LOOP for each ROW of ProcessedData Skipping by NumRows      
    
    for i=1:NumRows:Row
        
        % Incrementing Index1
        
        Index1 = Index1+1;
        
        % Storing Previous Values of i
        
        RowNumVector(1,Index1)=i;    
        
        % Getting Current NumRow
        
        NumRow = i;
        
        
        % Correcting for the 0 Hour field on the First Row of the ProcessedData Matrix
        
        if (i==1) % The First Row
            
            ProcessedData1(1,:)=ProcessedData(1,:); % Copying the First Row as it is 
            
            continue; % Advance to Next Iteration
            
        end
        
        % Getting the Previous Row Num Value
        
        Previous_RowNum = RowNumVector(1,(Index1-1));
        
        % Finding whether the Previous_RowNum is Integer or Fractional

        IntFrac_Indicator2 = rem(Previous_RowNum,1); % Can be used for Additio Weight
        
        % Computing Num Row Previous based on the Previous_RowNum Value
        
        if (IntFrac_Indicator2==0) % Previous_NumRow is Integer
            
            NumRow_Previous = Previous_RowNum+1;            
            
            StartVal_Multiplier = 1;
            
        elseif (IntFrac_Indicator2~=0) % Previou_NumRow is Fractional
            
            NumRow_Previous = ceil(Previous_RowNum);
            
            StartVal_Multiplier = 1-IntFrac_Indicator2;
            
        end
        
        % Finding whether the current Num Row is Integer or Fractional
        
        IntFrac_Indicator1 = rem(i,1); % Can be used for Addition Weight
        
        % Computing Next Num Row based on the IntFrac_Indicator1
        
        if (IntFrac_Indicator1==0) % Current Num Row is Integer
            
            NumRow_Next = NumRow;
            
            EndVal_Multiplier = 1;
            
        elseif (IntFrac_Indicator1~=0) % Current Num Row is Fractional
            
            NumRow_Next = ceil(NumRow);
            
            EndVal_Multiplier = IntFrac_Indicator1;
            
        end
        
        % Getting the Actual Indices of the Original Data Set
        
        Actual_Indices = [NumRow_Previous:NumRow_Next];
        
        % FOR LOOP for each DataCol
        
        for k=1:DataCols

            Indicator=AvgOrAdd(1,k); % For indication Values Should be Averaged or Added
            
            % Collecting the Desired Values from the Current Columns from the Original Data Set
            
            Actual_Values=ProcessedData(Actual_Indices,(k+4));
            
            % Using the Indicator for Either Averaging or Addition of the Values
            
            if (Indicator==1)  % Perform Addition on the Desired Values
                
                % Computing Weighted Sum
                
                Actual_Values(1,1) = StartVal_Multiplier*Actual_Values(1,1);  % Correcting Addition Biases due to Fractional NumRows
                
                Actual_Values(end,1) = EndVal_Multiplier*Actual_Values(end,1); % Correcting Addition Biases due to Fractional NumRows
                
                Addition_Value = sum(Actual_Values);                
                
                % Putting the Addition_Value in the Correect Cell of the ProcessedData1 matrix
                
                ProcessedData1(Index1,(k+4)) = Addition_Value;
                
            elseif (Indicator==0) % Perform Averaging on the Desired Values
                
                % Computing Average
                
                Averaged_Value = sum(Actual_Values)/length(Actual_Values);
                
                % Putting the Averaged_Value in the Correect Cell of the ProcessedData1 matrix
                
                ProcessedData1(Index1,(k+4)) = Averaged_Value;               
                
            end
            
            
        end    
       
               
     end
    
end

%% Writing the ProcessedData Matrix to an Excel File

% Creating a Proper File Name

Filename1 = strread(Filename,'%s','delimiter','_'); % Using '_' as Delimiter for removing File Extensions

Filename = Filename1{1,1};

filename = [Filename,'_Converted_File_MinutesResolution_',num2str(ResOriginal),'-Mins_To_',num2str(ResNew),'-Mins','.xlsx'];

% Wrting the Header

if (Headers==1)
    
    % Getting Headers
    
    [~ ,~,DataFile]=xlsread(Fullpathname,1);
    
    % Getting the Header Text
    
    Header1 = DataFile(1,5:(DataCols+4));    
    
    % Clearing the not needed DataFile Variable
    
    clearvars DataFile
    
    Header = {'Day', 'Month', 'Year', 'Time'};
    
    % Concatenating Headers derived frome the Original File
    
    Header = [Header, Header1];    
    
elseif (Headers==0)
    
    Header = {'Day', 'Month', 'Year', 'Time'};
    
end

sheet = 1;

xlRange = 'A1';

xlswrite(filename,Header,sheet,xlRange);

% Wrting the Cleaned Data

sheet = 1;

xlRange = 'A2';

xlswrite(filename,ProcessedData1,sheet,xlRange);


end

