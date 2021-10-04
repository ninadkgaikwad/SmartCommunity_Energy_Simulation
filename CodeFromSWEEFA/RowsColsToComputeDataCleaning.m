function [ Rows,Cols,TotDays ] = RowsColsToComputeDataCleaning( StartYear,StartMonth,StartDay,EndYear,EndMonth,EndDay,Res,DataCols,DateTimeCols )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Computing the Total number of Different Years in the Data
NumOfYears=EndYear-StartYear+1;

% Computing the Different Year Signature Values
for i=1:NumOfYears
    
    Year(1,i)=StartYear+i-1;
    
end

% Finding the Leap and Non-Leap Year
LeapYear= LeapYearFinder( Year ); 

% Initializing Day Counters
a=0;
b=0;
c=zeros(1,NumOfYears);

% Computing Number of days in the given Data Set
for j=1:NumOfYears
    
    if j==1 % Days for Start Year
        
        if NumOfYears==1
            
            [ StartDay1, EndDay1 ] = DaysToCompute( LeapYear(1,j),StartDay,StartMonth,EndDay,EndMonth );
            
            a=EndDay1-StartDay1+1; %Total Numbe of Days
            
        else
            
            [ StartDay1, EndDay1 ] = DaysToCompute( LeapYear(1,j),StartDay,StartMonth,31,12 );
            
            a=EndDay1-StartDay1+1; %Total Numbe of Days
        
        end
        
    elseif j==NumOfYears % Days for End Year
        
       [ StartDay1, EndDay1 ] = DaysToCompute( LeapYear(1,j),1,1,EndDay,EndMonth ) ;
       
       b=EndDay1-StartDay1+1; %Total Numbe of Days
        
    else % Days for all other Years  
        
        [ StartDay1, EndDay1 ] = DaysToCompute( LeapYear(1,j),1,1,31,12 );
        
        c(1,j)=EndDay1-StartDay1+1; %Total Numbe of Days
        
    end
    
end

% Total Number of Days in the Data Set
TotDays=a+b;

for k=1:NumOfYears
    
TotDays=TotDays+c(1,k);

end

% Data Points in ONE DAY (Resolution has to be in Minutes)
DataPoints=24*(60/Res);

% Total Data Points in the Given Data Set i.e the Total Number of Rows
Rows=TotDays*DataPoints;

% Total Number of Collumns in the Data Set
Cols=DataCols+DateTimeCols;



    
    










