function [ DateTimeMatrix,TotDataPoints,Time ] = StartEndCalender( StartYear,StartMonth,StartDay,TotDays,Res,DataCols )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% Initializing Start: Day, Month and Year
SD=StartDay;
SM=StartMonth;
SY=StartYear;

% Computing Nuber of Data Points in One day and Total number of Days: DayPoints and TotDataPoints & Time Vector
DayPoints=24*(60/Res);
TotDataPoints=TotDays*DayPoints;
Time=0:Res/60:(24-(Res/60));
Time=Time';

% Computing Total Number of Columns (Initial 4 Columns are for the Date time Signature)
TotCols=DataCols+4;

% Initializing DateTimeMatrix to Zeros and Count
DateTimeMatrix=zeros(TotDataPoints,TotCols);
Count=0;

% Intializing Month/Year Markers: Tn, Th, Tl and Yr 
Tn=0;
Th=0;
Tl=0;
Yr=0;

% Creating DateTimeMatrix using FOR LOOP
for i=1:TotDays
    
    % Computing Current value of Start Year: SY
    if i==1
        
        SY=SY;
        SM=SM;
        SD=SD;
        
    elseif (SD==31)&&(SM==12)
        
        SY=SY+1;
        SM=1;
        SD=1;
       
        % End of Year Marker
            Yr=1;
         

    end  
            % Finding if current Year is a Leap Year using a Pre-defined Function
        LP=LeapYearFinder(SY);

        % Resetting Start Day and Start Month: SD, SM 1
        if SD==31

            SD=1;
            SM=SM+1; 
            
            % End of Month Marker
            Tn=1;

        elseif SD==30

            for j=[4,6,9,11]

                if SM==j

                    SD=1;
                    SM=SM+1; 
                    
                    % End of Month Marker
                    Th=1;                    

                    break;

                end
            end

        elseif ((SD==28)&&(LP==0)&&(SM==2))||((SD==29)&&(LP==1)&&(SM==2))

            SD=1;
            SM=SM+1; 

            % End of Month Marker
            Tl=1;            
            
        end

DayIncrement=(i~=1)&&(((SD~=31)&&(SM==1||SM==3||SM==5||SM==7||SM==8||SM==10||SM==12))||((SD~=30)&&(SM==4||SM==6||SM==9||SM==11))||(((SD~=28)&&(LP==0)&&(SM==2))||((SD~=29)&&(LP==1)&&(SM==2))));

        if DayIncrement
            
            if (Tn==1)||(Th==1)||(Tl==1)||(Yr==1)
                
                SD=1;
                
            else
            
                SD=SD+1;
            
            end
        
        end
    
    % Updating DateMatrix for each Day
%         DateMatrix(i,:)=[SD,SM,SY];  

    for k=1:length(Time)
        
        Count=Count+1;
        DateTimeMatrix(Count,:)=[SD,SM,SY,Time(k,1),zeros(1,DataCols)];
        
    end
        
    % Resetting End of Month/Year Markers: Tn, Th and Tl 
    Tn=0;
    Th=0;
    Tl=0;
    Yr=0;
        
end 
    
    
        
end