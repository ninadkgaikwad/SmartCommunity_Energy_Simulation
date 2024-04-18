function [NewResFile] = PecanStreet_Low2HighRes(OriginalResolution,NewResolution,ProcessingType,WeatherFileLoad_Full)


%% Weather Data - Low Resolution to High Resolution - Zero Order Hold

[R,C] = size(WeatherFileLoad_Full);

FileRes=OriginalResolution;

%% Creating New Resolution File

if (ProcessingType==1) % Full File to be Processed
    
    % Getting Start and End DateTime
    StartDay=WeatherFileLoad_Full(1,1);
    StartMonth=WeatherFileLoad_Full(1,2);
    StartYear=WeatherFileLoad_Full(1,3);
    StartTime=WeatherFileLoad_Full(1,4);

    [ StartHr,StartMin,StartSec ] = DeciToHM( StartTime ); % Decimal Time to HMS

    EndDay=WeatherFileLoad_Full(R,1);
    EndMonth=WeatherFileLoad_Full(R,2);
    EndYear=WeatherFileLoad_Full(R,3);
    EndTime=WeatherFileLoad_Full(R,4);

    [ EndHr,EndMin,EndSec ] = DeciToHM( EndTime ); % Decimal Time to HMS   
    
    
elseif (ProcessingType==2) % Part of the File to be Processed we take User Input
    
    % Getting Start and End DateTime
    StartDay=1;
    StartMonth=1;
    StartYear=2017;
    StartTime=0;

    [ StartHr,StartMin,StartSec ] = DeciToHM( StartTime ); % Decimal Time to HMS

    EndDay=31;
    EndMonth=12;
    EndYear=2017;
    EndTime=23.5;

    [ EndHr,EndMin,EndSec ] = DeciToHM( EndTime ); % Decimal Time to HMS    
    
end

[ ~,StartIndex,EndIndex ] = DateTimeSeriesSlicer(WeatherFileLoad_Full,1,OriginalResolution,StartYear,EndYear,StartMonth,EndMonth,StartDay,EndDay,0,EndTime);

% Creating DateTime Object
Start_DateTime=datetime(StartYear,StartMonth,StartDay,StartHr,StartMin,StartSec);
End_DateTime=datetime(EndYear,EndMonth,EndDay,EndHr,EndMin,EndSec);

% Creating New Resolution Duration
NewResolution_Duration=duration(0,NewResolution,0);

% While Loop for New Resolution File
NewResFile=zeros(1,C); % Initializing5

DateTimeArray=Start_DateTime; % Initializing

Counter_NewTime=1; % Initilizing
Counter_OldTime=StartIndex; % Initializing

while (DateTimeArray(Counter_NewTime)<=End_DateTime)
    
 
    if (Counter_OldTime==1)
        
        % Decomposing New DateTime
        Day=DateTimeArray(Counter_NewTime).Day;
        Month=DateTimeArray(Counter_NewTime).Month;
        Year=DateTimeArray(Counter_NewTime).Year;
        
        Hour=DateTimeArray(Counter_NewTime).Hour;
        Minute=DateTimeArray(Counter_NewTime).Minute;
        Second=DateTimeArray(Counter_NewTime).Second;
        
        % Converting HMS to Decimal Time
        [ Time ] = HMToDeci( Hour,Minute,Second );
        
        % Updating NewRes File with OldRes File
        NewResFile(Counter_NewTime,:)=horzcat([Day,Month,Year,Time],WeatherFileLoad_Full(Counter_OldTime,5:end));
        
        % Incrementing Counter
        Counter_OldTime=Counter_OldTime+1 ;          
        
    else
        
        % Getting Current and Previous OldRes File DateTime
        CurrentDay=WeatherFileLoad_Full(Counter_OldTime,1);
        CurrentMonth=WeatherFileLoad_Full(Counter_OldTime,2);
        CurrentYear=WeatherFileLoad_Full(Counter_OldTime,3);
        
        CurrentTime=WeatherFileLoad_Full(Counter_OldTime,4);
        [ CurrentHr,CurrentMin,CurrentSec ] = DeciToHM( CurrentTime ); % Converting Decimal Time to HMS
        
        PreviousDay=WeatherFileLoad_Full(Counter_OldTime-1,1);
        PreviousMonth=WeatherFileLoad_Full(Counter_OldTime-1,2);
        PreviousYear=WeatherFileLoad_Full(Counter_OldTime-1,3);
        
        PreviousTime=WeatherFileLoad_Full(Counter_OldTime-1,4);
        [ PreviousHr,PreviousMin,PreviousSec ] = DeciToHM( PreviousTime ); % Converting Decimal Time to HMS      
         
        % Converting To DateTime Object
        CurrentDateTime_OldResFile=datetime(CurrentYear,CurrentMonth,CurrentDay,CurrentHr,CurrentMin,CurrentSec);
        PreviousDateTime_OldResFile=datetime(PreviousYear,PreviousMonth,PreviousDay,PreviousHr,PreviousMin,PreviousSec);
        
        if ((DateTimeArray(Counter_NewTime)<=CurrentDateTime_OldResFile)&&(DateTimeArray(Counter_NewTime)>PreviousDateTime_OldResFile))
            
            % Decomposing New DateTime
            Day=DateTimeArray(Counter_NewTime).Day;
            Month=DateTimeArray(Counter_NewTime).Month;
            Year=DateTimeArray(Counter_NewTime).Year;

            Hour=DateTimeArray(Counter_NewTime).Hour;
            Minute=DateTimeArray(Counter_NewTime).Minute;
            Second=DateTimeArray(Counter_NewTime).Second;

            % Converting HMS to Decimal Time
            [ Time ] = HMToDeci( Hour,Minute,Second );

            % Updating NewRes File with OldRes File
            NewResFile(Counter_NewTime,:)=horzcat([Day,Month,Year,Time],WeatherFileLoad_Full(Counter_OldTime,5:end));            

        else
            
            % Incrementing Counter
            Counter_OldTime=Counter_OldTime+1  ;           

            % Decomposing New DateTime
            Day=DateTimeArray(Counter_NewTime).Day;
            Month=DateTimeArray(Counter_NewTime).Month;
            Year=DateTimeArray(Counter_NewTime).Year;

            Hour=DateTimeArray(Counter_NewTime).Hour;
            Minute=DateTimeArray(Counter_NewTime).Minute;
            Second=DateTimeArray(Counter_NewTime).Second;

            % Converting HMS to Decimal Time
            [ Time ] = HMToDeci( Hour,Minute,Second );

            % Updating NewRes File with OldRes File
            NewResFile(Counter_NewTime,:)=horzcat([Day,Month,Year,Time],WeatherFileLoad_Full(Counter_OldTime,5:end));        
                             
            
        end
        

        
    end
    
    % Creating Next New Resolution DateTime
    DateTimeArray(Counter_NewTime+1)=DateTimeArray(Counter_NewTime)+NewResolution_Duration;

    % Incrementing Counter
    Counter_NewTime=Counter_NewTime+1  ;   
    
end



end

