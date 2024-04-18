function [ OriginalSeries,StartIndex,EndIndex ] = DateTimeSeriesSlicer(OriginalDataSeries,SeriesNum3,Res,StartYear,EndYear,StartMonth,EndMonth,StartDay,EndDay,StartTime,EndTime)

% Getting Size of the OriginalDataSeries

[r,c]=size(OriginalDataSeries);

% Creating Day Time Vector based on the File Resolution

DayVector=0:(Res/60):(24-(Res/60));

DayLen=length(DayVector);

% Finding Time Value Indices within the Day Vector

DiffStartTime=abs(StartTime-DayVector);
[MinST,IndexST]=min(DiffStartTime);

DiffEndTime=abs(EndTime-DayVector);
[MinET,IndexET]=min(DiffEndTime);

% Finding the Start Index

for i=1:DayLen:r
    
    if ((OriginalDataSeries(i,1)==StartDay)&&(OriginalDataSeries(i,2)==StartMonth)&&(OriginalDataSeries(i,3)==StartYear))
        
        StartIndex=i+(IndexST-1);
        
    end

end

% Finding the End Index

for i=1:DayLen:r
    
    if ((OriginalDataSeries(i,1)==EndDay)&&(OriginalDataSeries(i,2)==EndMonth)&&(OriginalDataSeries(i,3)==EndYear))
        
        EndIndex=i+(IndexET-1);
        
    end

end

% Getting the OriginalSeries

OriginalSeries=OriginalDataSeries(StartIndex:EndIndex,(4+SeriesNum3));

