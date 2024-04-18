function [n]=JulianDay(Day,Month,Year)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

% Initializing Day
StartMonthStartDay=0;

% Computing Leap Year or Not
[ LeapYear ] = LeapYearFinder( Year );


if LeapYear==0
   
   switch Month
       
       case 1
           StartMonthStartDay=1;
           
       case 2
           StartMonthStartDay=32;
           
       case 3
           StartMonthStartDay=60;
           
       case 4
           StartMonthStartDay=91;
           
       case 5
           StartMonthStartDay=121;
           
       case 6
           StartMonthStartDay=152;
           
       case 7
           StartMonthStartDay=182;
           
       case 8
           StartMonthStartDay=213;
           
       case 9
           StartMonthStartDay=244;
           
       case 10
           StartMonthStartDay=274;
           
       case 11
           StartMonthStartDay=305;
           
       case 12
           StartMonthStartDay=335;
           
          
   end
    
   

    
elseif LeapYear==1
   
   switch Month
       
       case 1
           StartMonthStartDay=1;
           
       case 2
           StartMonthStartDay=32;
           
       case 3
           StartMonthStartDay=61;
           
       case 4
           StartMonthStartDay=92;
           
       case 5
           StartMonthStartDay=122;
           
       case 6
           StartMonthStartDay=153;
           
       case 7
           StartMonthStartDay=183;
           
       case 8
           StartMonthStartDay=214;
           
       case 9
           StartMonthStartDay=245;
           
       case 10
           StartMonthStartDay=275;
           
       case 11
           StartMonthStartDay=306;
           
       case 12
           StartMonthStartDay=336;
           
          
   end
    
   
end

n=StartMonthStartDay+Day-1;



end

