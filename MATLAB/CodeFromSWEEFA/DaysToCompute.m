function [ StartDay, EndDay ] = DaysToCompute( LeapYear,StartDay,StartMonth,EndDay,EndMonth )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
StartMonthStartDay=0;

EndMonthStartDay=0;


if LeapYear==0
   
   switch StartMonth
       
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
    
   switch EndMonth
       
       case 1
           EndMonthStartDay=1;
           
       case 2
           EndMonthStartDay=32;
           
       case 3
           EndMonthStartDay=60;
           
       case 4
           EndMonthStartDay=91;
           
       case 5
           EndMonthStartDay=121;
           
       case 6
           EndMonthStartDay=152;
           
       case 7
           EndMonthStartDay=182;
           
       case 8
           EndMonthStartDay=213;
           
       case 9
           EndMonthStartDay=244;
           
       case 10
           EndMonthStartDay=274;
           
       case 11
           EndMonthStartDay=305;
           
       case 12
           EndMonthStartDay=335;
         
              
   end

    
elseif LeapYear==1
   
   switch StartMonth
       
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
    
   switch EndMonth
       
       case 1
           EndMonthStartDay=1;
           
       case 2
           EndMonthStartDay=32;
           
       case 3
           EndMonthStartDay=61;
           
       case 4
           EndMonthStartDay=92;
           
       case 5
           EndMonthStartDay=122;
           
       case 6
           EndMonthStartDay=153;
           
       case 7
           EndMonthStartDay=183;
           
       case 8
           EndMonthStartDay=214;
           
       case 9
           EndMonthStartDay=245;
           
       case 10
           EndMonthStartDay=275;
           
       case 11
           EndMonthStartDay=306;
           
       case 12
           EndMonthStartDay=336;
         
              
    end 
end

StartDay=StartMonthStartDay+StartDay-1;

EndDay=EndMonthStartDay+EndDay-1;

