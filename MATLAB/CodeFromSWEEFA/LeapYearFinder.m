function [ Leap ] = LeapYearFinder( Year )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% If Leap Year, then Leap=1;
%If Leap Year, the Leap=0;

for i=1:length(Year)

a=rem(Year(1,i),4); 
b=rem(Year(1,i),100); 
c=rem(Year(1,i),400);

if ((a==0)&&(b~=0))||(c==0)
    
    Leap(1,i)=1; 
    
else
    
    Leap(1,i)=0;
    
end

end

end

