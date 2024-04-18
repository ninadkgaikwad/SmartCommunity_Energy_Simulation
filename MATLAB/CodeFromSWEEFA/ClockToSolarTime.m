function [ ST,B,E] = ClockToSolarTime( n,hem,Ltm,L,CT)

%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

% Output is in Decimal Hours
len1=length(n);
len2=length(CT);

for i=1:len1
    
    B(1,i)=(360/364)*(n(1,i)-81);
    
    E(1,i)=(9.87*sin((pi/180)*(2*B(1,i))))-(7.53*cos((pi/180)*B(1,i)))-(1.5*sin((pi/180)*B(1,i)));
    
    for j=1:len2
    
    ST(1,j)=(CT(1,j))-((hem*(Ltm-L)*4)/60)+(E(1,i)/60);
    
    end
    
end


end




