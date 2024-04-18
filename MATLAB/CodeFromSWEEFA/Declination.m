function [ dec ] = Declination( n )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
len=length(n);

for i=1:len
    
    dec(1,i)=23.45*(sin((360/365)*(n(1,i)-81)*(pi/180)));
    
end

end

