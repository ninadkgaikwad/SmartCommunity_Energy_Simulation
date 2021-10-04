function [ Td ] = HMToDeci( hr,min,sec )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
MinD=min/60;
SecD=sec/3600;
Td=hr+MinD+SecD;


    
%Thm(1,i)=fprintf('The Time is %d hours : %d mins : %d secs\n',hr(1,i),min(1,i),sec(1,i));
    


end

