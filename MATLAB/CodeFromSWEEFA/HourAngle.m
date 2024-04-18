function [ H ] = HourAngle( Hp )
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here
hpsz=size(Hp);

len1=hpsz(1,1);
len2=hpsz(1,2);

for i=1:len1
    for j=1:len2
    
        H(i,j)=15*(12-(Hp(1,j)));
        
%         H(i,j)=(360/23.9344696)*(12-(Hp(1,j)));
    
    end
end


end

