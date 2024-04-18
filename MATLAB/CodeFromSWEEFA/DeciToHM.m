function [ hr,min,sec ] = DeciToHM( Td )
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here
len=length(Td);

for i=1:len
    
    hr(1,i)=fix(Td(1,i)/1);
    
    mmm(1,i)=rem(Td(1,i),1);
    
    mm(1,i)=mmm(1,i)*60;
    
    min(1,i)=fix(mm(1,i)/1);
    
    sss(1,i)=rem(mm(1,i),1);
    
    ss(1,i)=sss(1,i)*60;
    
    sec(1,i)=fix(ss(1,i)/1);
    
%Thm(1,i)=fprintf('The Time is %d hours : %d mins : %d secs\n',hr(1,i),min(1,i),sec(1,i));
    
end


end

