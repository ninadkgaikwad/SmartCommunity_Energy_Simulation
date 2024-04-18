function [ SunRise,SunSet,Indicator ] = SunRiseSet(L,dec )
%UNTITLED19 Summary of this function goes here
%   Detailed explanation goes here

%Output in Solar Time-Decimal Hours
len=length(dec);

Indicator=zeros(1,len);

for i=1:len
    
    Hsr=(180/pi)*(acos((-(tan(L*(pi/180)))*tan(dec(1,i)*(pi/180)))));    
    
    hsr=abs(Hsr);
    
    Q=((3.467)/(cos((pi/180)*L)*cos((pi/180)*dec(1,i))*sin((pi/180)*Hsr)));
    
    SunRise(1,i)=12-(hsr/15)-(Q/60);
    
    SunSet(1,i)=12+(hsr/15)+(Q/60);
    
%         SunRise(1,i)=12-(hsr/(360/23.9344696))-(Q/60);
%     
%     SunSet(1,i)=12+(hsr/(360/23.9344696))+(Q/60);
    
    Sr=abs(imag(SunRise(1,i)));
    
    Ss=abs(imag(SunSet(1,i)));
    
    if (Sr>0)&&(Ss>0)
        
        if (L>0)&&(dec(1,i)>=0)
            
            Indicator(1,i)=1; % 24 Hour Sunlight
        end
        
        if (L>0)&&(dec(1,i)<=0)
            
            Indicator(1,i)=-1; % 24 Hour Night
            
        end
        
        if (L<0)&&(dec(1,i)<=0)
            
            Indicator(1,i)=1; % 24 Hour Sunlight
            
        end
        
        if (L<0)&&(dec(1,i)>=0)
            
            Indicator(1,i)=-1; % 24 Hour Night
            
        end
    else
        
        Indicator(1,i)=0; % Sunrise and Sunset are present
        
    end
end


end

