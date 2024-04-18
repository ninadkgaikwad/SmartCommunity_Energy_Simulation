function [ beta, phi] = AltiAzi( dec ,L, H)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
%len1=length(dec);

len2=length(H);


    %for i=1:len1
        
        for j=1:len2
      
            %e=(cos(L(1,i)*(pi/180))*cos(dec(1,i)*(pi/180))*cos(H(1,i)*(pi/180)))+(sin(L(1,i)*(pi/180))*sin(dec(1,i)*(pi/180)))    
        
             beta(1,j)=(180/pi)*(asin((cos(L*(pi/180))*cos(dec*(pi/180))*cos(H(1,j)*(pi/180)))+(sin(L*(pi/180))*sin(dec*(pi/180)))));
            
                       
             azi1(1,j)=(180/pi)*(asin((cos(dec*(pi/180))*sin(H(1,j)*(pi/180)))/(cos(beta(1,j)*(pi/180)))));
    
             
             
             azi11(1,j)=abs(azi1(1,j));
    
             azi2(1,j)=180-azi11(1,j);
    
             azi22(1,j)=abs(azi2(1,j));
    
             x=cos(H(1,j)*(pi/180));
    
             y=(tan(dec*(pi/180)))/(tan(L*(pi/180)));
    
       
            if x>=y
                
                phi(1,j)=azi1(1,j);
                
            else
                
                if azi1(1,j)>=0
                    
                    phi(1,j)=azi2(1,j);
                    
                else
                    
                    phi(1,j)=-(azi2(1,j));
        
                end
                
            end
        
            
        
               
            
               end
               
               
            %end
        end
        
    
    
    




