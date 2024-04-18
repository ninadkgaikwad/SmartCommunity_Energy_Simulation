function [HEMSHouseRCModel_Output] = HEMS_HouseRCModel(HEMSHouse_Input,HEMSWeatherData_Output,Simulation_Params,HEMSPlant_Params,HEMSHouse_States)

% Author: Ninad Kiran Gaikwad
% Date: Jun/10/2019
% Description: House RC Model - Simulates the House RC Thermal Model

%% House RC Model - Simulates the House RC Thermal Model

%% Extracting Required Data from the Structs

%------------------------- HEMSHouse_Input--------------------------------%
R_w=HEMSHouse_Input.R_w;
R_attic=HEMSHouse_Input.R_attic;
R_roof=HEMSHouse_Input.R_roof;
R_im=HEMSHouse_Input.R_im;
R_win=HEMSHouse_Input.R_win;
C_w=HEMSHouse_Input.C_w;
C_attic=HEMSHouse_Input.C_attic;
C_im=HEMSHouse_Input.C_im;
C_in=HEMSHouse_Input.C_in;
C1=HEMSHouse_Input.C1;
C2=HEMSHouse_Input.C2;
C3=HEMSHouse_Input.C3;
Human_Num=HEMSHouse_Input.Human_Num;
Human_Heat=HEMSHouse_Input.Human_Heat;
Appliance_Heat=HEMSHouse_Input.Appliance_Heat;
Q_ac=HEMSHouse_Input.Q_ac;
Cp=HEMSHouse_Input.Cp;
V=HEMSHouse_Input.V;
Den_Air=HEMSHouse_Input.Den_Air;
C_oew=HEMSHouse_Input.C_oew;
SHGC=HEMSHouse_Input.SHGC;
Alpha_w=HEMSHouse_Input.Alpha_w;
Alpha_r=HEMSHouse_Input.Alpha_r;
Area_w=HEMSHouse_Input.Area_w;
Tilt_w=HEMSHouse_Input.Tilt_w;
Azi_w=HEMSHouse_Input.Azi_w;
Area_r=HEMSHouse_Input.Area_r;
Tilt_r=HEMSHouse_Input.Tilt_r;
Azi_r=HEMSHouse_Input.Azi_r;
Area_win=HEMSHouse_Input.Area_win;
Tilt_win=HEMSHouse_Input.Tilt_win;
Azi_win=HEMSHouse_Input.Azi_win;

%------------------------- Simulation_Params------------------------------%
StepSize=Simulation_Params.StepSize;

%---------------------- HEMSWeatherData_Output----------------------------%
Ws=HEMSWeatherData_Output.Ws;
T_am=HEMSWeatherData_Output.T_am;
DNI=HEMSWeatherData_Output.DNI;
DateTime_Matrix=HEMSWeatherData_Output.DateTime_Matrix;


%------------------------- Simulation_Params------------------------------%
hem=HEMSPlant_Params.hem;
Lat=HEMSPlant_Params.Lat;
Long=HEMSPlant_Params.Long;
Ltm=HEMSPlant_Params.Ltm;

%-------------------------- HEMSHouse_States------------------------------%
T_wall1=HEMSHouse_States.T_wall1;
T_ave1=HEMSHouse_States.T_ave1;
T_attic1=HEMSHouse_States.T_attic1;
T_im1=HEMSHouse_States.T_im1; 

u_k_hvac=HEMSHouse_States.u_k_hvac;

%% Basic Computations

Q_ihl=Human_Num*Human_Heat+Appliance_Heat; % [W] Includes Humans and all internal appliances which contribute to the Heat gain/loss of the house

WeatherData_Length=length(Ws);

% initialization of House Zone Temperatures
T_wall(1)=T_wall1;
T_ave(1)=T_ave1;
T_attic(1)=T_attic1;
T_im(1)=T_im1;

%% House Simulation

for ii=1:WeatherData_Length                  
    
    % Computing Q_venti and Q_infil
    Q_venti(ii) = Cp*V*Den_Air*(T_am(ii) - T_ave(ii));    
    Q_infil(ii) = Cp*(T_am(ii) - T_ave(ii))*C_oew*Ws(ii);

    % Computing View Factor - Window Roof Wall
    Day=DateTime_Matrix(ii,1);
    Month=DateTime_Matrix(ii,2);
    Year=DateTime_Matrix(ii,3);
    Time=DateTime_Matrix(ii,4); 

    [n]=JulianDay(Day,Month,Year); % Computing Julian Day
    [ dec ] = Declination( n ); % Computing Declination
    [ ST] = ClockToSolarTime( n,hem,Ltm,Long,Time); % Converting ClockTime to SolarTime
    if (ST<0) % Correcting For Negative SolarTime        
        ST=ST+24;        
    elseif (ST>=24) % Correcting For greater than 24 SolarTime 
        ST=ST-24;        
    end
    [ Ha ] = HourAngle( ST ); % Computing HourAngle
    [ beta, phi] = AltiAzi( dec ,Lat, Ha); % Computing Elevation and Azimuth angle of the Sun

    % View Factor Wall
    F_w(ii)=0;
    for jj=1:length(Area_w)       
        [ VF ] = ViewFactor(beta,phi,Tilt_w(jj),Azi_w(jj));        
        if ((VF<0) || (beta<0)) % In Shadow           
            VF=0;            
        end        
         F_w(ii)= F_w(ii)+(Area_w(jj)*VF);        
    end
    F_w(ii)=F_w(ii)/sum(Area_w);

    % View Factor Roof
    F_r(ii)=0;
    for jj=1:length(Area_r)        
        [ VF ] = ViewFactor(beta,phi,Tilt_r(jj),Azi_r(jj));        
        if ((VF<0) || (beta<0)) % In Shadow           variation of constants control system
            VF=0;            
        end        
         F_r(ii)= F_r(ii)+(Area_r(jj)*VF);        
    end
    F_r(ii)=F_r(ii)/sum(Area_r); 

    % View Factor Window
    F_win(ii)=0;
    for jj=1:length(Area_win)        
        [ VF ] = ViewFactor(beta,phi,Tilt_win(jj),Azi_win(jj));        
        if ((VF<0) || (beta<0)) % In Shadow           
            VF=0;            
        end        
         F_win(ii)= F_win(ii)+(Area_win(jj)*VF);        
    end
    F_win(ii)=F_win(ii)/sum(Area_win);    

    % Computing Convective Heat Transfer Coefficient
    h_c(ii) = 11.4 + (5.7*Ws(ii));

    % Computing Q_solarvariation of constants control system
    Q_solar(ii) = F_win(ii)*DNI(ii)*sum(Area_win)*SHGC;

    % Computing T_sol_w and T_sol_r
    T_sol_w(ii) = (Alpha_w/h_c(ii))*F_w(ii)*DNI(ii) + T_am(ii);    
    T_sol_r(ii) = (Alpha_r/h_c(ii))*F_r(ii)*DNI(ii) + T_am(ii); 

%     % Computing Euler Derivatives
%     D_T_wall = (((T_sol_w(ii) - T_wall(ii))*(1/(R_w/2))) - ((T_wall(ii) - T_ave(ii))*(1/(R_w/2))))*(1/C_w);
%     D_T_ave = (((T_wall(ii) - T_ave(ii))*(1/(R_w/2))) + ((T_attic(ii) - T_ave(ii))*(1/(R_attic))) + ((T_im(ii) - T_ave(ii))*(1/(R_im))) + ((T_am(ii) - T_ave(ii))*(1/(R_win))) + (C1*Q_ihl) + (C2*Q_ac) + (Q_venti(ii)) + (Q_infil(ii)))*(1/C_in); 
%     D_T_attic = (((T_sol_r(ii) - T_attic(ii))*(1/(R_roof))) - ((T_attic(ii) - T_ave(ii))*(1/(R_attic))))*(1/C_attic);   
%     D_T_im = (-((T_im(ii) - T_ave(ii))*(1/(R_im))) + (C3*Q_solar(ii)))*(1/C_im);
%     %D_T_im = (-(T_im(ii) - T_ave(ii))*(1/(R_im))) *(1/C_im); %;
% 
%     % Computing States
%     T_wall(ii+1) = T_wall(ii) + StepSize*(D_T_wall);
%     T_ave(ii+1) = T_ave(ii) + StepSize*(D_T_ave);
%     T_attic(ii+1) = T_attic(ii) + StepSize*(D_T_attic);
%     T_im(ii+1) = T_im(ii) + StepSize*(D_T_im); 
% 
%     Ta_Room(ii)=T_ave(ii);
% 
%     Ta_Outside(ii)=T_am(ii);
    
    % Creatig Current State Matrix
    X=[T_wall(ii);T_ave(ii);T_attic(ii);T_im(ii)];
    
    % Creating Input Matrix    
    U = [T_am(ii) ; T_sol_w(ii) ; T_sol_r(ii) ; Q_ihl ; -u_k_hvac*Q_ac ; Q_venti(ii) ; Q_infil(ii) ; Q_solar(ii)];
    
    % Creating R_2 and R_1
    R_2 = -(R_attic*R_im*R_win)-((R_w/2)*R_im*R_win)-((R_w/2)*R_attic*R_win)-((R_w/2)*R_attic*R_im);
    R_1 = (R_attic*R_im*R_win*(R_w/2));
    
    % Computing Ac Bc -- A B    
    Ac = [(-4/(C_w*R_w))       (2/(C_w*R_w))                         0                                     0         ;....
          (2/(C_in*R_w))     (R_2/(C_in*R_1))                 (1/(C_in*R_attic))                    (1/(C_in*R_im))  ; ....
          0                (1/(C_attic*R_attic))   ((-R_attic-R_roof)/(C_attic*R_roof*R_attic))            0         ;....
          0                   (1/(C_im*R_im))                        0                              (-1/(C_im*R_im))];
    
    Bc = [0                     (1/(C_w*(R_w/2)))        0                        0             0           0           0       0         ;....
          (1/(C_in*(R_win)))            0                0                    (C1/C_in)    (C2/C_in)    (1/C_in)    (1/C_in)    0         ;...
          0                             0             (1/(C_attic*R_roof))        0             0           0           0       0         ;...
          0                             0                0                        0             0           0           0       (C3/C_im)];
    
    A = expm(Ac*StepSize);

    B = Ac\(expm(Ac*StepSize)-eye(4))*(Bc);  
    
    % Computing Next State    
    X_new = A*X + B*U;   
    
    % Computing States
    T_wall(ii+1) = X_new(1,1);
    T_ave(ii+1) = X_new(2,1);
    T_attic(ii+1) = X_new(3,1);
    T_im(ii+1) = X_new(4,1);   
    
    Ta_Room(ii)=T_ave(ii+1);

    Ta_Outside(ii)=T_am(ii);    
                
end

%% Creating HEMSHouseRCModel_Output - Struct

HEMSHouseRCModel_Output=[]; % Empty Struct

HEMSHouseRCModel_Output.T_sol_w=T_sol_w;
HEMSHouseRCModel_Output.T_sol_r=T_sol_r;
HEMSHouseRCModel_Output.T_wall=T_wall;
HEMSHouseRCModel_Output.T_ave=T_ave;
HEMSHouseRCModel_Output.T_attic=T_attic;
HEMSHouseRCModel_Output.T_im=T_im;
HEMSHouseRCModel_Output.Ta_Room=Ta_Room;
HEMSHouseRCModel_Output.Ta_Outside=Ta_Outside;

end

