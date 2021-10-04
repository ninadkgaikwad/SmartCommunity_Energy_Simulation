function [ Ic,Ibc,Idc,Irc,CosInciAngle ] = FixedTilt( Ib,Id,C,beta,phi,tilt,phic,rho)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

CosInciAngle=(cos((pi/180)*(beta))*cos((pi/180)*(phi-phic))*sin((pi/180)*(tilt)))+(sin((pi/180)*(beta))*cos((pi/180)*(tilt))); % Incidence Angle

Ibc=Ib*CosInciAngle; % Beam Component on Collector

Idc=Id*((1+cos((pi/180)*(tilt)))/(2)); % Diffused Component on Collector

Irc=rho*Ib*(sin((pi/180)*(beta))+C)*((1-cos((pi/180)*(tilt)))/(2)); % Reflected Component on the Collector

Ic=Ibc+Idc+Irc; % Total Solar Irradiance on the collector

end

