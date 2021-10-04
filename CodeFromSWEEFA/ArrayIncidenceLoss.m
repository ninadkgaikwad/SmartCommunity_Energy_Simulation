function [ Iciam,Icsf ] = ArrayIncidenceLoss( Ic,CosInciAngle,bo,SF )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Fiam=1-(bo*((1/CosInciAngle)-1)); % Calculating Incidence Angle Modifier

Iciam=Ic*Fiam; % Solar Power Available after Incidence Angle Modification

Icsf=Iciam*(1-(SF/100)); % Solar Power Available For PV Module after Soiling Effect
end

