function [ CosInciAngle ] = ViewFactor(beta,phi,tilt,phic)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

CosInciAngle=(cos((pi/180)*(beta))*cos((pi/180)*(phi-phic))*sin((pi/180)*(tilt)))+(sin((pi/180)*(beta))*cos((pi/180)*(tilt))); % Incidence Angle


end

