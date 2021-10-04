function [ PVout,INVpin,INVpout,Pgrid,ArrayMismatchLoss,ShadingLoss,LIDLoss,OhmicLossP,InverterLoss,TransformerLossP,TrackerLossP  ] = PVoutputPower( Pmodtot, LID,LS,Arraymismat,Crys,Shading, OhmicLoss,TrackerL,INVeff,TransLoss  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% LID = Light Induced Degradation (for Crystalline Modules) (1-3%; Default=2%)
% LS = Light Soaking (For Thin Flim Modules)(3-5% or more; Default=3%)
% Arraymismat = Array Mismatch Factor (Default=2%)
% Crys = Is the PV Technology Crystalline (Crys==1) or Thin Film (Crys==0)??
% Shading = Shading Loss Factor (Default=1%)
% OhmicLoss = Array wiring loss (Default=3%)
% INVeff = Inverter Efficiency (%; Given in Inverter Datasheet)
% TransLoss = Transformer Loss Factor (Default=1%)

% Calculating Power Output from Array
if Crys==1
    PVout=Pmodtot*(1-((LID+Arraymismat+Shading)/100));
else
    PVout=Pmodtot*(1-((Arraymismat+Shading)/100)+((LS)/100));
end

% Calculating Power Input to Inverter
INVpin=PVout*(1-(OhmicLoss/100));

% Calculating Power Output from Inverter (%INVpout=INVpin*(1+(INVTemCF*(T-Tn)));)
INVpout=INVpin*(INVeff/100);

% Calculating Power Output from Inverter (%INVpout=INVpin*(1+(INVTemCF*(T-Tn)));)
TrackerLossPP=INVpout*(1-(TrackerL/100));

% Calculating Power Output to Grid through Transformer
Pgrid=TrackerLossPP*(1-(TransLoss/100));

% Calculating Power Losses 
ArrayMismatchLoss=Pmodtot*(Arraymismat/100);
ShadingLoss=Pmodtot*(Shading/100);
LIDLoss=Pmodtot*(LID/100);
OhmicLossP=PVout*(OhmicLoss/100);
InverterLoss=INVpin*(1-(INVeff/100));
TrackerLossP=INVpout*(TrackerL/100);
TransformerLossP=INVpout*(TransLoss/100);

end

