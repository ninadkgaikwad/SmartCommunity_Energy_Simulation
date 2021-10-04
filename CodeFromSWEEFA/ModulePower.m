function [ Pmodtot,Pmodin,Tm]=ModulePower( Pmod,ModTemCF,ModNum,Tn,Gn,Icsf,T,Ic,Uo,U1,Ws )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Values of U_{0} varied from 23.5 to 26.5 with a combined fit = 25 W/m^{2}K
% Values of U_{1} varied from 6.25 to 7.68 with a combined fit = 6.84 W/m^{2}K

Tm=T+((Ic)/(Uo+(U1*Ws))); % Faiman's Module Temperature Model

Pmodin=Pmod*(1+((ModTemCF/100)*(Tm-Tn)))*(Icsf/Gn); % Power generated in one Module

Pmodtot=ModNum*Pmodin; % Power generated in all Modules

end

