function [ Ib, Id,C] = BeamDiffGHI(n, GHI,beta )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

C=0.095+(0.04*sin((pi/180)*(360/365)*(n-100)));

Ib=GHI/(sin((pi/180)*beta)+C);

Id=C*Ib;



end

