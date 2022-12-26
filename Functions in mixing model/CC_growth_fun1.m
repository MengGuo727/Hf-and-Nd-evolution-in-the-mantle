function [Mc,Mdd,Mud,Krw] ...
    = CC_growth_fun1(t,ts,tmax,Mcp,kappa_g,Rp,Rs,kappa_r,Krw_s)
% this function simulates the crustal growth patterns and returns the net
% growth model, crustal generation, recycling, and reworking rates
%
% Meng Guo, Yale University
% Summer, 2019

% create empty matrixes to store results
Mc= nan(size(t));% net crustal growth
Mdd= nan(size(t));% crustal recycling rate
Mud= nan(size(t));% crustal generation rate 
Krw = nan(size(t));% crustal reworking rate

nt=length(t);

% present-day crustal reworking rate
Krw_p = Krw_s * Rp/Rs;

for i = 1:nt
  if t(i) < ts
    Mc(i) = 0;
    Mdd(i) = 0;
    Mud(i) = 0;
    Krw(i) = 0;
  else
    Mc(i) = (Mcp/(1-exp(-kappa_g*(tmax-ts))))*(1-exp(-kappa_g*(t(i)-ts)));
    Mdd(i) = Rs + ((Rp-Rs)/(1-exp(-kappa_r*(tmax-ts))))*(1-exp(-kappa_r*(t(i)-ts)));
    Mud(i) = (Mcp/(1-exp(-kappa_g*(tmax-ts))))*(kappa_g * ...
						exp(-kappa_g*(t(i)-ts)))+ Mdd(i);
    Krw(i) = Krw_s + ((Krw_p-Krw_s)/(1-exp(-kappa_r*(tmax-ts))))*...
        (1-exp(-kappa_r*(t(i)-ts)));
  end
end


