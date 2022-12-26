 function [Ti_backward, Q_backward, H_backward,V_backward, Z_backward] ...
     = Thermal_history(t,type,Q_tp,Qc_backward,Ti_tp,V_tp,rhom,...
     dTdP,Mc_backward,Mcp,H_BSE_tp,H_cc_tp)
% This function models the thermal evolution of the Earth.
% 
% Meng Guo, Yale University
% Summer 2019

% Input variables and their meanings:
% t: time, from 0 Ga to 4.567 Ga;
% type: the evolution type, 1 is for conventional scaling law (exponential 
% decrease), and 2 is for constant Q scaling law;
% Q_tp: present-day surface heat loss, in unit TW
% Ti_tp: present-day average internal temperature, in unit celsius degree 
% V_tp: present-day plate velocity, in unit cm/yr
% rhom: the average density of mantle, in unit kg/m3
% dTdP: adiabatic gradient in the shallow mantle, in unit K/Pa
% Mc: continental crust mass at each time-step, in unit kg
% Mcp: present-day continental crust mass, in unit kg
% H_BSE_tp: present-day bulk silicate earth heat production, in unit TW
% H_cc_tp: present-day continental crust heat production, in unit TW
% Output variables and their meanings: 
% H: internal heat production in the convectin mantle, in unit TW
% Q: surface heat loss, in unit TW
% Ti: average internal temperrature, in unit celsius degree
% Ur: Urey ratio Ur = H/Q, unitless
% V: plate velocity, in unit cm/yr
% Z: the initial depth of mantle melting, in unit meter

H_backward = nan(size(t)); 
Q_backward = nan(size(t)); 
Ti_backward = nan(size(t));
Ur_backward = nan(size(t));
V_backward = nan(size(t));
Z_backward = nan(size(t)); 
    
dt = t(2)-t(1);
nt = length(t);

% parameters to define a classical heat flow scaling law
a = 1.4013e-19;%  a is determined to yield Q of 36 TW at Ti of 1350 celsius degree
B = 6.52;% B determines the sensitivity of Q with respect to a change in Ti;
% when activation energy is 300kj/mol,B = 6.52 

%% Calculate H(t) vs t, using equation H(t)= H(tp)* sum[ci*pi*exp(lamada_i*t)/sum(ci* pi)]
% The heat production of convecting mantle H(t) = H_m(t) + (1 - Mc/Mcp)* H_cc(t)
H_m_tp = H_BSE_tp - H_cc_tp; % present-day mantle heat production
H_backward(1) = H_m_tp;% present-day mantle heat production

Q_Haden = Q_tp; % Q from 0 to 0.5Ga, in unit TW

% define decay constant for each heat-producing elements
lamada_K = 0.555; % in unit 1/Gyr 
lamada_U235 = 0.985; % in unit 1/Gyr
lamada_U238 = 0.155; % in unit 1/Gyr
lamada_Th = 0.0495; % in unit 1/Gyr

% define heat generation rate for each heat-producing elements
p_K = 2.79e-5; % in unit W/kg
p_U235 = 5.69e-4; % in unit W/kg
p_U238 = 9.37e-5; % in unit W/kg
p_Th = 2.69e-5; % in unit W/kg

% define concentration for each heat-producing elements U:Th:K = 1:4:(1.27
% × 10^4) for both continental crust (CC) and bulk silicate Earth (BSE)
Th_U = 4; % Th/U ratio
K_U = 1.27e4; % K/U ratio
K40_K = 1.28e-4; % 40K/K ratio
U238_U = 0.9927;% 238U/U ratio

c_U = 1; % what value used for c_U doesn't matter,
% because heat production will be adjusted to satisfy
% the assumed Urey ratio.

c_K = K_U * K40_K * c_U; % 40K concentration
c_U235 = (1-U238_U) * c_U; % 235U concentration
c_U238 = U238_U * c_U; % 238U concentration
c_Th = Th_U * c_U; % Th concentration, all Th are 232Th

tmp_sum = c_K * p_K + c_U235 * p_U235 + c_U238 * p_U238 + c_Th * p_Th;

for i = 2:nt
  frac_K(i) = c_K * p_K * exp(lamada_K * t(i))/ tmp_sum;
  frac_U235(i) = c_U235 * p_U235 * exp(lamada_U235 * t(i))/ tmp_sum;
  frac_U238(i) = c_U238 * p_U238 * exp(lamada_U238 * t(i))/ tmp_sum;
  frac_Th(i) = c_Th * p_Th * exp(lamada_Th * t(i))/ tmp_sum;
  H_backward(i) = (H_m_tp + (1 - Mc_backward(i)/Mcp)*H_cc_tp ) * ...
      (frac_K(i) + frac_U235(i) + frac_U238(i) + frac_Th(i));
end

%% Conventional scaling and Korenaga scaling (constant Q), calculate Q and Ti
Cm = 4.97e27; % heat capacity of the whole mantle, in unit J/K [Stacey, 1981]
% C = 7e27; % heat capacity of the whole Earth, in unit J/K [Stacey, 1981]
Ti_backward(1) = Ti_tp; % in unit °C, the present-day mantle potential temperature
Q_backward(1) = Q_tp; % present-day heat flux
    
t_haden = 501;% Hadean ends at 0.5Ga

factor1 = dt*1e9*365*24*60*60/Cm*1e12;

for i = 2:t_haden  
    Ti_backward(i) = Ti_backward(i-1)- factor1 * (H_backward(i-1)-Q_backward(i-1)+ Qc_backward(i-1));
    if type == 1
        Q_backward(i) = a * (Ti_backward(i)^B);
    else
        Q_backward(i) = Q_tp;
    end
end
    
for i = (t_haden+1):nt  
    Ti_backward(i) = Ti_backward(i-1)- factor1 * (H_backward(i-1)-Q_backward(i-1)+ Qc_backward(i-1));
    if type == 1
        Q_backward(i) = a * (Ti_backward(i)^B);
    else
        Q_backward(i) = Q_Haden;
    end
end

%% calculate Urey Ratio (Ur) 
for i = 1:nt;
  Ur_backward(i) = H_backward(i) / Q_backward(i);
end

%% Calculate the plate velocity (V) for Ar degassing model
V_backward(1) = V_tp;
for i = 2:nt;
  V_backward(i) = V_backward(1) * ((Q_backward(i)/Q_backward(1)) * (Ti_backward(1)/Ti_backward(i)))^2;% in unit cm/yr 
end

%% Calculate the initial depth of mantle melting (Z) for Ar degassing model
g = 9.8; % the gravitational acceleration , in unit m/s2
for i = 1:nt;
  Z_backward(i) = (Ti_backward(i) - 1150)/(rhom * g * (1.2e-7 - dTdP)); % in unit m
end
