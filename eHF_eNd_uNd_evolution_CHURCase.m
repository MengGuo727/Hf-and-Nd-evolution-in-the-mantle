%% eHf, e143Nd, and u142Nd evolution with finite mantle mixing
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1. Model magma ocean crystallization
% Physical condition in the mantle:
% We consider mantle crystalize from the bottom, into a lower mantle (LM) 
% and an upper mantle (UM). We consider the upper mantle has a composition 
% of 57% olivine (Ol), 14% Garnet (Gt), and 29% clinophyroxene (Cpx). The  
% lower mantle has a composition of 16% ferropericlase (FP), 79-75% Mg- 
% perovskite (MgPv), and 5-9% Ca-perovskite (CaPv). Eventually, upper 
% mantle will have  a depth of 660 km, and the rest of the mantle 
% (2900km - 660km) will be the lower mantle. The partitioning of Lu, Hf, 
% Sm, and Nd in Ol and FP are negligible.

% Step 1.1 Partition coefficients for mantle minerals
D_Sm_cpx = 0.0775;
D_Sm_Gt = 0.184;
D_Sm_CaPv = 21.1;
D_Sm_MgPv = 0.050;
D_Sm_Majorite = 0.078;

D_Nd_cpx = 0.0445;
D_Nd_Gt = 0.058;
D_Nd_CaPv = 16.7;
D_Nd_MgPv = 0.0161;
D_Nd_Majorite = 0.037;

D_Lu_cpx = 0.171;
D_Lu_Gt = 2.83;
D_Lu_CaPv = 13.5;
D_Lu_MgPv = 1.00;
D_Lu_Majorite = 0.78;

D_Hf_cpx = 0.0577;
D_Hf_Gt = 0.47;
D_Hf_CaPv = 1.70;
D_Hf_MgPv = 1.64;
D_Hf_Majorite = 0.166;

%% Step 1.2 Evolution of Sm/Nd and Lu/Hf ratios and concentration of elements 
% in LM and UM after crystallization of magma ocean

% The mass fraction of mineral in reservoir lower mantle
a_CaPv_LM = [0.05,0.06,0.07,0.08,0.09];
a_MgPv_LMB = [0.79,0.78,0.77,0.76,0.75];

% Calculate the bulk partition coefficient of each element in lower mantle.
[D_Sm_LM] = partition_coefficient(a_CaPv_LM,D_Sm_CaPv,a_MgPv_LMB,D_Sm_MgPv);
[D_Nd_LM] = partition_coefficient(a_CaPv_LM,D_Nd_CaPv,a_MgPv_LMB,D_Nd_MgPv);
[D_Lu_LM] = partition_coefficient(a_CaPv_LM,D_Lu_CaPv,a_MgPv_LMB,D_Lu_MgPv);
[D_Hf_LM] = partition_coefficient(a_CaPv_LM,D_Hf_CaPv,a_MgPv_LMB,D_Hf_MgPv);

% The ratio of Lu/Hf and Sm/Nd in inital BSE (R_liq_0)
R_LuHf_CHUR_t0 = 0.0332;% 0.0332,initial 176Lu/177Hf CHUR, Blichert-Toft et al., 1997
R_SmNd_CHUR_t0 = 0.1960;% 0.1960,initial 147Sm/144Nd CHUR, Bouvier et al., 2008
R_Sm142Nd_CHUR_t0 = 0.000336;% 0.000336,initial 146Sm/144Nd CHUR, Bouvier et al., 2008
C_Hf_BSE_t0 = 5.4286e-08;
C_Nd_BSE_t0 = 7.9138e-07;

% The mass fraction of residual melt (in wt%)
f_UM = 660/2900; % the residual melt became upper mantle after magma ocean crystalization
x = 1-f_UM;

fra_UM = f_UM;
fra_LM = x;

[Rliq_LuHf_LM_eq,Rliq_SmNd_LM_eq,Rliq_Sm142Nd_LM_eq,Rliq_LuHf_LM_fra,...
    Rliq_SmNd_LM_fra,Rliq_Sm142Nd_LM_fra,R_LuHf_LM_eq,R_SmNd_LM_eq,R_Sm142Nd_LM_eq,...
    R_LuHf_LM_fra,R_SmNd_LM_fra,R_Sm142Nd_LM_fra,...
    C_Hf_LM_eq,C_Nd_LM_eq,C_Hf_LM_fra,C_Nd_LM_fra] = ...
    ratio_concentrration_LM(R_LuHf_CHUR_t0,x,D_Lu_LM,D_Hf_LM,R_SmNd_CHUR_t0,...
    R_Sm142Nd_CHUR_t0,D_Sm_LM,D_Nd_LM,C_Hf_BSE_t0,C_Nd_BSE_t0);

% plot LM isotopic ratios
x_LuHf_BSE = [0.1967,0.1967]; x_SmNd_BSE = [0.1,0.24];
y_LuHf_BSE = [0,0.045]; y_SmNd_BSE = [0.0332,0.0332];
figure(1);
plot(R_SmNd_LM_eq(1,:),R_LuHf_LM_eq(1,:),'k-*','LineWidth',2); hold on;
plot(R_SmNd_LM_eq(2,:),R_LuHf_LM_eq(2,:),'c-*','LineWidth',2);
plot(R_SmNd_LM_eq(3,:),R_LuHf_LM_eq(3,:),'m-*','LineWidth',2);
plot(R_SmNd_LM_eq(4,:),R_LuHf_LM_eq(4,:),'g-*','LineWidth',2);
plot(R_SmNd_LM_eq(5,:),R_LuHf_LM_eq(5,:),'r-*','LineWidth',2);
plot(R_SmNd_LM_fra(1,:),R_LuHf_LM_fra(1,:),'k--o','LineWidth',2); 
plot(R_SmNd_LM_fra(2,:),R_LuHf_LM_fra(2,:),'c--o','LineWidth',2);
plot(R_SmNd_LM_fra(3,:),R_LuHf_LM_fra(3,:),'m--o','LineWidth',2);
plot(R_SmNd_LM_fra(4,:),R_LuHf_LM_fra(4,:),'g--o','LineWidth',2);
plot(R_SmNd_LM_fra(5,:),R_LuHf_LM_fra(5,:),'r--o','LineWidth',2);
plot(x_LuHf_BSE,y_LuHf_BSE,'k--','LineWidth',2); plot(x_SmNd_BSE,y_SmNd_BSE,'k--','LineWidth',2);
plot(Rliq_SmNd_LM_eq(1,:),Rliq_LuHf_LM_eq(1,:),'k-*','LineWidth',2);
plot(Rliq_SmNd_LM_eq(2,:),Rliq_LuHf_LM_eq(2,:),'c-*','LineWidth',2);
plot(Rliq_SmNd_LM_eq(3,:),Rliq_LuHf_LM_eq(3,:),'m-*','LineWidth',2);
plot(Rliq_SmNd_LM_eq(4,:),Rliq_LuHf_LM_eq(4,:),'g-*','LineWidth',2);
plot(Rliq_SmNd_LM_eq(5,:),Rliq_LuHf_LM_eq(5,:),'r-*','LineWidth',2);
plot(Rliq_SmNd_LM_fra(1,:),Rliq_LuHf_LM_fra(1,:),'k--o','LineWidth',2); 
plot(Rliq_SmNd_LM_fra(2,:),Rliq_LuHf_LM_fra(2,:),'c--o','LineWidth',2);
plot(Rliq_SmNd_LM_fra(3,:),Rliq_LuHf_LM_fra(3,:),'m--o','LineWidth',2);
plot(Rliq_SmNd_LM_fra(4,:),Rliq_LuHf_LM_fra(4,:),'g--o','LineWidth',2);
plot(Rliq_SmNd_LM_fra(5,:),Rliq_LuHf_LM_fra(5,:),'r--o','LineWidth',2);
xlim([0.1,0.24]); ylim([0.012,0.045]);
% draw Hadean mantle ratio range
xCenter = 0.217;
yCenter = 0.032;
xRadius = 0.005;
yRadius = 0.002;
theta = 0 : 0.01 : 2*pi;
x = xRadius * cos(theta) + xCenter;
y = yRadius * sin(theta) + yCenter;
plot(x, y, 'k','LineWidth', 3);
xlabel('^{147}Sm/^{144}Nd'); ylabel('^{176}Lu/^{177}Hf');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','Location','Northwest')

%% Step 1.3 Calculate the eHf, e143Nd, and u142Nd evolves with time in both LM and UM 
% Set time serie of the Earth's history
tmax = 4.567;% the length of the Earth's history[Ga]
dt = 0.001;% the length of each timestep
t = 0:dt:tmax;% the time serie of the Earth's history
nt = length(t);% the # of timesteps
t = t';

% decay constants
lambda_176Lu = 1.867e-2;% decay constant of 176Lu to 176Hf, in unit Gyr-1
lambda_147Sm = 6.54e-3;% decay constant of 147Sm to 143Nd, in unit Gyr-1
lambda_146Sm = 10.0873;% decay constant of 146Sm to 142Nd, in unit Gyr-1
Na = 6.022E+23; % Avogadro's constant

% the inital condition of the BSE
N_177Hf_BSE_t0 = 5.8060e+41;
N_176Hf_BSE_t0 = 1.6418e+41;
N_176Lu_BSE_t0 = R_LuHf_CHUR_t0 * N_177Hf_BSE_t0;

N_144Nd_BSE_t0 = 5.0212e+42;
N_143Nd_BSE_t0 = 2.5927e+42;
N_142Nd_BSE_t0 = 5.7732e+42;
N_147Sm_BSE_t0 = R_SmNd_CHUR_t0 * N_144Nd_BSE_t0;
N_146Sm_BSE_t0 = R_Sm142Nd_CHUR_t0 * N_144Nd_BSE_t0;

R_Hf_BSE_t0 = N_176Hf_BSE_t0 / N_177Hf_BSE_t0;
R_143Nd_BSE_t0 = N_143Nd_BSE_t0 / N_144Nd_BSE_t0;
R_142Nd_BSE_t0 = N_142Nd_BSE_t0 / N_144Nd_BSE_t0;

% Lu-Hf evoloves with time in the BSE
[N_176Lu_BSE,N_176Hf_BSE,N_177Hf_BSE,R_Hf_BSE] = ...
    BSE_EvolveWithTime(t,nt,N_176Lu_BSE_t0,N_176Hf_BSE_t0,...
    N_177Hf_BSE_t0,R_Hf_BSE_t0,lambda_176Lu);
% 147Sm-143Nd evoloves with time in the BSE
[N_147Sm_BSE,N_143Nd_BSE,N_144Nd_BSE,R_143Nd_BSE] = ...
    BSE_EvolveWithTime(t,nt,N_147Sm_BSE_t0,N_143Nd_BSE_t0,...
    N_144Nd_BSE_t0,R_143Nd_BSE_t0,lambda_147Sm);
% 146Sm-142Nd evoloves with time in the BSE
[N_146Sm_BSE,N_142Nd_BSE,N_144Nd_BSE,R_142Nd_BSE] = ...
    BSE_EvolveWithTime(t,nt,N_146Sm_BSE_t0,N_142Nd_BSE_t0,...
    N_144Nd_BSE_t0,R_142Nd_BSE_t0,lambda_146Sm);

%% Step 1.4 The evolution of eHf, e143Nd, and u142Nd in the BSE
% The evolution of eHf
M_BSE = 4.0359e+24;
M_UM_t0 = M_BSE * f_UM;
M_LM_t0 = M_BSE - M_UM_t0;

M_Hf = 178.49;% atmosic mass of Hf, Blichert-Toft et al. (1997)
N_177Hf_LM_t0 =  (N_177Hf_BSE_t0/(N_177Hf_BSE_t0+N_176Hf_BSE_t0)) * C_Hf_LM_eq(1) .* M_LM_t0 * Na *1e3 / M_Hf;
N_176Hf_LM_t0 =  (N_176Hf_BSE_t0/(N_177Hf_BSE_t0+N_176Hf_BSE_t0)) * C_Hf_LM_eq(1) .* M_LM_t0 * Na *1e3 / M_Hf;
N_176Lu_LM_t0 =  R_LuHf_CHUR_t0 *  N_177Hf_LM_t0; % R_LuHf_LM_eq(1,1) *  N_177Hf_LM_t0; 
N_177Hf_UM_t0 =  N_177Hf_BSE_t0 - N_177Hf_LM_t0;
N_176Hf_UM_t0 =  N_176Hf_BSE_t0 - N_176Hf_LM_t0;
N_176Lu_UM_t0 =  N_176Lu_BSE_t0 - N_176Lu_LM_t0; 

R_Hf_UM_t0 = N_176Hf_UM_t0 / N_177Hf_UM_t0;
R_Hf_LM_t0 = N_176Hf_LM_t0 / N_177Hf_LM_t0;
[N_176Lu_UM,N_176Hf_UM,N_177Hf_UM,R_Hf_UM] = ...
    BSE_EvolveWithTime(t,nt,N_176Lu_UM_t0,N_176Hf_UM_t0,...
    N_177Hf_UM_t0,R_Hf_UM_t0,lambda_176Lu);
[N_176Lu_LM,N_176Hf_LM,N_177Hf_LM,R_Hf_LM] = ...
    BSE_EvolveWithTime(t,nt,N_176Lu_LM_t0,N_176Hf_LM_t0,...
    N_177Hf_LM_t0,R_Hf_LM_t0,lambda_176Lu);
eHf_UM = nan(size(t)); eHf_LM = nan(size(t));
for i = 1:nt
    eHf_UM(i) = (R_Hf_UM(i) - R_Hf_BSE(i))/R_Hf_BSE(i) * 1e4;
    eHf_LM(i) = (R_Hf_LM(i) - R_Hf_BSE(i))/R_Hf_BSE(i) * 1e4;
end
figure(2); subplot(2,2,1);
plot(t,eHf_UM,'r-','linewidth',2); hold on;
plot(t,eHf_LM,'b-.','linewidth',2);
xlabel('Time (Gyr)','Fontsize',15);ylabel('\epsilonHf','Fontsize',15);box on; axis tight;

% The evolution of e143Nd
M_Nd = 144.24; % atmosic mass of Nd
N_144Nd_LM_t0 =  (N_144Nd_BSE_t0/(N_144Nd_BSE_t0+N_143Nd_BSE_t0+N_142Nd_BSE_t0)) * C_Nd_LM_eq(1) .* M_LM_t0 * Na *1e3 / M_Nd;
N_143Nd_LM_t0 =  (N_143Nd_BSE_t0/(N_144Nd_BSE_t0+N_143Nd_BSE_t0+N_142Nd_BSE_t0)) * C_Nd_LM_eq(1) .* M_LM_t0 * Na *1e3 / M_Nd;
N_147Sm_LM_t0 =  R_SmNd_CHUR_t0 *  N_144Nd_LM_t0; % R_SmNd_LM_eq(1,1) *  N_144Nd_LM_t0; 
N_144Nd_UM_t0 =  N_144Nd_BSE_t0 - N_144Nd_LM_t0; 
N_143Nd_UM_t0 =  N_143Nd_BSE_t0 - N_143Nd_LM_t0;
N_147Sm_UM_t0 =  N_147Sm_BSE_t0 - N_147Sm_LM_t0; 

R_143Nd_UM_t0 = N_143Nd_UM_t0 / N_144Nd_UM_t0;
R_143Nd_LM_t0 = N_143Nd_LM_t0 / N_144Nd_LM_t0;
[N_147Sm_UM,N_143Nd_UM,N_144Nd_UM,R_143Nd_UM] = ...
    BSE_EvolveWithTime(t,nt,N_147Sm_UM_t0,N_143Nd_UM_t0,...
    N_144Nd_UM_t0,R_143Nd_UM_t0,lambda_147Sm);
[N_147Sm_LM,N_143Nd_LM,N_144Nd_LM,R_143Nd_LM] = ...
    BSE_EvolveWithTime(t,nt,N_147Sm_LM_t0,N_143Nd_LM_t0,...
    N_144Nd_LM_t0,R_143Nd_LM_t0,lambda_147Sm);
e143Nd_UM = nan(size(t)); e143Nd_LM = nan(size(t));
for i = 1:nt
    e143Nd_UM(i) = (R_143Nd_UM(i) - R_143Nd_BSE(i))/R_143Nd_BSE(i) * 1e4;
    e143Nd_LM(i) = (R_143Nd_LM(i) - R_143Nd_BSE(i))/R_143Nd_BSE(i) * 1e4;
end
figure(2); subplot(2,2,2);
plot(t,e143Nd_UM,'r-','linewidth',2); hold on;
plot(t,e143Nd_LM,'b-.','linewidth',2);
xlabel('Time (Gyr)','Fontsize',15);ylabel('\epsilon^{143}Nd','Fontsize',15); box on; axis tight;

% The evolution of u142Nd
N_142Nd_LM_t0 =  (N_142Nd_BSE_t0/(N_144Nd_BSE_t0+N_143Nd_BSE_t0+N_142Nd_BSE_t0)) * C_Nd_LM_eq(1) .* M_LM_t0 * Na *1e3 / M_Nd;
N_146Sm_LM_t0 =  R_Sm142Nd_CHUR_t0 *  N_144Nd_LM_t0;% R_Sm142Nd_LM_eq(1,1) *  N_144Nd_LM_t0; 
N_142Nd_UM_t0 =  N_142Nd_BSE_t0 - N_142Nd_LM_t0; 
N_146Sm_UM_t0 =  N_146Sm_BSE_t0 - N_146Sm_LM_t0; 

R_142Nd_UM_t0 = N_142Nd_UM_t0 / N_144Nd_UM_t0;
R_142Nd_LM_t0 = N_142Nd_LM_t0 / N_144Nd_LM_t0;
[N_146Sm_UM,N_142Nd_UM,N_144Nd_UM,R_142Nd_UM] = ...
    BSE_EvolveWithTime(t,nt,N_146Sm_UM_t0,N_142Nd_UM_t0,...
    N_144Nd_UM_t0,R_142Nd_UM_t0,lambda_146Sm);
[N_146Sm_LM,N_142Nd_LM,N_144Nd_LM,R_142Nd_LM] = ...
    BSE_EvolveWithTime(t,nt,N_146Sm_LM_t0,N_142Nd_LM_t0,...
    N_144Nd_LM_t0,R_142Nd_LM_t0,lambda_146Sm);
u142Nd_UM = nan(size(t)); u142Nd_LM = nan(size(t));
for i = 1:nt
    u142Nd_UM(i) = (R_142Nd_UM(i) - R_142Nd_BSE(i))/R_142Nd_BSE(i) * 1e6;
    u142Nd_LM(i) = (R_142Nd_LM(i) - R_142Nd_BSE(i))/R_142Nd_BSE(i) * 1e6;
end
figure(2); subplot(2,2,3);
plot(t,u142Nd_UM,'r-','linewidth',2); hold on;
plot(t,u142Nd_LM,'b-.','linewidth',2);
xlabel('Time (Gyr)','Fontsize',15);ylabel('\mu^{142}Nd','Fontsize',15); box on; axis tight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2. Model crustal generation
% Set parameters values for the calcualtion of funtion "CC_growth_fun1.m"
kappa_r = 1;%  decay constant for crustal recycle rate [Gyr-1]
kappa_g = 5;% decay constant for crustal growth rate [Gyr-1]
Mcp = 2.09e22;% mass of present-day CC [kg]
Rs = 5e22;% recycling rates at crustal growth starting time (t=ts) [kg/Gyr],1e22
Rp = 4e21;% recycling rates at present-day time (t=tp) [kg/Gyr],1e20
ts = dt; % the starting time of crustal growth and recycle [Ga]; 0.7670 Ga to dt (3.8Ga vs instant growth)
Krw_factor = 0.5; % 0.1 to 0.8
Krw_s = Rs * Krw_factor;% initial Krw_factor

% Applying CC_growth_fun1
% Calculate the corresponding crustal growth pattern
[M_CC,Mdd,Mud,Krw_first] = CC_growth_fun1(t,ts,tmax,...
    Mcp,kappa_g,Rp,Rs,kappa_r,Krw_s);
% Calcualte formation age distribution
[F,S,m_tp,m,Krw] = ...
    Formation_surface_age_fun(t,Mud,Mdd,M_CC,Krw_first);
    
% Load in the observed formation & surface age distributions
% Formation age distribution data from Korenaga (2018a)
data_formationage = load('korenaga18a_Tunmix_orig.dat');
% Surface age distribution data from Roberts & Spencer (2015)
data_zircon_surf = load('korenaga18a_T_U_Pb.dat');
% use load_FandS_fun function to set data in the same dimension as the time series
[F_Jun_same,S_Jun_same] = load_FandS_fun(t,nt,data_formationage,data_zircon_surf);

% Plot crustal evolution
figure(9);
subplot(2,2,1);plot(t,M_CC,'r','linewidth',2); axis tight;
xlabel('Time (Gyr)','Fontsize',15);ylabel('Net growth of CC','Fontsize',15);title('Net CC growth','Fontsize',15);
subplot(2,2,2);plot(t,Mdd,'r','linewidth',2); axis tight;
xlabel('Time (Gyr)','Fontsize',15);ylabel('Recycling rate','Fontsize',15);title('CC recycling rate','Fontsize',15);
subplot(2,2,3);plot(t,Mud,'r','linewidth',2); axis tight;
xlabel('Time (Gyr)','Fontsize',15);ylabel('Growth rate','Fontsize',15);title('CC generation rate','Fontsize',15);
subplot(2,2,4);plot(t,F,'r','linewidth',2); hold on; plot(t,F_Jun_same,'k','linewidth',2); axis tight;
xlabel('Time (Gyr)','Fontsize',15);ylabel('Formation age','Fontsize',15);title('Formation age distribution','Fontsize',15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3. Model thermal evolution of the mantle
% Step 3.1 Model a normal thermal evolution through the Earth history
% Define the a priori ranges for independent variables of thermal history function
H_BSE_tp_min = 13; H_BSE_tp_max = 19; H_BSE_tp = 16;% heat production in the BSE, in unit TW
H_cc_tp_min = 5; H_cc_tp_max = 10; H_cc_tp = 7.5;% heat production in the continental crust, in unit TW
Q_total_tp_min = 43; Q_total_tp_max = 49; Q_total_tp = 46;% total heat flux, in unit TW
Qc_tp_max = 15; Qc_tp_min = 5; Qc_tp = 10;% core heat flux, in unit TW (Lay et al.,2008)
d_Qc_max = 5; d_Qc_min = 2; d_Qc = 3.5;% difference between the initial and the present-day Qc, in unit TW, O'Rourke et al. (2017)
% constants for performing thermal evolution model
p_K = 2.79e-5;% heat production for K, in unit W/kg
p_U235 = 5.69e-4;% unit heat production for 235U, in unit W/kg
p_U238 = 9.37e-5;% unit heat production for 238U, in unit W/kg
p_Th = 2.69e-5;% unit heat production for Th, in unit W/kg
Th_U = 4;% Th/U ratio
K_U = 1.27e4;% K/U ratio
K40_K = 1.28e-4;% 40K/K ratio
U238_U = 0.9927;% 238U/U ratio
mK_40 = 40;% K's atomic weight
K_factor = K_U*K40_K;% 40K relative to U
heat_factor = U238_U*p_U238 + (1-U238_U)*p_U235 + Th_U*p_Th + K_factor*p_K;
Tkg_atoms = 1e12*1e3*Na/mK_40;% convert unit from Tkg to number of atoms
Ti_tp = 1350; % present-day mantle potential temperature, in unit degree C
V_tp = 5; % present-day plate velocity, in unit cm/yr (Parsons,1981)
dTdP = 1.54e-8; % dT/dp, in unit K/Pa (Korenaga et al., 20020)
type = 2;% set the type to be constant Q scaling law
% Calculate the dependent variables
Q_tp = Q_total_tp - H_cc_tp;
K40_BSE_tp = (K_factor * H_BSE_tp) / heat_factor * Tkg_atoms;
K40_CC_observe = (K_factor * H_cc_tp) / heat_factor * Tkg_atoms;
[Qc,Qc_backward] = Qc_backward_fun(Qc_tp,d_Qc,nt,t);

rhom_MM = 3300; % the average density of mantle, in unit kg/m3

% Calculate Thermal evolution
% Since we are calculating H backward in time, Mc need to be backward in time as well
Mc_backward = flipud(M_CC);
[Ti_backward, Q_backward, H_backward,V_backward, Z_backward] = ...
    Thermal_history(t,type,Q_tp,Qc_backward,Ti_tp,V_tp,...
    rhom_MM,dTdP,Mc_backward,Mcp,H_BSE_tp,H_cc_tp);
% Change the results to be forward in time for later calculation
Ti = flipud(Ti_backward);
Q = flipud(Q_backward);
H = flipud(H_backward);
V = flipud(V_backward);
Z = flipud(Z_backward);

V_OM = V;
% Step 3.2 Rewrite the initial plate velocity to be 10 times higher
% according to the mantle condition after magma ocean solidification
V(1) = V(1) * 10;
for i = 2:500
    V(i) = V(1) - (i-1)*((V(1) - V(500))/(500-1));
end

% load the mantle potential temperature data from Herzberg et al.(2010)
data_Tp = xlsread('Herz data.xlsx');
% set the anchor points for misfit calculation
[Tp_anchorHerz1,Tp_anchorHerz2,Tp_anchorHerz3,Tp_anchorHerz4,...
    t_anchorHerz1,t_anchorHerz2,t_anchorHerz3,t_anchorHerz4,...
    t_Herz,Tp_Herz] = load_Tp_fun(data_Tp);

sz = 60; % scatter size
figure(10);
subplot(2,2,1); hold off;
plot(t,Ti,'k','linewidth',2); hold on;
scatter(t_Herz,Tp_Herz,sz,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','c',...
    'LineWidth',1.5);
xlabel('Time (Gyr)'); ylabel('Mantle Potential Temperature');

subplot(2,2,2); hold off;
plot(t,H,'r','linewidth',2); hold on;
plot(t,Q,'b','linewidth',2);
xlabel('Time (Gyr)');ylabel('H and Q (TW)');

subplot(2,2,3); hold off;
plot(t,V,'linewidth',2);
xlabel('Time (Gyr)');ylabel('Plate velocity');

subplot(2,2,4); hold off;
plot(t,Z,'linewidth',2); 
xlabel('Time (Gyr)');ylabel('Melting depth');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4.Determine the mantle overturn time scale and streching factor according
% to the mantle strain rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % THREE DIFFERENT ALPHA
% V_average = mean(V) * 1e4; % in unit Km/Gyr, V was in unit cm/yr
% t_overturn = round(2*2900/V_average/dt); % # of timestep for LM, RCC, DM, and CR to remelt
% V_KmPerGyr = V * 1e4; % in unit Km/Gyr, V was in unit cm/yr
% strain_rate_normal = 2.5 * 0.08 * V_KmPerGyr / 2900; % in Gyr-1
% alpha = exp(strain_rate_normal * dt); % streching factor, unitless
% 
% factor = ones(size(t));
% factor(1:500) = 0.25;
% factor(2500:nt) = 0.95;
% strain_rate_normal_DM = 2.5 *factor * 0.08 .* V_KmPerGyr / 2900; % in Gyr-1
% alpha_DM = exp(strain_rate_normal_DM * dt); % streching factor, unitless
% 
% V_KmPerGyr_OM = V * 1e4;% in unit Km/Gyr, V was in unit cm/yr
% strain_rate_normal_OM = 2/1.5 * 0.08 *  V_KmPerGyr_OM / 2900; % in Gyr-1
% alpha_OM = exp(strain_rate_normal_OM * dt); % streching factor, unitless
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TWO DIFFERENT ALPHA
V_average = mean(V) * 1e4; % in unit Km/Gyr, V was in unit cm/yr
t_overturn = round(2*2900/V_average/dt); % # of timestep for LM, RCC, DM, and CR to remelt
V_KmPerGyr = V * 1e4; % in unit Km/Gyr, V was in unit cm/yr
strain_rate_normal = 2 * 0.08 * V_KmPerGyr / 2900; % in Gyr-1
alpha = exp(strain_rate_normal * dt); % streching factor, unitless
alpha_OM = alpha; % exp(strain_rate_normal_OM * dt); % streching factor, unitless

factor = ones(size(t));
factor(1:500) = 0.25;
factor(2500:nt) = 0.95;
strain_rate_normal_DM = 2 * factor * 0.08 .* V_KmPerGyr / 2900; % in Gyr-1
alpha_DM = exp(strain_rate_normal_DM * dt); % streching factor, unitless
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 5. Calculate isotopic abundance in CC and MM assuming instant mantle mixing
% and slove for elements' bulk partition coefficients during frist stage mantle melting
% Step 5.1 176Lu-176Hf system
% Step 5.1.1 Instant mixing model
Mmp = 4.015e24; % mass of present-day mantle [kg]
C_Hf_BSE_tp_min = 0.188e-6; C_Hf_BSE_tp_max = 0.266e-6;% 227ppb +- 39 Hf in Tp BSE, Lyubetskay and Korenaga (2007a)
M_177Hf= 177;% atomic mass of 177Hf
Na = 6.022E+23; % Avogadro's constant
Fra_177Hf_tp = 0.18627;% fraction of 177Hf to all Hf isotopes, 18.627% Blichert-Toft et al. (1997)
Fra_all = 1;% fraction of all Hf to all Hf isotopes, 1
M_Lu = 174.967;
C_Hf_CC_tp_max = 4.81e-6; C_Hf_CC_tp_min = 2.59e-6;% 3.7ppm Hf in Tp CC, Rudnick and Gao (2014)
C_Lu_CC_tp_max = 0.39e-6; C_Lu_CC_tp_min = 0.21e-6;% 0.3ppm Lu in Tp CC, Rudnick and Gao (2014)
M_176Lu = 176;% atomic mass of 176Lu
Fra_176Lu_tp = 2.5839e-2;% fraction of 176Lu to all Lu isotopes, 2.5839% Blichert-Toft et al. (1997)
M_176Hf = 176;% atomic mass of 176Hf
Fra_176Hf_tp = 5.2558e-2;% fraction of 176Hf to all Hf isotopes, 5.2558% Blichert-Toft et al. (1997)
% Step 5.1.2 The partition coefficients of Hf and Lu, using euqtion (2) of Hofmann (2003)
F = 0.009;% assumed melting fraction to calcualte first-stage melting's partition coefficient
C_Lu_BSE_tp_max = 63.2e-9; C_Lu_BSE_tp_min = 44.8e-9;% 54 +- 9.2ppb, Lu concentration in the BSE, Lyubetskay and Korenaga (2007a)
D_Hf1_max = (C_Hf_BSE_tp_max - C_Hf_CC_tp_min*F) / (C_Hf_CC_tp_min * (1-F));
D_Hf1_min = (C_Hf_BSE_tp_min - C_Hf_CC_tp_max*F) / (C_Hf_CC_tp_max * (1-F));
D_Lu1_max = (C_Lu_BSE_tp_max - C_Lu_CC_tp_min*F) / (C_Lu_CC_tp_min * (1-F));
D_Lu1_min = (C_Lu_BSE_tp_min - C_Lu_CC_tp_max*F) / (C_Lu_CC_tp_max * (1-F));
D_Lu_inuse = 0.1068; % for rapid growth
% D_Lu_inuse = 0.1168;% for gradual growth
D_Hf_inuse = 0.0684;

% Step 5.2 147Sm-143Nd system
% Step 5.2.1 Instant mixing model
C_Nd_BSE_tp_max= 1.164e-6; C_Nd_BSE_tp_min = 0.824e-6;% 994ppb +- 170 Nd in Tp BSE, Lyubetskay and Korenaga (2007a)
M_144Nd= 144;% atomic mass of 144Nd
Fra_144Nd_tp = 0.23801;%  23.801% Nd are 144Nd at present-day
Fra_all = 1;% fraction of all Nd to all Nd isotopes, 1
C_Nd_CC_tp_max = 23.175e-6; C_Nd_CC_tp_min = 15.75e-6;% 22.5ppm +- 30% Nd in Tp CC, Rudnick and Gao (2014)
C_Sm_CC_tp_max = 6.11e-6; C_Sm_CC_tp_min = 3.29e-6;% 4.7ppm +- 30% Sm in Tp CC, Rudnick and Gao (2014)
M_147Sm = 147;% atomic mass of 147Sm
M_Sm = 150.365;
Fra_147Sm_tp = 0.14995;% 14.995% Sm are 147Sm,
M_143Nd = 143;% atomic mass of 143Nd
Fra_143Nd_tp = 0.12182;% 12.182% Nd are 143Nd at present-day
% Step 5.2.2 calculate Sm in BSE and CC at tp and t0
C_Sm_BSE_tp_max = 379e-9; C_Sm_BSE_tp_min = 269e-9;% 324 +- 55ppb, Sm concentration in the BSE, Lyubetskay and Korenaga (2007a)
D_Sm1_max = (C_Sm_BSE_tp_max- C_Sm_CC_tp_min*F) / (C_Sm_CC_tp_min * (1-F));
D_Sm1_min = (C_Sm_BSE_tp_min - C_Sm_CC_tp_max*F) / (C_Sm_CC_tp_max * (1-F));
D_Nd1_max = (C_Nd_BSE_tp_max - C_Nd_CC_tp_min*F) / (C_Nd_CC_tp_min * (1-F));
D_Nd1_min = (C_Nd_BSE_tp_min - C_Nd_CC_tp_max*F) / (C_Nd_CC_tp_max * (1-F));
D_Sm_inuse = 0.0353;
D_Nd_inuse = 0.0300; % for rapid growth
% D_Nd_inuse = 0.0288; % for gradual growth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 6. Calculate isotopic abundance in CC and MM assuming finite mantle mixing
% Step 6.1 146Sm-142Nd system
l_critical = 10000; % in unit m, can be changed to fit observation
% should be in the order of several km. change the l_0
M_142Nd = 142;
M_146Sm = 146;
F1_tp = 0.10; % melt fraction at present-day
F2 = 0.05; % fraction of second partial melting to generate Mud and residual-crust
rhom_CC = 2900; % mantle density, in unit kg/m3
n_sample = 50;
l_samplingbox = 100000; % in unit m
M_BSE = ones(size(t)) * M_BSE;
n_test = round(t(nt)/0.25); % the # of times of random sampling in the mantle

%% 146Sm-142Nd
[C_Sm142Nd_BSE_max,N_Sm142Nd_dCC_max,N_Sm142Nd_dCC_total_max,N_Sm142Nd_UM_max,N_Sm142Nd_LM_max,N_Sm142Nd_dCC_ini_max,...
    N_Sm142Nd_DM_max,N_Sm142Nd_CR_max,N_Sm142Nd_DM_total_max,N_Sm142Nd_OM_max,N_Sm142Nd_CR_total_max,N_Sm142Nd_BSE_max,...
    N_Sm142Nd_WM_max,N_Sm142Nd_PC_max,N_Sm142Nd_RCC_max,N_Sm142Nd_OM_consumed_max,...
    N_Sm142Nd_RCC_total_max,N_Sm142Nd_PC_ini_max,N_Sm142Nd_RCC_ini_max,N_Sm142Nd_DM_ini_max,N_Sm142Nd_CR_ini_max,...
    M_Reservoir_max,M_LM_max,M_UM_max,l_LM_max,l_UM_max,M_subReservoir_total_max,N_Sm142Nd_OM_consumed_total_max,F1_max,...
    R_142Nd_WM_max,R_142Nd_OM_max,R_142Nd_reservoir_max,e142Nd_WM_max,l_DM_number_max,l_CR_number_max,l_RCC_number_max,...
    l_DM_flag_matrix_max,l_CR_flag_matrix_max,l_RCC_flag_matrix_max,F_consumed_max,e142Nd_sample_max] = ...
    Crust_Extraction_with_Finite_Mantle_Mixing_Cloud(N_144Nd_BSE,N_146Sm_BSE,N_142Nd_BSE,...
    M_BSE,t,nt,dt,Na,M_144Nd,M_142Nd,M_146Sm,Mud,Z,F1_tp,F2,rhom_MM,alpha,alpha_OM,alpha_DM,lambda_146Sm,...
    M_CC,l_critical,Mdd,M_LM_t0,M_UM_t0,N_146Sm_LM_t0,N_146Sm_UM_t0,N_144Nd_LM_t0,...
    N_144Nd_UM_t0,N_142Nd_LM_t0,N_142Nd_UM_t0,t_overturn,D_Sm_inuse,D_Nd_inuse,R_142Nd_BSE,n_sample,l_samplingbox,ts,n_test);

data_u142Nd = xlsread('u142Nd.xlsx');
x_u142Nd_positive = data_u142Nd(:,1);
y_u142Nd_positive = data_u142Nd(:,2);
x_u142Nd_negative = data_u142Nd(:,3);
y_u142Nd_negative = data_u142Nd(:,4);

u142Nd_WM_max = e142Nd_WM_max * 100;
R_142Nd_dCC_max = N_Sm142Nd_dCC_total_max(:,2) ./ N_Sm142Nd_dCC_total_max(:,1);
e142Nd_dCC_max = (R_142Nd_dCC_max - R_142Nd_BSE) ./ R_142Nd_BSE * 1e4;
u142Nd_dCC_max = e142Nd_dCC_max*100;

f = nan(length(e142Nd_sample_max),1);
for i = 1:length(e142Nd_sample_max)
    f(i) = 0.05*rand;
end
t_sample = nan(size(e142Nd_sample_max));

for i = 1:n_test
   t_sample(i,:) = t(i*0.25/dt ).* ones(length(e142Nd_sample_max),1)+f;
end

u142Nd_sample_max = e142Nd_sample_max * 100;

figure(14);hold on;
scatter(x_u142Nd_positive,y_u142Nd_positive,sz,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','c',...
    'LineWidth',1.5);
scatter(x_u142Nd_negative,y_u142Nd_negative,sz,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','c',...
    'LineWidth',1.5);
plot(t,u142Nd_WM_max,'k','linewidth',2);
plot(t,u142Nd_dCC_max,'k-.','lineWidth',2);
scatter(t_sample(1,:),u142Nd_sample_max(1,:),'r+','LineWidth',1);
scatter(t_sample(2,:),u142Nd_sample_max(2,:),'r+','LineWidth',1);
scatter(t_sample(3,:),u142Nd_sample_max(3,:),'r+','LineWidth',1);
scatter(t_sample(4,:),u142Nd_sample_max(4,:),'r+','LineWidth',1);
scatter(t_sample(5,:),u142Nd_sample_max(5,:),'r+','LineWidth',1);
scatter(t_sample(6,:),u142Nd_sample_max(6,:),'r+','LineWidth',1);
scatter(t_sample(7,:),u142Nd_sample_max(7,:),'r+','LineWidth',1);
scatter(t_sample(8,:),u142Nd_sample_max(8,:),'r+','LineWidth',1);
scatter(t_sample(9,:),u142Nd_sample_max(9,:),'r+','LineWidth',1);
scatter(t_sample(10,:),u142Nd_sample_max(10,:),'r+','LineWidth',1);
scatter(t_sample(11,:),u142Nd_sample_max(11,:),'r+','LineWidth',1);
scatter(t_sample(12,:),u142Nd_sample_max(12,:),'r+','LineWidth',1);
scatter(t_sample(13,:),u142Nd_sample_max(13,:),'r+','LineWidth',1);
scatter(t_sample(14,:),u142Nd_sample_max(14,:),'r+','LineWidth',1);
scatter(t_sample(15,:),u142Nd_sample_max(15,:),'r+','LineWidth',1);
scatter(t_sample(16,:),u142Nd_sample_max(16,:),'r+','LineWidth',1);
scatter(t_sample(17,:),u142Nd_sample_max(17,:),'r+','LineWidth',1);
scatter(t_sample(18,:),u142Nd_sample_max(18,:),'r+','LineWidth',1);
y = zeros(size(t));
plot(t,y,'k:','linewidth',2)
xlabel('Time (Gyr)','Fontsize',12);ylabel('\mu^{142}Nd','Fontsize',12); box on;
axis tight;

%% Lu-Hf system
[C_LuHf_BSE_max,N_LuHearlef_dCC_max,N_LuHf_dCC_total_max,N_LuHf_UM_max,N_LuHf_LM_max,N_LuHf_dCC_ini_max,...
    N_LuHf_DM_max,N_LuHf_CR_max,N_LuHf_DM_total_max,N_LuHf_OM_max,N_LuHf_CR_total_max,N_LuHf_BSE_max,...
    N_LuHf_WM_max,N_LuHf_PC_max,N_LuHf_RCC_max,N_LuHf_OM_consumed_max,...
    N_LuHf_RCC_total_max,N_LuHf_PC_ini_max,N_LuHf_RCC_ini_max,N_LuHf_DM_ini_max,N_LuHf_CR_ini_max,...
    M_Reservoir_max,M_LM_max,M_UM_max,l_LM_max,l_UM_max,M_subReservoir_total_max,N_LuHf_OM_consumed_total_max,F1_max,...
    R_Hf_WM_max,R_Hf_OM_max,R_Hf_reservoir_max,e176Hf_WM_max,l_DM_number_max,l_CR_number_max,l_RCC_number_max,...
    l_DM_flag_matrix_max,l_CR_flag_matrix_max_max,l_RCC_flag_matrix_max,F_consumed_max,eHf_sample_max] = ...
Crust_Extraction_with_Finite_Mantle_Mixing_Cloud(N_177Hf_BSE,N_176Lu_BSE,N_176Hf_BSE,...
    M_BSE,t,nt,dt,Na,M_177Hf,M_176Hf,M_176Lu,Mud,Z,F1_tp,F2,rhom_MM,alpha,alpha_OM,alpha_DM,lambda_176Lu,...
    M_CC,l_critical,Mdd,M_LM_t0,M_UM_t0,N_176Lu_LM_t0,N_176Lu_UM_t0,N_177Hf_LM_t0,...
    N_177Hf_UM_t0,N_176Hf_LM_t0,N_176Hf_UM_t0,t_overturn,D_Lu_inuse,D_Hf_inuse,R_Hf_BSE,n_sample,l_samplingbox,ts,n_test);

data_e167Hf = xlsread('e167Hf.xlsx');
x_e167Hf_positive = data_e167Hf(:,1); x_e167Hf_positive = t(nt)- x_e167Hf_positive;
y_e167Hf_positive = data_e167Hf(:,2); 

data_Kemp10 = xlsread('Kemp10.xlsx');
x_eHf_Kemp10 = data_Kemp10(:,1) / 1000; x_eHf_Kemp10 = t(nt)- x_eHf_Kemp10;
y_eHf_Kemp10 = data_Kemp10(:,2); err = data_Kemp10(:,3);

R_Hf_dCC_max = N_LuHf_dCC_total_max(:,2) ./ N_LuHf_dCC_total_max(:,1);
e167Hf_dCC_max = (R_Hf_dCC_max - R_Hf_BSE) ./ R_Hf_BSE * 1e4;

x_38 = [3.8 0];x_38 = t(nt) - x_38;y_38 = [0 15.2];
x_45 = [0 t(nt)];y_45 = [0 15.2];
figure(16); hold on;
scatter(x_e167Hf_positive,y_e167Hf_positive,sz,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','c',...
    'LineWidth',1.5);
scatter(t_sample(1,:),eHf_sample_max(1,:),'r+','LineWidth',1);
scatter(t_sample(2,:),eHf_sample_max(2,:),'r+','LineWidth',1);
scatter(t_sample(3,:),eHf_sample_max(3,:),'r+','LineWidth',1);
scatter(t_sample(4,:),eHf_sample_max(4,:),'r+','LineWidth',1);
scatter(t_sample(5,:),eHf_sample_max(5,:),'r+','LineWidth',1);
scatter(t_sample(6,:),eHf_sample_max(6,:),'r+','LineWidth',1);
scatter(t_sample(7,:),eHf_sample_max(7,:),'r+','LineWidth',1);
scatter(t_sample(8,:),eHf_sample_max(8,:),'r+','LineWidth',1);
scatter(t_sample(9,:),eHf_sample_max(9,:),'r+','LineWidth',1);
scatter(t_sample(10,:),eHf_sample_max(10,:),'r+','LineWidth',1);
scatter(t_sample(11,:),eHf_sample_max(11,:),'r+','LineWidth',1);
scatter(t_sample(12,:),eHf_sample_max(12,:),'r+','LineWidth',1);
scatter(t_sample(13,:),eHf_sample_max(13,:),'r+','LineWidth',1);
scatter(t_sample(14,:),eHf_sample_max(14,:),'r+','LineWidth',1);
scatter(t_sample(15,:),eHf_sample_max(15,:),'r+','LineWidth',1);
scatter(t_sample(16,:),eHf_sample_max(16,:),'r+','LineWidth',1);
scatter(t_sample(17,:),eHf_sample_max(17,:),'r+','LineWidth',1);
scatter(t_sample(18,:),eHf_sample_max(18,:),'r+','LineWidth',1);
y = zeros(size(t));
plot(t,e176Hf_WM_max,'k','linewidth',2); 
plot(t,e167Hf_dCC_max,'k-.','linewidth',2);
plot(x_38,y_38,'k-','linewidth',2);
plot(x_45,y_45,'k-','linewidth',2);
plot(t,y,'k:','linewidth',2); 
xlabel('Time (Gyr)','Fontsize',16);ylabel('\epsilonHf','Fontsize',16); box on; axis tight;

%% 147Sm-143Nd system
[C_Sm143Nd_BSE_max,N_Sm143Nd_dCC_max,N_Sm143Nd_dCC_total_max,N_Sm143Nd_UM_max,N_Sm143Nd_LM_max,N_Sm143Nd_dCC_ini_max,...
    N_Sm143Nd_DM_max,N_Sm143Nd_CR_max,N_Sm143Nd_DM_total_max,N_Sm143Nd_OM_max,N_Sm143Nd_CR_total_max,N_Sm143Nd_BSE_max,...
    N_Sm143Nd_WM_max,N_Sm143Nd_PC_max,N_Sm143Nd_RCC_max,N_Sm143Nd_OM_consumed_max,...
    N_Sm143Nd_RCC_total_max,N_Sm143Nd_PC_ini_max,N_Sm143Nd_RCC_ini_max,N_Sm143Nd_DM_ini_max,N_Sm143Nd_CR_ini_max,...
    M_Reservoir,M_LM,M_UM,l_LM,l_UM,M_subReservoir_total,N_Sm143Nd_OM_consumed_total_max,F1,...
    R_143Nd_WM_max,R_143Nd_OM_max,R_143Nd_reservoir_max,e143Nd_WM_max,l_DM_number,l_CR_number,l_RCC_number,...
    l_DM_flag_matrix,l_CR_flag_matrix,l_RCC_flag_matrix,F_consumed,e143Nd_sample_max] = ...
Crust_Extraction_with_Finite_Mantle_Mixing_Cloud(N_144Nd_BSE,N_147Sm_BSE,N_143Nd_BSE,...
    M_BSE,t,nt,dt,Na,M_144Nd,M_143Nd,M_147Sm,Mud,Z,F1_tp,F2,rhom_MM,alpha,alpha_OM,alpha_DM,lambda_147Sm,...
    M_CC,l_critical,Mdd,M_LM_t0,M_UM_t0,N_147Sm_LM_t0,N_147Sm_UM_t0,N_144Nd_LM_t0,...
    N_144Nd_UM_t0,N_143Nd_LM_t0,N_143Nd_UM_t0,t_overturn,D_Sm_inuse,D_Nd_inuse,R_143Nd_BSE,n_sample,l_samplingbox,ts,n_test);

data_e143Nd = xlsread('e143Nd.xlsx');
x_e143Nd = data_e143Nd(:,1);
y_e143Nd = data_e143Nd(:,2);

R_143Nd_dCC_max = N_Sm143Nd_dCC_total_max(:,2) ./ N_Sm143Nd_dCC_total_max(:,1);
e143Nd_dCC_max = (R_143Nd_dCC_max - R_143Nd_BSE) ./ R_143Nd_BSE * 1e4;

x_38 = [3.8 0];x_38 = t(nt)-x_38;y_38 = [0, 9];
x_45 = [0 t(nt)];y_45 = [0, 9];
figure(20);hold on;
scatter(x_e143Nd,y_e143Nd,sz,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','c',...
    'LineWidth',1.5); 
scatter(t_sample(1,:),e143Nd_sample_max(1,:),'r+','LineWidth',1);
scatter(t_sample(2,:),e143Nd_sample_max(2,:),'r+','LineWidth',1);
scatter(t_sample(3,:),e143Nd_sample_max(3,:),'r+','LineWidth',1);
scatter(t_sample(4,:),e143Nd_sample_max(4,:),'r+','LineWidth',1);
scatter(t_sample(5,:),e143Nd_sample_max(5,:),'r+','LineWidth',1);
scatter(t_sample(6,:),e143Nd_sample_max(6,:),'r+','LineWidth',1);
scatter(t_sample(7,:),e143Nd_sample_max(7,:),'r+','LineWidth',1);
scatter(t_sample(8,:),e143Nd_sample_max(8,:),'r+','LineWidth',1);
scatter(t_sample(9,:),e143Nd_sample_max(9,:),'r+','LineWidth',1);
scatter(t_sample(10,:),e143Nd_sample_max(10,:),'r+','LineWidth',1);
scatter(t_sample(11,:),e143Nd_sample_max(11,:),'r+','LineWidth',1);
scatter(t_sample(12,:),e143Nd_sample_max(12,:),'r+','LineWidth',1);
scatter(t_sample(13,:),e143Nd_sample_max(13,:),'r+','LineWidth',1);
scatter(t_sample(14,:),e143Nd_sample_max(14,:),'r+','LineWidth',1);
scatter(t_sample(15,:),e143Nd_sample_max(15,:),'r+','LineWidth',1);
scatter(t_sample(16,:),e143Nd_sample_max(16,:),'r+','LineWidth',1);
scatter(t_sample(17,:),e143Nd_sample_max(17,:),'r+','LineWidth',1);
scatter(t_sample(18,:),e143Nd_sample_max(18,:),'r+','LineWidth',1);
y = zeros(size(t));
plot(t,e143Nd_WM_max,'k','linewidth',2);
plot(t,e143Nd_dCC_max,'k-.','linewidth',2);
plot(x_38,y_38,'k-','linewidth',2);
plot(x_45,y_45,'k-','linewidth',2);
plot(t,y,'k:','linewidth',2); box on; 
xlabel('Time (Gyr)','Fontsize',16);ylabel('\epsilon^{143}Nd','Fontsize',16)
axis tight;
