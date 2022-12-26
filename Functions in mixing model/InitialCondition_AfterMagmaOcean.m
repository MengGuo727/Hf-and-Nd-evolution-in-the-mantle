%% Simulate the inital solidificaiton condition of mantle right after magma ocean
clear all;
close all;
clc;

%% Physical condition in the mantle
% We consider mantle crystalize from the bottom, into a lower mantle and 
% an upper mantle. For a small amount of mantle to crystalize (height of h), 
% it becomes 60% solid, and 40% liquid, but the liquid will be trapped in 
% the crystallized matrix and eventually solidify to be part of the lower 
% mantle. We consider the upper mantle has a composition of 57% olivine (Ol),
% 14% Garnet (Gt), and 29% clinophyroxene (Cpx). The lower mantle has a 
% composition of 16% ferropericlase (FP), 79-75% Mg-perovskite (MgPv), and
% 5-9% Ca-perovskite (CaPv). Eventually, upper mantle will have a depth of
% 660 km, and the rest of the mantle (2900km - 660km) will be the lower
% mantle.

% UM: upper mantle
% LM: lower mantle

%% Step 1. Calcualte the Lu, Hf, Sm, and Nd's partition coefficients in lower mantle
% Partition coefficients for different mantle minerals (Table S2 of Caro et al., 2005)
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

% The mass fraction of mineral in reservoir lower mantle
a_CaPv_LM = [0.05,0.06,0.07,0.08,0.09];
a_FP_LM = 0.16;
a_MgPv_LMB = [0.79,0.78,0.77,0.76,0.75];

% Calculate the bulk partition coefficient of each element in lower mantle
[D_Sm_LM] = partition_coefficient(a_CaPv_LM,D_Sm_CaPv,a_MgPv_LMB,D_Sm_MgPv);
[D_Nd_LM] = partition_coefficient(a_CaPv_LM,D_Nd_CaPv,a_MgPv_LMB,D_Nd_MgPv);
[D_Lu_LM] = partition_coefficient(a_CaPv_LM,D_Lu_CaPv,a_MgPv_LMB,D_Lu_MgPv);
[D_Hf_LM] = partition_coefficient(a_CaPv_LM,D_Hf_CaPv,a_MgPv_LMB,D_Hf_MgPv);

%% Step 2. Calcualte BSE inital daughter isotope's concentration
% BSE initial Hf concentration
Na = 6.022E+23; % Avogadro's constant
N_Hf177_BSE_tp = 5.80599494648837e+41;
m_177Hf = 177;
M_177Hf_BSE_t0 = N_Hf177_BSE_tp * m_177Hf / 1000 / Na;
M_BSE = 4.0359e+24;
C_177Hf_BSE_t0 = M_177Hf_BSE_t0 / M_BSE;
% BSE initial Nd concentration
Nd144_BSE_tp = 0.238 * 1.250e-6 * M_BSE * Na *1e3 /144;%in unit # of atoms
m_144Nd = 144;
M_144Nd_BSE_t0 = Nd144_BSE_tp * m_144Nd / 1000 / Na;
C_144Nd_BSE_t0 = M_144Nd_BSE_t0 / M_BSE;

%% Step 3. Calcualte the isotopic ratio during lower mantle crystallization
n_crystalize = 1000; % number of crystallization event during magma ocean solidification
h_mantle = 2900; % depth of the mantle, in unit km
h_UM = 660; % depth of upper mantle, in unit km
h = (h_mantle - h_UM) / n_crystalize; % depth each crystalized LM material, including trapped melt
R_LuHf_CHUR_t0 = 0.0332;% 0.0332,initial 176Lu/177Hf CHUR, Blichert-Toft et al., 1997
R_SmNd_CHUR_t0 = 0.1960;% 0.1960,initial 147Sm/144Nd CHUR, Bouvier et al., 2008

type_crystallization = 2; % 1 is equilibrium crystallization, 2 is fractional crystallization
f_crystal = 0.6; % fraction of solid during crystalize each h amount of mantle
[R_liqid_LuHf_f1,R_solid_LuHf_f1,R_liqid_SmNd_f1,R_solid_SmNd_f1,...
    R_average_LuHf_LM_f1,R_average_SmNd_LM_f1,R_average_LuHf_UM_f1,R_average_SmNd_UM_f1,...
    C_Nd_solid_f1, C_Nd_liquid_f1, C_Hf_solid_f1, C_Hf_liquid_f1,x_f1,f_f1] = ...
    LM_crystallization(n_crystalize,f_crystal,R_LuHf_CHUR_t0,R_SmNd_CHUR_t0,...
    type_crystallization,D_Sm_LM,D_Nd_LM,D_Lu_LM,D_Hf_LM,h,h_mantle,...
    C_144Nd_BSE_t0,C_177Hf_BSE_t0,1);

f_crystal = 1; % fraction of solid during crystalize each h amount of mantle
[R_liqid_LuHf_f2,R_solid_LuHf_f2,R_liqid_SmNd_f2,R_solid_SmNd_f2,...
    R_average_LuHf_LM_f2,R_average_SmNd_LM_f2,R_average_LuHf_UM_f2,R_average_SmNd_UM_f2,...
    C_Nd_solid_f2, C_Nd_liquid_f2, C_Hf_solid_f2, C_Hf_liquid_f2,x_f2,f_f2] = ...
    LM_crystallization(n_crystalize,f_crystal,R_LuHf_CHUR_t0,R_SmNd_CHUR_t0,...
    type_crystallization,D_Sm_LM,D_Nd_LM,D_Lu_LM,D_Hf_LM,h,h_mantle,...
    C_144Nd_BSE_t0,C_177Hf_BSE_t0,1);

type_crystallization = 1; % 1 is equilibrium crystallization, 2 is fractional crystallization
f_crystal = 1; % fraction of solid during crystalize each h amount of mantle
[R_liqid_LuHf_eq,R_solid_LuHf_eq,R_liqid_SmNd_eq,R_solid_SmNd_eq,...
    R_average_LuHf_LM_eq,R_average_SmNd_LM_eq,R_average_LuHf_UM_eq,R_average_SmNd_UM_eq,...
    C_Nd_solid_eq, C_Nd_liquid_eq, C_Hf_solid_eq, C_Hf_liquid_eq,x_eq,f_eq] = ...
    LM_crystallization(n_crystalize,f_crystal,R_LuHf_CHUR_t0,R_SmNd_CHUR_t0,...
    type_crystallization,D_Sm_LM,D_Nd_LM,D_Lu_LM,D_Hf_LM,h,h_mantle,...
    C_144Nd_BSE_t0,C_177Hf_BSE_t0,1);

% completely fractional crystallization
type_crystallization = 2; % 1 is equilibrium crystallization, 2 is fractional crystallization
f_crystal = 1; % fraction of solid during crystalize each h amount of mantle
[R_liqid_LuHf_fra,R_solid_LuHf_fra,R_liqid_SmNd_fra,R_solid_SmNd_fra,...
    R_average_LuHf_LM_fra,R_average_SmNd_LM_fra,R_average_LuHf_UM_fra,R_average_SmNd_UM_fra,...
    C_Nd_solid_fra, C_Nd_liquid_fra, C_Hf_solid_fra, C_Hf_liquid_fra,x_fra,f_fra] = ...
    LM_crystallization(n_crystalize,f_crystal,R_LuHf_CHUR_t0,R_SmNd_CHUR_t0,...
    type_crystallization,D_Sm_LM,D_Nd_LM,D_Lu_LM,D_Hf_LM,h,h_mantle,...
    C_144Nd_BSE_t0,C_177Hf_BSE_t0,2);

%% Plot the results
h_LM = nan(size(R_average_SmNd_UM_eq(:,5)));
h_LM(1) = 0;
length_result = length(R_average_SmNd_UM_eq(:,5));
for i = 2:length_result 
    h_LM(i) = h_LM(i-1) + h;
end

figure(1);subplot(2,2,1);grid on;
plot(h_LM,R_average_LuHf_LM_f1(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_LM_f1(:,2),'LineWidth',2);
plot(h_LM,R_average_LuHf_LM_f1(:,3),'LineWidth',2);
plot(h_LM,R_average_LuHf_LM_f1(:,4),'LineWidth',2);
plot(h_LM,R_average_LuHf_LM_f1(:,5),'LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{176}Lu/^{177}Hf in crystalized lower mantle');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','Location','Northwest')
subplot(2,2,2);grid on;
plot(h_LM,R_average_LuHf_UM_f1(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_UM_f1(:,2),'LineWidth',2);
plot(h_LM,R_average_LuHf_UM_f1(:,3),'LineWidth',2);
plot(h_LM,R_average_LuHf_UM_f1(:,4),'LineWidth',2);
plot(h_LM,R_average_LuHf_UM_f1(:,5),'LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{176}Lu/^{177}Hf in upper mantle');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','Location','Southwest')
subplot(2,2,3);grid on;
plot(h_LM,R_average_SmNd_LM_f1(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_LM_f1(:,2),'LineWidth',2);
plot(h_LM,R_average_SmNd_LM_f1(:,3),'LineWidth',2);
plot(h_LM,R_average_SmNd_LM_f1(:,4),'LineWidth',2);
plot(h_LM,R_average_SmNd_LM_f1(:,5),'LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{147}Sm/^{144}Nd in crystalized lower mantle');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','Location','Northwest')
subplot(2,2,4);grid on;
plot(h_LM,R_average_SmNd_UM_f1(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_UM_f1(:,2),'LineWidth',2);
plot(h_LM,R_average_SmNd_UM_f1(:,3),'LineWidth',2);
plot(h_LM,R_average_SmNd_UM_f1(:,4),'LineWidth',2);
plot(h_LM,R_average_SmNd_UM_f1(:,5),'LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{147}Sm/^{144}Ndin upper mantle');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','Location','Southwest')

C_Nd_t0 = nan(n_crystalize+1,1);
C_Hf_t0 = nan(n_crystalize+1,1);
R_SmNd_t0 = nan(n_crystalize+1,1);
R_LuHf_t0 = nan(n_crystalize+1,1);
for i = 1:n_crystalize
    C_Nd_t0(i) = C_144Nd_BSE_t0;
    C_Hf_t0(i) = C_177Hf_BSE_t0;
    R_SmNd_t0(i)  = R_SmNd_CHUR_t0;
    R_LuHf_t0(i)  = R_LuHf_CHUR_t0;
end

figure(2);subplot(2,2,1);grid on;
plot(h_LM,R_average_LuHf_LM_f2(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_LM_f2(:,2),'LineWidth',2);
plot(h_LM,R_average_LuHf_LM_f2(:,3),'LineWidth',2);
plot(h_LM,R_average_LuHf_LM_f2(:,4),'LineWidth',2);
plot(h_LM,R_average_LuHf_LM_f2(:,5),'LineWidth',2);
plot(h_LM,R_LuHf_t0,'g-','LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{176}Lu/^{177}Hf in crystalized lower mantle');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','CHUR','Location','Northwest')
subplot(2,2,2);grid on;
plot(h_LM,R_average_LuHf_UM_f2(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_UM_f2(:,2),'LineWidth',2);
plot(h_LM,R_average_LuHf_UM_f2(:,3),'LineWidth',2);
plot(h_LM,R_average_LuHf_UM_f2(:,4),'LineWidth',2);
plot(h_LM,R_average_LuHf_UM_f2(:,5),'LineWidth',2);
plot(h_LM,R_LuHf_t0,'g-','LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{176}Lu/^{177}Hf in upper mantle');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','CHUR','Location','Southwest')
subplot(2,2,3);grid on;
plot(h_LM,R_average_SmNd_LM_f2(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_LM_f2(:,2),'LineWidth',2);
plot(h_LM,R_average_SmNd_LM_f2(:,3),'LineWidth',2);
plot(h_LM,R_average_SmNd_LM_f2(:,4),'LineWidth',2);
plot(h_LM,R_average_SmNd_LM_f2(:,5),'LineWidth',2);
plot(h_LM,R_SmNd_t0,'g-','LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{147}Sm/^{144}Nd in crystalized lower mantle');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','CHUR','Location','Northwest')
subplot(2,2,4);grid on;
plot(h_LM,R_average_SmNd_UM_f2(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_UM_f2(:,2),'LineWidth',2);
plot(h_LM,R_average_SmNd_UM_f2(:,3),'LineWidth',2);
plot(h_LM,R_average_SmNd_UM_f2(:,4),'LineWidth',2);
plot(h_LM,R_average_SmNd_UM_f2(:,5),'LineWidth',2);
plot(h_LM,R_SmNd_t0,'g-','LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{147}Sm/^{144}Nd in upper mantle');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','CHUR','Location','Southwest')

figure(3);subplot(2,2,1);grid on;
plot(h_LM,R_average_LuHf_LM_f1(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_LM_f2(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_LM_eq(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_LM_fra(:,1),'--','LineWidth',2); hold on;
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{176}Lu/^{177}Hf in crystalized lower mantle');
legend('Fractional crystallization, f = 60%','Fractional crystallization, f = 100%',...
    'Equilibirum crystallization','Completely fractional crystallization',...
    'Location','Northeast')
subplot(2,2,2);grid on;
plot(h_LM,R_average_LuHf_UM_f1(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_UM_f2(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_UM_eq(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_LuHf_UM_fra(:,1),'--','LineWidth',2); hold on;
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{176}Lu/^{177}Hf in crystalized upper mantle');
legend('Fractional crystallization, f = 60%','Fractional crystallization, f = 100%',...
    'Equilibirum crystallization','Completely fractional crystallization',...
    'Location','Southwest')

subplot(2,2,3);grid on;
plot(h_LM,R_average_SmNd_LM_f1(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_LM_f2(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_LM_eq(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_LM_fra(:,1),'--','LineWidth',2); hold on;
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{147}Sm/^{144}Nd in crystalized lower mantle');
legend('Fractional crystallization, f = 60%','Fractional crystallization, f = 100%',...
    'Equilibirum crystallization','Completely fractional crystallization',...
    'Location','Northeast')

subplot(2,2,4);grid on;
plot(h_LM,R_average_SmNd_UM_f1(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_UM_f2(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_UM_eq(:,1),'LineWidth',2); hold on;
plot(h_LM,R_average_SmNd_UM_fra(:,1),'--','LineWidth',2); hold on;
xlabel('Depth of crystalized lower mantle');
ylabel('Average ^{147}Sm/^{144}Nd in crystalized upper mantle');
legend('Fractional crystallization, f = 60%','Fractional crystallization, f = 100%',...
    'Equilibirum crystallization','Completely fractional crystallization',...
    'Location','Southwest')

%% Test the consistancy with Caro et al. (2005) Fig. 2

% Plot the evolution of Sm/Nd and Lu/Hf of residual melt during 
% crystallization of a magma ocean at lower-mantle pressures
x_LuHf_BSE = [0.1967,0.1967];
y_LuHf_BSE = [0,0.045];
x_SmNd_BSE = [0.1,0.24];
y_SmNd_BSE = [0.0332,0.0332];
data_Caro = xlsread('caro et al fig1.xlsx');
x_caro = data_Caro(:,1);y_caro = data_Caro(:,2);

figure(4);
plot(R_average_SmNd_UM_eq(:,1),R_average_LuHf_UM_eq(:,1),'k-','LineWidth',2); hold on;
plot(R_average_SmNd_UM_eq(:,2),R_average_LuHf_UM_eq(:,2),'c-','LineWidth',2);
plot(R_average_SmNd_UM_eq(:,3),R_average_LuHf_UM_eq(:,3),'m-','LineWidth',2);
plot(R_average_SmNd_UM_eq(:,4),R_average_LuHf_UM_eq(:,4),'g-','LineWidth',2);
plot(R_average_SmNd_UM_eq(:,5),R_average_LuHf_UM_eq(:,5),'r-','LineWidth',2);
plot(R_average_SmNd_UM_f2(:,1),R_average_LuHf_UM_f2(:,1),'k--','LineWidth',2); 
plot(R_average_SmNd_UM_f2(:,2),R_average_LuHf_UM_f2(:,2),'c--','LineWidth',2);
plot(R_average_SmNd_UM_f2(:,3),R_average_LuHf_UM_f2(:,3),'m--','LineWidth',2);
plot(R_average_SmNd_UM_f2(:,4),R_average_LuHf_UM_f2(:,4),'g--','LineWidth',2);
plot(R_average_SmNd_UM_f2(:,5),R_average_LuHf_UM_f2(:,5),'r--','LineWidth',2);
plot(x_LuHf_BSE,y_LuHf_BSE,'k--','LineWidth',2);
plot(x_SmNd_BSE,y_SmNd_BSE,'k--','LineWidth',2);
scatter(x_caro,y_caro,90,'filled','b');
xlim([0.1,0.24]); ylim([0.012,0.045]);
xlabel('^{147}Sm/^{144}Nd');
ylabel('^{176}Lu/^{177}Hf');
legend('5% CaPv in LM','6% CaPv in LM','7% CaPv in LM','8% CaPv in LM',...
    '9% CaPv in LM','Location','Northwest')

%% Plot concentration and ratio in liquid and solid during each crystallization
figure(6); 
C_Nd_t0 = nan(n_crystalize+1,1);
C_Hf_t0 = nan(n_crystalize+1,1);
R_SmNd_t0 = nan(n_crystalize+1,1);
R_LuHf_t0 = nan(n_crystalize+1,1);
for i = 1:n_crystalize
    C_Nd_t0(i) = C_144Nd_BSE_t0;
    C_Hf_t0(i) = C_177Hf_BSE_t0;
    R_SmNd_t0(i)  = R_SmNd_CHUR_t0;
    R_LuHf_t0(i)  = R_LuHf_CHUR_t0;
end
subplot(2,2,1)
plot(h_LM,C_Nd_liquid_f2(:,3),'r--','LineWidth',2);hold on;
plot(h_LM,C_Nd_solid_f2(:,3),'b-','LineWidth',2);
plot(h_LM,C_Nd_t0,'g-','LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('concentration of Nd');
legend('liquid','solid','C_{0}')

subplot(2,2,2)
plot(h_LM,C_Hf_liquid_f2(:,3),'r--','LineWidth',2);hold on;
plot(h_LM,C_Hf_solid_f2(:,3),'b-','LineWidth',2);
plot(h_LM,C_Hf_t0,'g-','LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('concentration of Hf');
legend('liquid','solid','C_{0}')

subplot(2,2,3)
plot(h_LM,R_liqid_SmNd_f2(:,3),'r--','LineWidth',2);hold on;
plot(h_LM,R_solid_SmNd_f2(:,3),'b-','LineWidth',2);
plot(h_LM,R_SmNd_t0,'g-','LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Ratio of Sm/Nd');
legend('liquid','solid','R_{t0}')

subplot(2,2,4)
plot(h_LM,R_liqid_LuHf_f2(:,3),'r--','LineWidth',2);hold on;
plot(h_LM,R_solid_LuHf_f2(:,3),'b-','LineWidth',2);
plot(h_LM,R_LuHf_t0,'g-','LineWidth',2);
xlabel('Depth of crystalized lower mantle');
ylabel('Ratio of Lu/Hf');
legend('liquid','solid','R_{t0}')
