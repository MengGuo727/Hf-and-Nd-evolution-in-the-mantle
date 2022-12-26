function[M_BSE_tp,M_BSE,M_MM,Mc_f,N_Hf177_BSE_tp,N_Hf_BSE_tp,N_Hf177_CC_tp,...
    N_Hf_CC_tp,N_Hf177_CC,N_Hf_CC,N_Hf177_MM,N_Hf177_BSE,N_Lu176_CC_tp,N_Lu176_BSE_t0,...
    N_Lu176_BSE_tp,N_Lu176_CC_t0,N_Lu176_BSE,N_Lu176_CC,N_Lu176_MM,N_Hf176_CC_tp,...
    N_Hf176_BSE_t0,N_Hf176_CC_t0,N_Hf176_BSE_tp,N_Hf176_BSE,N_Hf176_CC,...
    N_Hf176_MM,R_Hf_BSE,R_Hf_CC,R_Hf_MM] = ...
Prepare_isotopeInfo(t,nt,Mmp,Mcp,M_CC,C_Hf_BSE_tp,Fra_Hf177_tp,Na,M_Hf177,Fra_allHf,M_Hf,...
    C_Hf_CC_tp,C_Lu_CC_tp,M_Lu176,Fra_Lu176_tp,R_LuHf_CHUR_t0,lambda_Lu176,M_Hf176,...
    Fra_Hf176_tp,R_Hf_CHUR_t0)

% This function calculates the isotopic concentration/abundance in CC,MM,
% and BSE under the assumption of instant mixing in the mantle
% 
% Meng Guo, Yale University
% Summer 2020

%% Variables names and their meanings:
% t: time serie
% nt: length of time serie
% Na: Avogadro's constant
% lambda_Lu176: decay constant of 176Lu
% N_Isotope_Reservoir_Time: this term mean the # of atoms of an isotope in
   % a reservoir at certaint ime. Example: "N_Hf176_BSE_tp" means # of 176Hf
   % in BSE at present-day
% C_Isotope_Reservoir_Time: this term mean the concentration of an isotope in
   % a reservoir at certaint ime. Example: "C_Hf176_BSE_tp" means the 
   % concentration of 176Hf in BSE at present-day
% R_Isotope_Reservoir_Time: this term mean the ratio of two isotopes in
   % a reservoir. Example: "R_LuHf_BSE_t0" means the isotioic ratio of
   % 176Lu to 177Hf in BSE at the beginning of solar system
% R_Hf_CHUR_t0: initial 176Hf/177Hf, CHUR
% R_LuHf_CHUR_t0: initial 176Lu/177Hf, CHUR
% Fra_Isotope_Time: this terms mean the fraction of an isotope abundance to its 
   % element abundance. Example: "Fra_Hf177_tp" means the fraction of 177Hf
   % to all Hf at present-day
% M_Isotope_Reservoir_Time: this term mean the mass of an isotope/element in
   % a reservoir at certaint ime. Example: "M_Hf176_BSE_tp" means the mass of 
   % 176Hf in BSE at present-day
% Mmp: present-day mantle mass, in unit kg
% Mcp: present-day continental crust mass, in unit kg

%% Step 1. Calculate BSE and MM resevoirs masses
M_BSE_tp = Mmp + Mcp;% mass of BSE, in unit kg
M_BSE = nan(size(t)); 
for i = 1:nt
    M_BSE(i) = M_BSE_tp;
end

% calculate corresponding whole mantle mass
M_MM = nan(size(t)); 
for i = 1:nt
    M_MM(i) = M_BSE(i)/2 - M_CC(i); 
end
% depleted mantle mass = M_WholeMantlew - M_PM, where M_PM is half of M_BSE
M_DepletedM = nan(size(t)); 
for i = 1:nt
    M_DepletedM(i) = M_BSE(i)/2 - M_CC(i); 
end

%% Step 2. Calculate the stable isotope of daughter element 177Hf for CC, MM, BSE
Mc_f = M_CC/Mcp; % continental crust mass normalized to presetn-day mass
% calcualte # of177Hf atoms in the present-day BSE
[N_Hf177_BSE_tp] = C2N(C_Hf_BSE_tp,Fra_Hf177_tp,M_BSE_tp,Na,M_Hf177);% in unit # of atoms
% calcualte # of Hf atoms in the present-day BSE
[N_Hf_BSE_tp] = C2N(C_Hf_BSE_tp,Fra_allHf,M_BSE_tp,Na,M_Hf);% in unit # of atoms
% calcualte # of177Hf atoms in the present-day CC
[N_Hf177_CC_tp] = C2N(C_Hf_CC_tp,Fra_Hf177_tp,Mcp,Na,M_Hf177);% in unit # of atoms
% calcualte # of Hf atoms in the present-day CC
[N_Hf_CC_tp] = C2N(C_Hf_CC_tp,Fra_allHf,Mcp,Na,M_Hf);% in unit # of atoms

N_Hf177_CC = nan(size(t));
N_Hf_CC = nan(size(t));
N_Hf177_MM = nan(size(t));
N_Hf177_BSE = nan(size(t));
for i = 1:nt
    N_Hf177_BSE(i) = N_Hf177_BSE_tp;
    N_Hf177_CC(i) = N_Hf177_CC_tp * Mc_f(i);
    N_Hf_CC(i) = N_Hf_CC_tp * Mc_f(i);
    N_Hf177_MM(i) = N_Hf177_BSE(i)- N_Hf177_CC(i);
end

%% Step 3. Calculate the parent isotope 176Lu for CC, MM, BSE
[N_Lu176_CC_tp] = C2N(C_Lu_CC_tp,Fra_Lu176_tp,Mcp,Na,M_Lu176);% in unit # of atoms
N_Lu176_BSE_t0 = R_LuHf_CHUR_t0 * N_Hf177_BSE(1); 
N_Lu176_BSE_tp = N_Lu176_BSE_t0 * exp(-lambda_Lu176 * t(nt));% in unit # of atoms
N_Lu176_CC_t0 =  N_Lu176_CC_tp  * exp(lambda_Lu176 * t(nt));% in unit # of atoms

N_Lu176_BSE = nan(size(t));
N_Lu176_BSE(1)= N_Lu176_BSE_t0 ;
N_Lu176_CC = nan(size(t));
N_Lu176_MM = nan(size(t));
for i = 2:nt
    N_Lu176_BSE(i) = N_Lu176_BSE_t0 * exp(-lambda_Lu176*t(i));
end
for i = 1:nt
    N_Lu176_CC(i) = N_Lu176_CC_t0 * exp(-lambda_Lu176*t(i))* Mc_f(i);
    N_Lu176_MM(i) = N_Lu176_BSE(i)- N_Lu176_CC(i);
end


%% Step 4. Calulate the daughter isotope 176Hf in CC, MM, BSE
[N_Hf176_CC_tp] = C2N(C_Hf_CC_tp,Fra_Hf176_tp,Mcp,Na,M_Hf176);% 5.2558% Hf are 176Hf, 3.7ppm Hf in Tp CC, in unit # of atoms
N_Hf176_BSE_t0 = R_Hf_CHUR_t0* N_Hf177_BSE(1);
N_Hf176_CC_t0 =  N_Hf176_CC_tp - N_Lu176_CC_tp*(exp(lambda_Lu176 * t(nt))-1);
N_Hf176_BSE_tp =  N_Hf176_BSE_t0 + N_Lu176_BSE_tp*(exp(lambda_Lu176 * t(nt))-1);

N_Hf176_BSE = nan(size(t));
N_Hf176_BSE(1) = N_Hf176_BSE_t0;
N_Hf176_CC = nan(size(t));
N_Hf176_CC(1) = N_Hf176_CC_t0;
N_Hf176_MM = nan(size(t));
for i = 2:nt
    N_Hf176_BSE(i) = N_Hf176_BSE_t0 + N_Lu176_BSE_t0*(1-exp(-lambda_Lu176 * t(i)));
end
for i = 1:nt
    N_Hf176_CC(i) = (N_Hf176_CC_t0 + N_Lu176_CC_t0*(1-exp(-lambda_Lu176 * t(i))))* Mc_f(i);
    N_Hf176_MM(i) = N_Hf176_BSE(i)- N_Hf176_CC(i);
end

%% Step 4. Calulate 176Hf/177Hf ratio
R_Hf_CC = nan(size(t));
R_Hf_MM= nan(size(t));
R_Hf_BSE = nan(size(t));
for i = 1:nt
    if (M_CC(i) == 0)
        R_Hf_CC(i) =0;
    else
        R_Hf_CC(i) = N_Hf176_CC(i) / N_Hf177_CC(i);
    end
    R_Hf_MM(i) = N_Hf176_MM(i) / N_Hf177_MM(i);
    R_Hf_BSE(i) = N_Hf176_BSE(i) / N_Hf177_BSE(i);
end