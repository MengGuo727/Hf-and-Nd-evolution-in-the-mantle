function[M_BSE_tp,M_BSE,M_MM,Mc_f,N_Hf177_BSE_tp,N_Hf_BSE_tp,N_Hf177_CC_tp,...
    N_Hf_CC_tp,N_Hf177_CC,N_Hf_CC,N_Hf177_MM,N_Hf177_BSE,N_Lu176_CC_tp,N_Lu176_BSE_t0,...
    N_Lu176_BSE_tp,N_Lu176_CC_t0,N_Lu176_BSE,N_Lu176_CC,N_Lu176_MM,N_Hf176_CC_tp,...
    N_Hf176_BSE_t0,N_Hf176_CC_t0,N_Hf176_BSE_tp,N_Hf176_BSE,N_Hf176_CC,N_Hf176_MM,...
    R_Hf_CC,R_Hf_MM,R_Hf_BSE,e176Hf_CC,e176Hf_MM,e176Hf_BSE] = ...
Instant_mixing(t,nt,Mmp,Mcp,M_CC,C_Hf_BSE_tp,Fra_Hf177_tp,Na,M_Hf177,Fra_allHf,M_Hf,...
    C_Hf_CC_tp,C_Lu_CC_tp,M_Lu176,Fra_Lu176_tp,R_LuHf_CHUR_t0,lambda_Lu176,M_Hf176,...
    Fra_Hf176_tp,R_Hf_CHUR_t0)

% This function calculates instant mixing in the mantle
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

%% Step 1. Calculate isotopic abundance in CC, MM, and BSE assume instant mixing
[M_BSE_tp,M_BSE,M_MM,Mc_f,N_Hf177_BSE_tp,N_Hf_BSE_tp,N_Hf177_CC_tp,...
    N_Hf_CC_tp,N_Hf177_CC,N_Hf_CC,N_Hf177_MM,N_Hf177_BSE,N_Lu176_CC_tp,N_Lu176_BSE_t0,...
    N_Lu176_BSE_tp,N_Lu176_CC_t0,N_Lu176_BSE,N_Lu176_CC,N_Lu176_MM,N_Hf176_CC_tp,...
    N_Hf176_BSE_t0,N_Hf176_CC_t0,N_Hf176_BSE_tp,N_Hf176_BSE,N_Hf176_CC,...
    N_Hf176_MM,R_Hf_BSE,R_Hf_CC,R_Hf_MM] = ...
Prepare_isotopeInfo(t,nt,Mmp,Mcp,M_CC,C_Hf_BSE_tp,Fra_Hf177_tp,Na,M_Hf177,Fra_allHf,M_Hf,...
    C_Hf_CC_tp,C_Lu_CC_tp,M_Lu176,Fra_Lu176_tp,R_LuHf_CHUR_t0,lambda_Lu176,M_Hf176,...
    Fra_Hf176_tp,R_Hf_CHUR_t0);

%% Step 2. Calcualte the ratio of Lu/Hf and eHf through time
% Calculate the e176Hf assume no mantle mixing
e176Hf_CC = Epsilon(R_Hf_CC, R_Hf_BSE,1);
e176Hf_MM = Epsilon(R_Hf_MM, R_Hf_BSE,1);
e176Hf_BSE = Epsilon(R_Hf_BSE, R_Hf_BSE,1);