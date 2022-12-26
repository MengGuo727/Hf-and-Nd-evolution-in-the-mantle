function[M_BSE_tp,M_BSE,M_MM,Mc_f,N_Stable_CC_tp,N_DaughterElement_CC_tp,N_Stable_CC,...
    N_DaughterElement_CC,N_Stable_MM,N_Parent_CC_tp,...
    N_Parent_CC_t0,N_Parent_CC,N_Parent_MM,N_Daughter_CC_tp,...
    N_Daughter_CC_t0,N_Daughter_CC,N_Daughter_MM,R_BSE,R_CC,R_MM] = ...
Prepare_isotopeInfo_T0toTp(t,nt,Mmp,Mcp,M_CC,Fra_Stable_tp,Na,M_Stable,Fra_all,M_DaughterElement,...
    C_DaughterElement_CC_tp,C_ParentElement_CC_tp,M_Parent,Fra_Parent_tp,lambda_Parent,M_Daughter,...
    Fra_Daughter_tp,N_Stable_BSE,N_Parent_BSE,N_Daughter_BSE)

% This function calculates the isotopic concentration/abundance in CC,MM,
% and BSE under the assumption of instant mixing in the mantle
% 
% Meng Guo, Yale University
% Summer 2020

%% Variables names and their meanings:
% t: time serie
% nt: length of time serie
% Na: Avogadro's constant
% lambda_Parent: decay constant of 176ParentElement
% N_Isotope_Reservoir_Time: this term mean the # of atoms of an isotope in
   % a reservoir at certaint ime. Example: "N_Daughter_CC_tp" means # of 176DaughterElement
   % in continental crust at present-day
% C_Isotope_Reservoir_Time: this term mean the concentration of an isotope in
   % a reservoir at certaint ime. Example: "C_Daughter_BSE_tp" means the 
   % concentration of 176DaughterElement in BSE at present-day
% R_Isotope_Reservoir_Time: this term mean the ratio of two isotopes in
   % a reservoir. Example: "R_ParentElementDaughterElement_BSE_t0" means the isotioic ratio of
   % 176ParentElement to 177DaughterElement in BSE at the beginning of solar system
% Fra_Isotope_Time: this terms mean the fraction of an isotope abundance to its 
   % element abundance. Example: "Fra_Stable_tp" means the fraction of 177DaughterElement
   % to all DaughterElement at present-day
% M_Isotope_Reservoir_Time: this term mean the mass of an isotope/element in
   % a reservoir at certaint ime. Example: "M_Daughter_BSE_tp" means the mass of 
   % 176DaughterElement in BSE at present-day
% Mmp: present-day mantle mass, in unit kg
% Mcp: present-day continental crust mass, in unit kg

% Stable: stalbe isotope; Parent: parent isotope; Daughter: daughter isotop
% ParentElement: all the isotope of parent element; DaughterElement: all the isotope of daughter element; 
% EX.: for Lu-Hf system, Stable:177Hf; Parent: 176Lu; Daughter: 176Hf; ParentElement: Lu; DaughterElement:Hf

%% Step 1. Calculate BSE and MM resevoirs masses
M_BSE_tp = Mmp + Mcp;% mass of BSE, in unit kg
M_BSE = nan(size(t)); 
for i = 1:nt
    M_BSE(i) = M_BSE_tp;
end

% calculate corresponding whole mantle mass
M_MM = nan(size(t)); 
for i = 1:nt
    M_MM(i) = M_BSE(i) - M_CC(i); 
end

%% Step 2. Calculate the stable isotope of daughter element 177DaughterElement for CC, MM, BSE
Mc_f = M_CC/Mcp; % continental crust mass normalized to presetn-day mass
% calcualte # of177DaughterElement atoms in the present-day CC
[N_Stable_CC_tp] = C2N(C_DaughterElement_CC_tp,Fra_Stable_tp,Mcp,Na,M_Stable);% in unit # of atoms
% calcualte # of DaughterElement atoms in the present-day CC
[N_DaughterElement_CC_tp] = C2N(C_DaughterElement_CC_tp,Fra_all,Mcp,Na,M_DaughterElement);% in unit # of atoms

N_Stable_CC = nan(size(t));
N_DaughterElement_CC = nan(size(t));
N_Stable_MM = nan(size(t));
for i = 1:nt
    N_Stable_CC(i) = N_Stable_CC_tp * Mc_f(i);
    N_DaughterElement_CC(i) = N_DaughterElement_CC_tp * Mc_f(i);
    N_Stable_MM(i) = N_Stable_BSE(i)- N_Stable_CC(i);
end

%% Step 3. Calculate the parent isotope 176ParentElement for CC, MM, BSE
[N_Parent_CC_tp] = C2N(C_ParentElement_CC_tp,Fra_Parent_tp,Mcp,Na,M_Parent);% in unit # of atoms
N_Parent_CC_t0 =  N_Parent_CC_tp  * exp(lambda_Parent * t(nt));% in unit # of atoms
N_Parent_CC = nan(size(t));
N_Parent_MM = nan(size(t));
for i = 1:nt
    N_Parent_CC(i) = N_Parent_CC_t0 * exp(-lambda_Parent*t(i))* Mc_f(i);
    N_Parent_MM(i) = N_Parent_BSE(i)- N_Parent_CC(i);
end


%% Step 4. Calulate the daughter isotope 176DaughterElement in CC, MM, BSE
[N_Daughter_CC_tp] = C2N(C_DaughterElement_CC_tp,Fra_Daughter_tp,Mcp,Na,M_Daughter);% 5.2558% DaughterElement are 176DaughterElement, 3.7ppm DaughterElement in Tp CC, in unit # of atoms
N_Daughter_CC_t0 =  N_Daughter_CC_tp - N_Parent_CC_tp*(exp(lambda_Parent * t(nt))-1);
N_Daughter_CC = nan(size(t));
N_Daughter_CC(1) = N_Daughter_CC_t0;
N_Daughter_MM = nan(size(t));
for i = 1:nt
    N_Daughter_CC(i) = (N_Daughter_CC_t0 + N_Parent_CC_t0*(1-exp(-lambda_Parent * t(i))))* Mc_f(i);
    N_Daughter_MM(i) = N_Daughter_BSE(i)- N_Daughter_CC(i);
end

%% Step 4. Calulate 176DaughterElement/177DaughterElement ratio
R_CC = nan(size(t));
R_MM= nan(size(t));
R_BSE = nan(size(t));
for i = 1:nt
    if (M_CC(i) == 0)
        R_CC(i) =0;
    else
        R_CC(i) = N_Daughter_CC(i) / N_Stable_CC(i);
    end
    R_MM(i) = N_Daughter_MM(i) / N_Stable_MM(i);
    R_BSE(i) = N_Daughter_BSE(i) / N_Stable_BSE(i);
end