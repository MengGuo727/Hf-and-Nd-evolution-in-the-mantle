% This function calculates crust mantle differentiation with finite mixing
% in the mantle.
function[C_BSE,N_dCC,N_dCC_total,N_UM,N_LM,N_dCC_ini,N_DM,N_CR,N_DM_total,...
    N_OM,N_CR_total,N_BSE,N_WM,N_PC,N_RCC,N_OM_consumed,...
    N_RCC_total,N_PC_ini,N_RCC_ini,N_DM_ini,N_CR_ini,M_Reservoir,M_LM,M_UM,...
    l_LM,l_UM,M_subReservoir_total,N_OM_consumed_total,F1,R_WM,R_OM,...
    R_reservoir,e_WM,l_DM_number,l_CR_number,l_RCC_number,...
    l_DM_flag_matrix,l_CR_flag_matrix,l_RCC_flag_matrix,F_consumed,e_sample] = ...
Crust_Extraction_with_Finite_Mantle_Mixing_Cloud(N_Stable_BSE,N_Parent_BSE,N_Daughter_BSE,...
    M_BSE,t,nt,dt,Na,M_Stable,M_Daughter,M_Parent,Mud,...
    Z,F1_tp,F2_0,rhom_MM,alpha,alpha_OM,alpha_DM,lambda_Parent,M_CC,l_critical,Mdd,...
    M_LM_t0,M_UM_t0,N_Parent_LM_t0,N_Parent_UM_t0,N_Stable_LM_t0,N_Stable_UM_t0,...
    N_Daughter_LM_t0,N_Daughter_UM_t0,t_overturn,D_Parent1,D_Daughter1,R_BSE,n_sample,l_samplingbox,ts,n_test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables names and their meanings:
% t: time serie
% nt: length of time serie
% dt: length of timestep
% Na: Avogadro's constant
% F1_tp: present-day mantle partial melting fraction
% F2: fraction of second partial melting to generate dCC and residue-crust
% Z: mantle melting depth, in unit km
% rho_DM: mantle density, in unit kg/m3
% alpha: DMand CR blob flatten constant, can be changed to fit observation
% lambda_Parent: decay constant of 176Lu

% N_Isotope_Reservoir_Time: this term means the # of atoms of an isotope in
% a reservoir at certaint ime. Example: "N_Daughter_BSE_tp" means # of 176Hf
% in BSE at present-day

% C_Isotope_Reservoir_Time: this term means the concentration of an isotope
% in a reservoir at certaint ime. Example: "C_Daughter_BSE_tp" means the
% concentration of 176Hf in BSE at present-day

% R_Isotope_Reservoir_Time: this term means the ratio of two isotopes in
% a reservoir. Example: "R_LuHf_BSE_t0" means the isotioic ratio of
% 176Lu to 177Hf in BSE at the beginning of solar system

% M_Isotope_Reservoir_Time: this term means the mass of an isotope/element
% in a reservoir at certaint ime. Example: "M_Daughter_BSE_tp" means the mass of
% 176Hf in BSE at present-day

% dCC: the generation rate of continental crust [kg/Gyr]

% l_critical: the critical length of blob; For a blob has l<l_critical, it can remelt

% All the concentration or # of atoms matrixes are constructed as:
  % first column: Stable, second column: Daughter, and third column: Parent
  % for example: C_BSE = [C_Stable_BSE,C_Daughter_BSE,C_Parent_BSE]

% The mass or length matrixes for different reservoirs are constructed as:
  % M_Reservoir = [M_BSE,M_WM,M_WM_consumed,M_OM,M_DM,M_PC,M_CR,M_dCC,M_RCC]
  % "BSE" means the bulk silicate Earth
  % "WM" means the whole mantle composition that can remelt
  % "OM" means the primitive mantle composition that can remelt
  % "DM" means the depleted mantle (the crystalized mantle after first melting stage)
  % "PC" means the proto crust (the liquid melt after first melting stage)
  % "CR" means the crustal residue (the crystalized mantle after second melting stage)
  % "dCC" means the newly generated continental crust  (the liquid melt after second melting stage)
  % "RCC" means the recycled continental crust
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1. Calculate melting fraction for frist stage (F1) using melting depth (Z)
Z_tp = Z(nt); % 10% melting at present-day
F1 = Z.* F1_tp/Z_tp;% knowing today's F1(tp) = 10% and Z(tp) = Z(nt)
F2 = nan(size(t));
for k = 1:nt
    if (Mud(k) == 0)
        F1(k) = 0;
        F2(k) = 0;
    else
        F2(k) = F2_0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2. Calculate the elements concentrations in the BSE throughout Earth history
[C_Stable_BSE] = N2C(nt,N_Stable_BSE,M_Stable,Na,M_BSE,t,1);
[C_Daughter_BSE] = N2C(nt,N_Daughter_BSE,M_Daughter,Na,M_BSE,t,1);
[C_Parent_BSE] = N2C(nt,N_Parent_BSE,M_Parent,Na,M_BSE,t,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3. Assume initial melting is followed by F2= 5% secondary melting to
% generate dCC from OM,calculate the reservoir mass and isotipic concentration
% for PC, DM, CC and then CR.

% Step 3.1 Calculate the mass of each reservoir
M_dCC = Mud * dt; % mass of newly generated continenntal crust at each timestep [kg] 
M_PC = M_dCC ./ F2; % mass of proto-crust generated at each timestep [kg]
M_CR = M_PC - M_dCC; % mass of residual crust generated at each timestep [kg]
M_WM_consumed = M_PC ./ F1; % mass of consumed primitive mantle at each timestep [kg]
M_DM = M_WM_consumed - M_PC; % mass of depleted mantle blob generated at each timestep [kg]
M_RCC = Mdd * dt; % mass of recycled continental crust for each timestep
for k = 1:nt
    if (M_dCC(k) == 0)
        M_PC(k) = 0;
        M_WM_consumed(k) = 0;
        M_DM(k) = 0;
        M_CR(k) = 0;
    end
end

% Step 3.2 Create empty matrix to store the isotopic concentration in each reservoir
N_Stable_PC = zeros(size(t)); N_Parent_PC = zeros(size(t)); N_Daughter_PC = zeros(size(t));
N_Stable_DM = zeros(size(t)); N_Parent_DM = zeros(size(t)); N_Daughter_DM = zeros(size(t));
C_Stable_WM = zeros(size(t)); C_Parent_WM = zeros(size(t)); C_Daughter_WM = zeros(size(t));
N_Stable_WM = zeros(size(t)); N_Parent_WM = zeros(size(t)); N_Daughter_WM = zeros(size(t)); 
N_Stable_CR = zeros(size(t)); N_Parent_CR = zeros(size(t)); N_Daughter_CR = zeros(size(t));
N_Daughter_RCC = zeros(size(t)); N_Parent_RCC = zeros(size(t)); N_Stable_RCC = zeros(size(t));
N_Stable_OM = zeros(size(t)); N_Parent_OM = zeros(size(t)); N_Daughter_OM = zeros(size(t));

l_DM = zeros(size(t)); l_CR = zeros(size(t)); l_RCC = zeros(size(t));
l_DM_flag = zeros(size(t)); l_CR_flag = zeros(size(t)); l_RCC_flag = zeros(size(t));
l_RCC_flag_matrix = zeros(4568,4568); l_DM_flag_matrix = zeros(4568,4568); l_CR_flag_matrix = zeros(4568,4568);
% e_DM = nan(4568,4568); e_CR = nan(4568,4568); e_RCC = nan(4568,4568); 
l_DM_number = zeros(4568,4568); l_CR_number = zeros(4568,4568); l_RCC_number = zeros(4568,4568);
M_WM = zeros(size(t)); M_OM = zeros(size(t)); M_DM_total = zeros(size(t));
M_LM = zeros(size(t)); M_UM = zeros(size(t));l_LM = zeros(size(t));l_UM = zeros(size(t));
N_Stable_LM = zeros(size(t)); N_Parent_LM = zeros(size(t)); N_Daughter_LM = zeros(size(t));
N_Stable_UM = zeros(size(t)); N_Parent_UM = zeros(size(t)); N_Daughter_UM = zeros(size(t));
M_CR_total = zeros(size(t)); M_RCC_total = zeros(size(t));
N_Stable_WM_consumed = zeros(size(t)); N_Daughter_WM_consumed = zeros(size(t)); N_Parent_WM_consumed = zeros(size(t));
N_Stable_OM_consumed = zeros(size(t)); N_Daughter_OM_consumed = zeros(size(t)); N_Parent_OM_consumed = zeros(size(t));
N_Stable_OM_consumed_total = zeros(size(t)); N_Daughter_OM_consumed_total = zeros(size(t)); N_Parent_OM_consumed_total = zeros(size(t)); 
N_Stable_dCC = zeros(size(t)); N_Daughter_dCC = zeros(size(t)); N_Parent_dCC = zeros(size(t));
M_WM_consumed_total = zeros(size(t));
N_Stable_dCC_consumed = zeros(size(t)); N_Parent_dCC_consumed = zeros(size(t)); N_Daughter_dCC_consumed= zeros(size(t));
N_Stable_DM_consumed = zeros(size(t)); N_Parent_DM_consumed = zeros(size(t)); N_Daughter_DM_consumed= zeros(size(t));
N_Stable_CR_consumed = zeros(size(t)); N_Parent_CR_consumed = zeros(size(t)); N_Daughter_CR_consumed= zeros(size(t));
N_Stable_RCC_consumed = zeros(size(t)); N_Parent_RCC_consumed = zeros(size(t)); N_Daughter_RCC_consumed= zeros(size(t));

N_Daughter_PC_ini=zeros(size(t)); N_Stable_PC_ini=zeros(size(t)); N_Parent_PC_ini=zeros(size(t));
N_Daughter_RCC_ini=zeros(size(t)); N_Stable_RCC_ini=zeros(size(t)); N_Parent_RCC_ini=zeros(size(t));
N_Daughter_dCC_ini=zeros(size(t)); N_Stable_dCC_ini=zeros(size(t)); N_Parent_dCC_ini=zeros(size(t));
N_Daughter_DM_ini=zeros(size(t)); N_Stable_DM_ini=zeros(size(t)); N_Parent_DM_ini=zeros(size(t));
N_Daughter_CR_ini=zeros(size(t)); N_Stable_CR_ini=zeros(size(t)); N_Parent_CR_ini=zeros(size(t));

N_Stable_PC_matrix = zeros(nt,nt); N_Parent_PC_matrix = zeros(nt,nt); N_Daughter_PC_matrix = zeros(nt,nt);
N_Stable_DM_matrix = zeros(nt,nt); N_Parent_DM_matrix = zeros(nt,nt); N_Daughter_DM_matrix = zeros(nt,nt);
N_Stable_CR_matrix = zeros(nt,nt); N_Parent_CR_matrix = zeros(nt,nt); N_Daughter_CR_matrix = zeros(nt,nt);
N_Daughter_RCC_matrix = zeros(nt,nt); N_Parent_RCC_matrix = zeros(nt,nt); N_Stable_RCC_matrix = zeros(nt,nt);
N_Daughter_dCC_matrix = zeros(nt,nt); N_Parent_dCC_matrix = zeros(nt,nt); N_Stable_dCC_matrix = zeros(nt,nt);

C_Stable_PC_matrix = zeros(nt,nt); C_Parent_PC_matrix = zeros(nt,nt); C_Daughter_PC_matrix = zeros(nt,nt);
C_Stable_DM_matrix = zeros(nt,nt); C_Parent_DM_matrix = zeros(nt,nt); C_Daughter_DM_matrix = zeros(nt,nt);
C_Stable_CR_matrix = zeros(nt,nt); C_Parent_CR_matrix = zeros(nt,nt); C_Daughter_CR_matrix = zeros(nt,nt);
C_Daughter_RCC_matrix = zeros(nt,nt); C_Parent_RCC_matrix = zeros(nt,nt); C_Stable_RCC_matrix = zeros(nt,nt);
C_Daughter_dCC_matrix = zeros(nt,nt); C_Parent_dCC_matrix = zeros(nt,nt); C_Stable_dCC_matrix = zeros(nt,nt);
C_Stable_LM = zeros(size(t)); C_Parent_LM = zeros(size(t)); C_Daughter_LM = zeros(size(t));
C_Stable_UM = zeros(size(t)); C_Parent_UM = zeros(size(t)); C_Daughter_UM = zeros(size(t));

% Step 3.3 At first timestep, we consider mantle just experiences magma ocean
% solidification. Both upper and lower mantle can participate in melting, 
% which are considered to be the whole remeltable mantle composition (WM)
M_LM(1) = M_LM_t0; M_UM(1) = M_UM_t0;
l_LM(1) = 2240000;
l_UM(1) = 660000;
N_Stable_LM(1) = N_Stable_LM_t0; N_Parent_LM(1) = N_Parent_LM_t0; N_Daughter_LM(1) = N_Daughter_LM_t0;
N_Stable_UM(1) = N_Stable_UM_t0; N_Parent_UM(1) = N_Parent_UM_t0; N_Daughter_UM(1) = N_Daughter_UM_t0;
N_Stable_WM(1) =  N_Stable_UM(1) + N_Stable_LM(1);
N_Parent_WM(1) =  N_Parent_UM(1) + N_Parent_LM(1);
N_Daughter_WM(1) =  N_Daughter_UM(1) + N_Daughter_LM(1);
M_WM(1) = M_UM(1) + M_LM(1); % both lower and upper mantle can participate in remelting
C_Stable_WM(1) = N2C(nt,N_Stable_WM(1),M_Stable,Na,M_WM(1),t,2);
C_Parent_WM(1) = N2C(nt,N_Parent_WM(1),M_Parent,Na,M_WM(1),t,2);
C_Daughter_WM(1) = N2C(nt,N_Daughter_WM(1),M_Daughter,Na,M_WM(1),t,2);
R_WM = nan(size(t));
R_WM(1) = N_Daughter_WM(1) / N_Stable_WM(1);
e_WM = nan(size(t));
e_WM(1) = Epsilon(R_WM(1), R_BSE(1),1);
% calcualte mass of each subreservoir after magma ocean solidification
M_OM(1) = M_UM(1) + M_LM(1);
N_Stable_OM(1)= N_Stable_UM(1) +N_Stable_LM(1);
N_Parent_OM(1)= N_Parent_UM(1) + N_Parent_LM(1);
N_Daughter_OM(1)= N_Daughter_UM(1) + N_Daughter_LM(1);
R_OM = nan(size(t));
R_OM(1) = N_Daughter_OM(1) / N_Stable_OM(1);

N_Stable_DM_total = zeros(size(t));N_Stable_DM_total(1)= N_Stable_DM(1);
N_Stable_CR_total = zeros(size(t));N_Stable_CR_total(1)= N_Stable_CR(1);
N_Stable_dCC_total = zeros(size(t));N_Stable_dCC_total(1)= N_Stable_dCC(1);
N_Stable_RCC_total = zeros(size(t));N_Stable_RCC_total(1)= N_Stable_RCC(1);

N_Parent_DM_total = zeros(size(t));N_Parent_DM_total(1)= N_Parent_DM(1);
N_Parent_CR_total = zeros(size(t));N_Parent_CR_total(1)= N_Parent_CR(1);
N_Parent_dCC_total = zeros(size(t));N_Parent_dCC_total(1)= N_Parent_dCC(1);
N_Parent_RCC_total = zeros(size(t));N_Parent_RCC_total(1)= N_Parent_RCC(1);

N_Daughter_DM_total = zeros(size(t));N_Daughter_DM_total(1)= N_Daughter_DM(1);
N_Daughter_CR_total = zeros(size(t));N_Daughter_CR_total(1)= N_Daughter_CR(1);
N_Daughter_dCC_total = zeros(size(t));N_Daughter_dCC_total(1)= N_Daughter_dCC(1);
N_Daughter_RCC_total = zeros(size(t));N_Daughter_RCC_total(1)= N_Daughter_RCC(1);

F_consumed = zeros(size(t));

fraction = 1;

D_Daughter2 = D_Daughter1; D_Parent2 = D_Parent1; % Assume D2 = D1

M_RCC(round(ts/dt+1)) = 0; 

e_sample = nan(n_test,n_sample);  
M_sample = nan(n_test,n_sample);  
N_Stable_sample = nan(n_test,n_sample);  
N_Daughter_sample = nan(n_test,n_sample);  
k = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4. Track the number of elements in each subreservoir during crustal growth
for i = 2:nt
    
    if (i == 644)
        i = i;
    end
    
    % Step 4.1 Calculate the fraction of total remeltable mantle consumed at this time step
    F_consumed(i) = M_WM_consumed(i) / M_WM(i-1);
    % in case of very rapid crustal generation, fix F1 to avoid over consumption of mantle
    if (F_consumed(i) > 1 )
        F_consumed(i) = 1;
        M_WM_consumed(i) = M_WM(i-1);
        F1(i) = M_PC(i) / M_WM_consumed(i);
    end
    
    % Step 4.2 Calculate the proto-crust (PC) isotope concentrations using F1,
    %  D1, and the isotope concentration of whole remeltable mantle (WM). 
    % Equation: C(l) = C(0)/(F+D-DF)
    if (F1(i) > 0)
        C_Stable_PC = C_Stable_WM(i-1)/(F1(i) + D_Daughter1 - D_Daughter1*F1(i));
        C_Parent_PC = C_Parent_WM(i-1)/(F1(i) + D_Parent1 - D_Parent1*F1(i));
        C_Daughter_PC = C_Daughter_WM(i-1)/(F1(i) + D_Daughter1 - D_Daughter1*F1(i));
    else % in case of no crustal generation at this time step 
        C_Stable_PC=0;
        C_Parent_PC=0;
        C_Daughter_PC=0;
    end
    % Step 4.3 change concentration to number of atoms
    [N_Stable_PC(i)] = C2N(C_Stable_PC,fraction,M_PC(i),Na,M_Stable);
    [N_Parent_PC(i)] = C2N(C_Parent_PC,fraction,M_PC(i),Na,M_Parent);
    [N_Daughter_PC(i)] = C2N(C_Daughter_PC,fraction,M_PC(i),Na,M_Daughter);
    N_Daughter_PC_ini(i) = N_Daughter_PC(i);
    N_Stable_PC_ini(i) = N_Stable_PC(i);
    N_Parent_PC_ini(i) = N_Parent_PC(i);
    
    % Check: calculate the decay of isotopes in the sub proto-crust reservoir from this time step before
    for j = 1:(i-1)
        N_Daughter_PC(j) = N_Daughter_PC(j) + N_Parent_PC(j)*(1-exp(-lambda_Parent*dt));
        N_Parent_PC(j) = N_Parent_PC(j)*exp(-lambda_Parent*dt);
    end   
    % Save in matrix to tack throguh time
    N_Stable_PC_matrix(1:i,i) = N_Stable_PC(1:i); 
    N_Parent_PC_matrix(1:i,i) = N_Parent_PC(1:i); 
    N_Daughter_PC_matrix(1:i,i) = N_Daughter_PC(1:i);
    C_Stable_PC_matrix(1:i,i) = N2C(i,N_Stable_PC(1:i),M_Stable,Na,M_PC(1:i),t(1:i),1);
    C_Parent_PC_matrix(1:i,i) = N2C(i,N_Parent_PC(1:i),M_Parent,Na,M_PC(1:i),t(1:i),1);
    C_Daughter_PC_matrix(1:i,i) = N2C(i,N_Daughter_PC(1:i),M_Daughter,Na,M_PC(1:i),t(1:i),1);

    % Step 4.3 Calculate the # of atoms in the remaining original primative mantle (OM = LM + UM) 
    % Step 4.3.1 calculate the isotopes consumed during mantle melting
    if (M_LM(i-1) > 0)
        N_Stable_LM_consumed = N_Stable_LM(i-1)*F_consumed(i);
        N_Daughter_LM_consumed = N_Daughter_LM(i-1)*F_consumed(i);
        N_Parent_LM_consumed = N_Parent_LM(i-1)*F_consumed(i);
        M_LM(i) = M_LM(i-1)*(1-F_consumed(i));
        l_LM(i) = l_LM(i-1)*(1-F_consumed(i))/alpha_OM(i);
        % Step 4.3.2  calculate the decay of primitive mantle
        N_Daughter_LM(i) = N_Daughter_LM(i-1) + (N_Parent_LM(i-1)*(1-exp(-lambda_Parent*dt)));
        N_Parent_LM(i) = N_Parent_LM(i-1) * exp(-lambda_Parent*dt);
        % Step 4.3.3 calculate the leftover isotopes in primitive mantle after the generation of CC
        N_Stable_LM(i) = N_Stable_LM(i-1) - N_Stable_LM_consumed;
        N_Daughter_LM(i) = N_Daughter_LM(i) - N_Daughter_LM_consumed;
        N_Parent_LM(i) = N_Parent_LM(i) - N_Parent_LM_consumed;
    else
        M_LM(i) = 0;
        N_Stable_LM(i) = 0;
        N_Daughter_LM(i) = 0;
        N_Parent_LM(i) = 0;
        l_LM(i) = 0;
    end
    if (M_LM(i) == 0)
        N_Stable_LM(i) = 0;
        N_Daughter_LM(i) = 0;
        N_Parent_LM(i) = 0;
        l_LM(i) = 0;
    end
    
    if (M_UM(i-1) > 0)
        N_Stable_UM_consumed = N_Stable_UM(i-1)*F_consumed(i);
        N_Daughter_UM_consumed = N_Daughter_UM(i-1)*F_consumed(i);
        N_Parent_UM_consumed = N_Parent_UM(i-1)*F_consumed(i);
        M_UM(i) = M_UM(i-1)*(1-F_consumed(i));
        l_UM(i) = l_UM(i-1)*(1-F_consumed(i))/alpha_OM(i);
        % Step 4.3.2  calculate the decay of primitive mantle
        N_Daughter_UM(i) = N_Daughter_UM(i-1) + (N_Parent_UM(i-1)*(1-exp(-lambda_Parent*dt)));
        N_Parent_UM(i) = N_Parent_UM(i-1) * exp(-lambda_Parent*dt);
        % Step 4.3.3 calculate the leftover isotopes in primitive mantle after the generation of CC
        N_Stable_UM(i) = N_Stable_UM(i-1) - N_Stable_UM_consumed;
        N_Daughter_UM(i) = N_Daughter_UM(i) - N_Daughter_UM_consumed;
        N_Parent_UM(i) = N_Parent_UM(i) - N_Parent_UM_consumed;
    else
        M_UM(i) = 0;
        N_Stable_UM(i) = 0;
        N_Daughter_UM(i) = 0;
        N_Parent_UM(i) = 0;
        l_UM(i) = 0;
    end   
    if (M_UM(i) == 0)
        N_Stable_UM(i) = 0;
        N_Daughter_UM(i) = 0;
        N_Parent_UM(i) = 0;
        l_UM(i) = 0;
    end
    
    % Step 4.3.4 calculate # of atoms in OM that can remelt
    N_Stable_OM(i) =  N_Stable_LM(i) + N_Stable_UM(i);
    N_Daughter_OM(i) =  N_Daughter_LM(i) + N_Daughter_UM(i);
    N_Parent_OM(i) =  N_Parent_LM(i) + N_Parent_UM(i);
    M_OM(i) = M_LM(i) + M_UM(i);
    C_Stable_LM(i) = N2C(i,N_Stable_LM(i),M_Stable,Na,M_LM(i),t(i),2);
    C_Parent_LM(i) = N2C(i,N_Parent_LM(i),M_Parent,Na,M_LM(i),t(i),2);
    C_Daughter_LM(i) = N2C(i,N_Daughter_LM(i),M_Daughter,Na,M_LM(i),t(i),2);
    C_Stable_UM(i) = N2C(i,N_Stable_UM(i),M_Stable,Na,M_UM(i),t(i),2);
    C_Parent_UM(i) = N2C(i,N_Parent_UM(i),M_Parent,Na,M_UM(i),t(i),2);
    C_Daughter_UM(i) = N2C(i,N_Daughter_UM(i),M_Daughter,Na,M_UM(i),t(i),2);
    
    % Step 4.4 Calcualte how much elements left in each remeltable DM, CR, 
    % and CRR subreservoir
    % Step 4.4.1 calculate the isotopes consumed during mantle melting
    for j = 1:(i-1)
        % DM subreservoirs
        if (l_DM_flag(j) == 1 && M_DM(j) > 0)
            N_Stable_DM_consumed(j) = N_Stable_DM(j)*F_consumed(i);
            N_Daughter_DM_consumed(j) = N_Daughter_DM(j)*F_consumed(i);
            N_Parent_DM_consumed(j) = N_Parent_DM(j)*F_consumed(i);
            M_DM(j) = M_DM(j)*(1-F_consumed(i));
            l_DM(j) = l_DM(j)*(1-F_consumed(i));
        elseif (l_DM_flag(j) == 1 && M_DM(j) == 0)
            N_Stable_DM_consumed(j) = 0;
            N_Daughter_DM_consumed(j) = 0;
            N_Parent_DM_consumed(j) = 0;
            M_DM(j) = 0;
            l_DM(j) = 0;
        else
            N_Stable_DM_consumed(j) =0;
            N_Daughter_DM_consumed(j) = 0;
            N_Parent_DM_consumed(j) = 0;
            M_DM(j) = M_DM(j);
            l_DM(j) = l_DM(j);
        end
        % CR subreservoirs
        if (l_CR_flag(j) == 1 && M_CR(j) > 0)
            N_Stable_CR_consumed(j) = N_Stable_CR(j)*F_consumed(i);
            N_Daughter_CR_consumed(j) = N_Daughter_CR(j)*F_consumed(i);
            N_Parent_CR_consumed(j) = N_Parent_CR(j)*F_consumed(i);
            M_CR(j) = M_CR(j)*(1-F_consumed(i));
            l_CR(j) = l_CR(j)*(1-F_consumed(i));
        elseif (l_CR_flag(j) == 1 && M_CR(j) == 0)
            N_Stable_CR_consumed(j) = 0;
            N_Daughter_CR_consumed(j) = 0;
            N_Parent_CR_consumed(j) = 0;
            M_CR(j) = 0;
            l_CR(j) = 0;
        else
            N_Stable_CR_consumed(j) =0;
            N_Daughter_CR_consumed(j) = 0;
            N_Parent_CR_consumed(j) = 0;
            M_CR(j) = M_CR(j);
            l_CR(j) = l_CR(j);
        end
        % RCC subreservoirs
        if (l_RCC_flag(j) == 1 && M_RCC(j) > 0)
            N_Stable_RCC_consumed(j) = N_Stable_RCC(j)*F_consumed(i);
            N_Daughter_RCC_consumed(j) = N_Daughter_RCC(j)*F_consumed(i);
            N_Parent_RCC_consumed(j) = N_Parent_RCC(j)*F_consumed(i);
            M_RCC(j) = M_RCC(j)*(1-F_consumed(i));
            l_RCC(j) = l_RCC(j)*(1-F_consumed(i));
        elseif (l_RCC_flag(j) == 1 && M_RCC(j) == 0)
            N_Stable_RCC_consumed(j) = 0;
            N_Daughter_RCC_consumed(j) = 0;
            N_Parent_RCC_consumed(j) = 0;
            M_RCC(j) = 0;
            l_RCC(j) = 0;
        else
            N_Stable_RCC_consumed(j) =0;
            N_Daughter_RCC_consumed(j) = 0;
            N_Parent_RCC_consumed(j) = 0;
            M_RCC(j) = M_RCC(j);
            l_RCC(j) = l_RCC(j);
        end
    end
    % Step 4.4.2 calculate the decay of isotopes in the sub crustal residue reservoir
    for j = 1:(i-1)
        if (M_DM(j) > 0)
            % DM subreservoirs
            N_Daughter_DM(j) = N_Daughter_DM(j) + N_Parent_DM(j)*(1-exp(-lambda_Parent*dt));
            N_Parent_DM(j) = N_Parent_DM(j)*exp(-lambda_Parent*dt);
        end
        % CR subreservoirs
        if (M_CR(j) > 0)
            N_Daughter_CR(j) = N_Daughter_CR(j) + N_Parent_CR(j)*(1-exp(-lambda_Parent*dt));
            N_Parent_CR(j) = N_Parent_CR(j)*exp(-lambda_Parent*dt);
        end
        % RCC subreservoirs
        if (M_RCC(j) > 0)
            N_Daughter_RCC(j) = N_Daughter_RCC(j) + N_Parent_RCC(j)*(1-exp(-lambda_Parent*dt));
            N_Parent_RCC(j) = N_Parent_RCC(j)*exp(-lambda_Parent*dt);
        end
    end
    % Step 4.4.3 calculate the leftover sotopes in the sub crustal residue reservoir
    for j = 1:(i-1)
        % DM subreservoirs
        if (l_DM_flag(j) == 1 && M_DM(j) > 0)
            N_Stable_DM(j) = N_Stable_DM(j) - N_Stable_DM_consumed(j);
            N_Daughter_DM(j) = N_Daughter_DM(j) - N_Daughter_DM_consumed(j);
            N_Parent_DM(j) = N_Parent_DM(j) - N_Parent_DM_consumed(j);          
        end
        % CR subreservoirs
        if (l_CR_flag(j) == 1 && M_CR(j) > 0)
            N_Stable_CR(j) = N_Stable_CR(j) - N_Stable_CR_consumed(j);
            N_Daughter_CR(j) = N_Daughter_CR(j) - N_Daughter_CR_consumed(j);
            N_Parent_CR(j) = N_Parent_CR(j) - N_Parent_CR_consumed(j);
        end
        % RCC subreservoirs
        if (l_RCC_flag(j) == 1 && M_RCC(j) > 0)
            N_Stable_RCC(j) = N_Stable_RCC(j) - N_Stable_RCC_consumed(j);
            N_Daughter_RCC(j) = N_Daughter_RCC(j) - N_Daughter_RCC_consumed(j);
            N_Parent_RCC(j) = N_Parent_RCC(j) - N_Parent_RCC_consumed(j);
        end
    end
    
    % Step 4.5 Calculate the depleted mantle (DM) isotopic concentration using mass balance
    % relationship between consumed whole mantle, generated PC and DM.
    % Equation: N_DM = N_WM_consumed - N_PC
    M_WM_consumed_total(i) = M_WM_consumed_total(i-1) + M_WM_consumed(i);
    N_Stable_WM_consumed(i) = C2N(C_Stable_WM(i-1),fraction,M_WM_consumed(i),Na,M_Stable);
    N_Stable_DM(i) = N_Stable_WM_consumed(i) - N_Stable_PC(i);
    N_Parent_WM_consumed(i) = C2N(C_Parent_WM(i-1),fraction,M_WM_consumed(i),Na,M_Parent);
    N_Parent_DM(i) = N_Parent_WM_consumed(i) - N_Parent_PC(i);
    N_Daughter_WM_consumed(i) = C2N(C_Daughter_WM(i-1),fraction,M_WM_consumed(i),Na,M_Daughter);
    N_Daughter_DM(i) = N_Daughter_WM_consumed(i) - N_Daughter_PC(i);
    if ( M_WM_consumed(i) == 0)
        N_Stable_DM(i) = 0;
        N_Parent_DM(i) = 0;
        N_Daughter_DM(i) = 0;
    end
    N_Daughter_DM_ini(i) = N_Daughter_DM(i);
    N_Stable_DM_ini(i) = N_Stable_DM(i);
    N_Parent_DM_ini(i) = N_Parent_DM(i);
    % the total # of atoms in DM
    N_Stable_DM_total(i) = sum(N_Stable_DM(1:i));
    N_Parent_DM_total(i) = sum(N_Parent_DM(1:i));
    N_Daughter_DM_total(i) = sum(N_Daughter_DM(1:i));
    M_DM_total(i) = sum(M_DM(1:i));
    % Save in matrix to tack throguh time
    N_Stable_DM_matrix(1:i,i) = N_Stable_DM(1:i);
    N_Parent_DM_matrix(1:i,i) = N_Parent_DM(1:i);
    N_Daughter_DM_matrix(1:i,i) = N_Daughter_DM(1:i);
    C_Stable_DM_matrix(1:i,i) = N2C(i,N_Stable_DM(1:i),M_Stable,Na,M_DM(1:i),t(1:i),1);
    C_Parent_DM_matrix(1:i,i) = N2C(i,N_Parent_DM(1:i),M_Parent,Na,M_DM(1:i),t(1:i),1);
    C_Daughter_DM_matrix(1:i,i) = N2C(i,N_Daughter_DM(1:i),M_Daughter,Na,M_DM(1:i),t(1:i),1);

    % Step 4.6 Calculate the continental crust (dCC) isotopic concentration 
    % Step 4.6.1 Calculate the newly generated continental crust (dCC) isotopic
    % concentration using PC, F2 and D2. 
    % Equation: C(l) = C(0)/(F+D-DF)
    F_Enrichment_Daughter = F2(i) + D_Daughter2 - D_Daughter2 * F2(i);
    F_Enrichment_Parent = F2(i) + D_Parent2 - D_Parent2 * F2(i);    
    C_Stable_dCC = C_Stable_PC / F_Enrichment_Daughter;
    [N_Stable_dCC(i)] = C2N(C_Stable_dCC,fraction,M_dCC(i),Na,M_Stable);
    C_Parent_dCC = C_Parent_PC / F_Enrichment_Parent;
    [N_Parent_dCC(i)] = C2N(C_Parent_dCC,fraction,M_dCC(i),Na,M_Parent);
    C_Daughter_dCC = C_Daughter_PC / F_Enrichment_Daughter;
    [N_Daughter_dCC(i)] = C2N(C_Daughter_dCC,fraction,M_dCC(i),Na,M_Daughter);
    N_Daughter_dCC_ini(i) = N_Daughter_dCC(i);
    N_Stable_dCC_ini(i) = N_Stable_dCC(i);
    N_Parent_dCC_ini(i) = N_Parent_dCC(i);
    % Step 4.6.2 Calculate the decay within the continental crust
    for j = 1:(i-1)
        N_Daughter_dCC(j) = N_Daughter_dCC(j) + N_Parent_dCC(j)*(1-exp(-lambda_Parent*dt));
        N_Parent_dCC(j) = N_Parent_dCC(j)*exp(-lambda_Parent*dt);
    end
    % Step 4.6.3 Calculate the total # of atoms in dCC from a step before
    N_Stable_dCC_total(i-1) = sum(N_Stable_dCC(1:i-1));
    N_Parent_dCC_total(i-1) = sum(N_Parent_dCC(1:i-1));
    N_Daughter_dCC_total(i-1) = sum(N_Daughter_dCC(1:i-1));

    % Step 4.7 Calculate crustal residue (CR) concentration using mass balance
    % Equatizon: N_CR = N_PC - N_dCC;
    N_Stable_CR(i) = N_Stable_PC(i) - N_Stable_dCC(i);
    N_Parent_CR(i) = N_Parent_PC(i) - N_Parent_dCC(i);
    N_Daughter_CR(i) = N_Daughter_PC(i) - N_Daughter_dCC(i);
    N_Daughter_CR_ini(i) = N_Daughter_CR(i);
    N_Stable_CR_ini(i) = N_Stable_CR(i);
    N_Parent_CR_ini(i) = N_Parent_CR(i);
    % the total # of atoms in CR
    N_Stable_CR_total(i) = sum(N_Stable_CR(1:i));
    N_Parent_CR_total(i) = sum(N_Parent_CR(1:i));
    N_Daughter_CR_total(i) = sum(N_Daughter_CR(1:i));   
    M_CR_total(i) = sum(M_CR(1:i));
    % Save in matrix to tack throguh time
    N_Stable_CR_matrix(1:i,i) = N_Stable_CR(1:i);
    N_Parent_CR_matrix(1:i,i) = N_Parent_CR(1:i); 
    N_Daughter_CR_matrix(1:i,i) = N_Daughter_CR(1:i);
    C_Stable_CR_matrix(1:i,i) = N2C(i,N_Stable_CR(1:i),M_Stable,Na,M_CR(1:i),t(1:i),1);
    C_Parent_CR_matrix(1:i,i) = N2C(i,N_Parent_CR(1:i),M_Parent,Na,M_CR(1:i),t(1:i),1);
    C_Daughter_CR_matrix(1:i,i) = N2C(i,N_Daughter_CR(1:i),M_Daughter,Na,M_CR(1:i),t(1:i),1);

    % Step 4.8  Calcualte the concentration of recycled crust (RCC) and its
    % effect on the concentration of continental crust (dCC)
    % Step 4.8.1 calcualte the concentraion in the RCC  at this time step
    if (M_CC(i)>0)
        N_Daughter_RCC(i) = N_Daughter_dCC_total(i-1) * M_RCC(i)/M_CC(i);
        N_Stable_RCC(i) = N_Stable_dCC_total(i-1) * M_RCC(i)/M_CC(i);
        N_Parent_RCC(i) = N_Parent_dCC_total(i-1) * M_RCC(i)/M_CC(i);
    end
    N_Daughter_RCC_ini(i) = N_Daughter_RCC(i);
    N_Stable_RCC_ini(i) = N_Stable_RCC(i);
    N_Parent_RCC_ini(i) = N_Parent_RCC(i);   
    % Step 4.8.3 calculate the isotopes left in dCC
    % calculate the atoms lost due to recycling
    for j = 1:(i-1)
        if (N_Stable_dCC_total(i-1) > 0)
            N_Stable_dCC_consumed(j) = N_Stable_dCC(j) * (N_Stable_RCC(i) / N_Stable_dCC_total(i-1));
        else
            N_Stable_dCC_consumed(j) = 0;
        end
        if (N_Parent_dCC_total(i-1)>0)
            N_Parent_dCC_consumed(j) = N_Parent_dCC(j) * (N_Parent_RCC(i) / N_Parent_dCC_total(i-1));
        else
            N_Parent_dCC_consumed(j) = 0;
        end
        if (N_Daughter_dCC_total(i-1)>0)
            N_Daughter_dCC_consumed(j) = N_Daughter_dCC(j) * (N_Daughter_RCC(i) / N_Daughter_dCC_total(i-1));
        else
            N_Daughter_dCC_consumed(j) = 0;
        end
    end
    %  calculate the amount of isotopes left in the sub-CC reservoirs
    for j = 1:(i-1)
        N_Daughter_dCC(j) = N_Daughter_dCC(j) - N_Daughter_dCC_consumed(j);
        N_Parent_dCC(j) = N_Parent_dCC(j) - N_Parent_dCC_consumed(j);
        N_Stable_dCC(j) = N_Stable_dCC(j) - N_Stable_dCC_consumed(j);      
    end
    % Step 4.9 Calculate the total # of atoms in each reservoirs
    % the total # of atoms in RCC
    N_Daughter_RCC_total(i) = sum(N_Daughter_RCC(1:i));
    N_Stable_RCC_total(i) = sum(N_Stable_RCC(1:i));
    N_Parent_RCC_total(i) = sum(N_Parent_RCC(1:i));
    M_RCC_total(i) = sum(M_RCC(1:i));
    % Save in matrix to tack throguh time
    N_Stable_RCC_matrix(1:i,i) = N_Stable_RCC(1:i);
    N_Parent_RCC_matrix(1:i,i) = N_Parent_RCC(1:i);
    N_Daughter_RCC_matrix(1:i,i) = N_Daughter_RCC(1:i);
    C_Stable_RCC_matrix(1:i,i) = N2C(i,N_Stable_RCC(1:i),M_Stable,Na,M_RCC(1:i),t(1:i),1);
    C_Parent_RCC_matrix(1:i,i) = N2C(i,N_Parent_RCC(1:i),M_Parent,Na,M_RCC(1:i),t(1:i),1);
    C_Daughter_RCC_matrix(1:i,i) = N2C(i,N_Daughter_RCC(1:i),M_Daughter,Na,M_RCC(1:i),t(1:i),1);
    % the total # of atoms in dCC
    N_Stable_dCC_total(i) = sum(N_Stable_dCC(1:i));
    N_Parent_dCC_total(i) = sum(N_Parent_dCC(1:i));
    N_Daughter_dCC_total(i) = sum(N_Daughter_dCC(1:i));
    % Save in matrix to tack throguh time
    N_Stable_dCC_matrix(1:i,i) = N_Stable_dCC(1:i);
    N_Parent_dCC_matrix(1:i,i) = N_Parent_dCC(1:i);
    N_Daughter_dCC_matrix(1:i,i) = N_Daughter_dCC(1:i);
    C_Stable_dCC_matrix(1:i,i) = N2C(i,N_Stable_dCC(1:i),M_Stable,Na,M_dCC(1:i),t(1:i),1);
    C_Parent_dCC_matrix(1:i,i) = N2C(i,N_Parent_dCC(1:i),M_Parent,Na,M_dCC(1:i),t(1:i),1);
    C_Daughter_dCC_matrix(1:i,i) = N2C(i,N_Daughter_dCC(1:i),M_Daughter,Na,M_dCC(1:i),t(1:i),1);

    % Step 4.10 Calculate mantle blobs lengths for each subreservoir (DM, CR, RCC)
    % the lengths of blobs at this time step
    if M_DM(i)>0
        l_DM(i) = Z(i);
        l_DM_flag(i) = 0;
    else
        l_DM(i) = 0;
        l_DM_flag(i) = 0;
    end
    if M_CR(i)>0
        l_CR(i) = 10e3;
        l_CR_flag(i) = 0;
    else
        l_CR(i) = 0;
        l_CR_flag(i) = 0;
    end
    if M_RCC(i)>0
        l_RCC(i) = 5e3;
        l_RCC_flag(i) = 0;
    else
        l_RCC(i) = 0;
        l_RCC_flag(i) = 0;
    end
    % the lengths of blobs before this time step
    for j = 1:(i-1)
        l_DM(j) = l_DM(j)/alpha_DM(j);
        if (l_DM(j) == 0)
            l_DM_flag(j) = 0;
        elseif (l_DM(j) <= l_critical && j < (i-t_overturn)) % (l_DM(j) <= l_critical/2 && j < (i-t_overturn)) 
            l_DM_flag(j) = 1;
        else
            l_DM_flag(j) = 0;
        end
        l_CR(j) = l_CR(j)/alpha(j);
        if (l_CR(j) == 0 )
            l_CR_flag(j) = 0;
        elseif (l_CR(j) <= l_critical && j < (i-t_overturn))
            l_CR_flag(j) = 1;
        else
            l_CR_flag(j) = 0;
        end   
        l_RCC(j) = l_RCC(j)/alpha(j);
        if (l_RCC(j) == 0)
            l_RCC_flag(j) = 0;
        elseif (l_RCC(j) <= l_critical && j < (i-t_overturn)) 
            l_RCC_flag(j) = 1;      
            
        else
            l_RCC_flag(j) = 0;
        end
    end
  
    l_DM_number(:,i) = l_DM;
    l_DM_flag_matrix(:,i) = l_DM_flag;
    l_CR_number(:,i) = l_CR;
    l_CR_flag_matrix(:,i) = l_CR_flag;
    l_RCC_number(:,i) = l_RCC;
    l_RCC_flag_matrix(:,i) = l_RCC_flag;
        
    % Step 4.11 Calculate the composition of material in the mantle that can remelt
    N_Stable_WM(i) = N_Stable_OM(i) + sum(N_Stable_DM.*l_DM_flag) + sum(N_Stable_CR.*l_CR_flag)...
        + sum(N_Stable_RCC.*l_RCC_flag);
    N_Parent_WM(i) = N_Parent_OM(i) + sum(N_Parent_DM.*l_DM_flag) + sum(N_Parent_CR.*l_CR_flag)...
        + sum(N_Parent_RCC.*l_RCC_flag);
    N_Daughter_WM(i) = N_Daughter_OM(i) + sum(N_Daughter_DM.*l_DM_flag) + sum(N_Daughter_CR.*l_CR_flag)...
        + sum(N_Daughter_RCC.*l_RCC_flag);

    M_WM(i) = M_OM(i) + sum(M_DM.*l_DM_flag) + sum(M_CR.*l_CR_flag) + ...
        sum(M_RCC.*l_RCC_flag);     
    C_Stable_WM(i) = N2C(nt,N_Stable_WM(i),M_Stable,Na,M_WM(i),t,2);
    C_Parent_WM(i) = N2C(nt,N_Parent_WM(i),M_Parent,Na,M_WM(i),t,2);
    C_Daughter_WM(i) = N2C(nt,N_Daughter_WM(i),M_Daughter,Na,M_WM(i),t,2);

    % Step 4.12 Calculate epsilon in the mantle
    R_OM(i) = N_Daughter_OM(i) / N_Stable_OM(i);
    R_WM(i) = N_Daughter_WM(i) / N_Stable_WM(i); 
    e_WM(i) = Epsilon(R_WM(i), R_OM(i),1);
    
% Step 5 Random sampling to get mantle isotopic signals
if (mod(t(i),0.25) == 0)
        k = k+1;
        M_WM_sample = M_BSE(i) - M_CC(i);
        [e_sample(k,:),M_sample(k,:),N_Stable_sample(k,:),...
            N_Daughter_sample(k,:)] = ...
            Random_Sampling_inMantle(n_sample,l_samplingbox,rhom_MM,M_WM_sample,...
            R_BSE(i),l_UM(i),l_LM(i),i,alpha,M_LM(i),N_Stable_LM(i),...
            N_Daughter_LM(i),M_UM(i),N_Stable_UM(i),N_Daughter_UM(i),...
            M_DM(1:i),l_DM(1:i),l_DM_flag(1:i),N_Stable_DM(1:i),N_Daughter_DM(1:i),...
            M_CR(1:i),l_CR(1:i),l_CR_flag(1:i),N_Stable_CR(1:i),N_Daughter_CR(1:i),...
            M_RCC(1:i),l_RCC(1:i),l_RCC_flag(1:i),N_Stable_RCC(1:i),N_Daughter_RCC(1:i),M_Stable,M_Daughter);
end   
end

% Step 6. Store the similar results in same matrix
C_BSE = [C_Stable_BSE,C_Daughter_BSE,C_Parent_BSE];
N_dCC = [N_Stable_dCC,N_Daughter_dCC,N_Parent_dCC];
N_dCC_total = [N_Stable_dCC_total,N_Daughter_dCC_total,N_Parent_dCC_total];
N_dCC_ini = [N_Stable_dCC_ini,N_Daughter_dCC_ini,N_Parent_dCC_ini];
N_DM = [N_Stable_DM,N_Daughter_DM,N_Parent_DM];
N_CR = [N_Stable_CR,N_Daughter_CR,N_Parent_CR];
N_DM_total = [N_Stable_DM_total,N_Daughter_DM_total,N_Parent_DM_total];
N_OM = [N_Stable_OM,N_Daughter_OM,N_Parent_OM];
N_CR_total = [N_Stable_CR_total,N_Daughter_CR_total,N_Parent_CR_total];
N_BSE = [N_Stable_BSE,N_Daughter_BSE,N_Parent_BSE];
N_WM = [N_Stable_WM,N_Daughter_WM,N_Parent_WM];
N_PC = [N_Stable_PC,N_Daughter_PC,N_Parent_PC];
N_RCC = [N_Stable_RCC,N_Daughter_RCC,N_Parent_RCC];
N_OM_consumed = [N_Stable_OM_consumed,N_Daughter_OM_consumed,N_Parent_OM_consumed];
N_RCC_total = [N_Stable_RCC_total,N_Daughter_RCC_total,N_Parent_RCC_total];
N_PC_ini = [N_Stable_PC_ini,N_Daughter_PC_ini,N_Parent_PC_ini]; % need to plot PC as individual
N_RCC_ini = [N_Stable_RCC_ini,N_Daughter_RCC_ini,N_Parent_RCC_ini];
N_DM_ini = [N_Stable_DM_ini,N_Daughter_DM_ini,N_Parent_DM_ini];
N_CR_ini = [N_Stable_CR_ini,N_Daughter_CR_ini,N_Parent_CR_ini];
M_Reservoir = [M_LM,M_UM,M_DM,M_CR,M_dCC,M_RCC];
M_subReservoir_total = [M_LM,M_UM,M_DM_total,M_CR_total,M_RCC_total];
N_OM_consumed_total = [N_Stable_OM_consumed_total,N_Daughter_OM_consumed_total,N_Parent_OM_consumed_total];
R_reservoir = [R_BSE,R_WM,R_OM];
N_UM = [N_Stable_UM,N_Daughter_UM,N_Parent_UM];
N_LM = [N_Stable_LM,N_Daughter_LM,N_Parent_LM];
