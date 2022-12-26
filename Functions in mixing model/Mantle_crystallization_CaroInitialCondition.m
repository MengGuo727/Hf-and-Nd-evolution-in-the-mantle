function[Rliq_LuHf_LM_eq,Rliq_SmNd_LM_eq,Rliq_LuHf_LM_fra,Rliq_SmNd_LM_fra,...
    R_liqid_LuHf,R_solid_LuHf,R_liqid_SmNd,R_solid_SmNd,R_liqid_Sm142Nd,...
    R_solid_Sm142Nd,R_average_LuHf_LM,R_average_SmNd_LM,R_average_Sm142Nd_LM,...
    R_average_LuHf_UM,R_average_SmNd_UM,R_average_Sm142Nd_UM,...
    C_Nd_solid, C_Nd_liquid, C_Hf_solid, C_Hf_liquid,x,f] = ...
    Mantle_crystallization_CaroInitialCondition(n_crystalize,f_crystal,R_LuHf_CHUR_t0,R_SmNd_CHUR_t0,...
    R_Sm142Nd_CHUR_t0,type_crystallization,D_Sm_LM,D_Nd_LM,D_Lu_LM,D_Hf_LM,h,h_mantle,...
    C_144Nd_BSE_t0,C_177Hf_BSE_t0,type_fraction,f)

%% Calcualte the ratio of Lu/Hf and Sm/Nd in lower mantle after crystalization
length_CaPv = length(D_Lu_LM);
length_f = length(f);
Rliq_LuHf_LM_eq = nan(length_CaPv,length_f);
Rliq_SmNd_LM_eq = nan(length_CaPv,length_f); 
Rliq_LuHf_LM_fra = nan(length_CaPv,length_f);
Rliq_SmNd_LM_fra = nan(length_CaPv,length_f); 
for i = 1:length_CaPv
    for j = 1:length_f
        Rliq_LuHf_LM_eq(i,j) = ratio_residualliquid(R_LuHf_CHUR_t0,f(j),D_Lu_LM(i),D_Hf_LM(i),1);
        Rliq_SmNd_LM_eq(i,j) = ratio_residualliquid(R_SmNd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),1);
        Rliq_LuHf_LM_fra(i,j) = ratio_residualliquid(R_LuHf_CHUR_t0,f(j),D_Lu_LM(i),D_Hf_LM(i),2);
        Rliq_SmNd_LM_fra(i,j) = ratio_residualliquid(R_SmNd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),2);
    end
end


%%
n_crystalize = n_crystalize + 1;
length_D = length(D_Sm_LM);
R_liqid_LuHf = nan(n_crystalize,length_D);
R_solid_LuHf  = nan(n_crystalize,length_D);
R_liqid_SmNd = nan(n_crystalize,length_D);
R_solid_SmNd  = nan(n_crystalize,length_D);
R_liqid_Sm142Nd = nan(n_crystalize,length_D);
R_solid_Sm142Nd  = nan(n_crystalize,length_D);
R_average_LuHf_LM = nan(n_crystalize,length_D);
R_average_SmNd_LM = nan(n_crystalize,length_D);
R_average_Sm142Nd_LM = nan(n_crystalize,length_D);
R_average_LuHf_UM = nan(n_crystalize,length_D);
R_average_SmNd_UM = nan(n_crystalize,length_D);
R_average_Sm142Nd_UM = nan(n_crystalize,length_D);
R_previous_LuHf = nan(n_crystalize,length_D);
R_previous_SmNd = nan(n_crystalize,length_D);
R_previous_Sm142Nd = nan(n_crystalize,length_D);
C_Nd_solid = nan(n_crystalize,length_D);
C_Nd_liquid = nan(n_crystalize,length_D);
C_Hf_solid = nan(n_crystalize,length_D);
C_Hf_liquid = nan(n_crystalize,length_D);
x = nan(n_crystalize,1);% the mass fraction of solid during crystallization
f = nan(n_crystalize,1);% the mass fraction of liquid during crystallization, f= 1-x
C_previous_Nd = nan(n_crystalize,length_D);
C_previous_Hf = nan(n_crystalize,length_D);

% at first timestep, only have BSE, no differentiation yet
R_liqid_LuHf(1,:) = R_LuHf_CHUR_t0;
R_solid_LuHf(1,:)  = R_LuHf_CHUR_t0;
R_liqid_SmNd(1,:) = R_SmNd_CHUR_t0;
R_solid_SmNd(1,:)  = R_SmNd_CHUR_t0;
R_liqid_Sm142Nd(1,:) = R_Sm142Nd_CHUR_t0;
R_solid_Sm142Nd(1,:)  = R_Sm142Nd_CHUR_t0;
R_average_LuHf_LM(1,:) = R_LuHf_CHUR_t0;
R_average_SmNd_LM(1,:) = R_SmNd_CHUR_t0;
R_average_Sm142Nd_LM(1,:) = R_Sm142Nd_CHUR_t0;
R_average_LuHf_UM(1,:) = R_LuHf_CHUR_t0;
R_average_SmNd_UM(1,:) = R_SmNd_CHUR_t0;
R_average_Sm142Nd_UM(1,:) = R_Sm142Nd_CHUR_t0;
R_previous_LuHf(1:2,:) = R_LuHf_CHUR_t0;
R_previous_SmNd(1:2,:) = R_SmNd_CHUR_t0;
R_previous_Sm142Nd(1:2,:) = R_Sm142Nd_CHUR_t0;
C_Nd_solid(1,:) = C_144Nd_BSE_t0;
C_Nd_liquid(1,:) = C_144Nd_BSE_t0;
C_Hf_solid(1,:) = C_177Hf_BSE_t0;
C_Hf_liquid(1,:) = C_177Hf_BSE_t0;
x(1) = 0;% the mass fraction of solid during crystallization
f(1) = 0;% the mass fraction of liquid during crystallization, f= 1-x
C_previous_Nd(1:2,:) = C_144Nd_BSE_t0;
C_previous_Hf(1:2,:) = C_177Hf_BSE_t0;

for i = 2 : n_crystalize
    if (type_crystallization == 2)% fractional crystallization
        x(i) = h(i-1)*f_crystal / (h_mantle - sum(h(1:i-2)));
    else % equilibrium crystallization
        x(i) = (sum(h(1:i-1))*f_crystal) / h_mantle;
    end
    f(i) = 1 - x(i);
    for j = 1 : length_D
        % step 1. calculate the isotopic ratio in the crystalized solid
        % reservoir (LM) and the residue liquid (UM)
        R_liqid_LuHf(i,j) = ratio_residualliquid...
            (R_previous_LuHf(i,j),x(i),D_Lu_LM(j),D_Hf_LM(j),type_fraction);
        R_solid_LuHf(i,j) = ratio_crys_reservoir...
            (R_previous_LuHf(i,j),x(i),D_Lu_LM(j),D_Hf_LM(j),type_fraction);
        R_liqid_SmNd(i,j) = ratio_residualliquid...
            (R_previous_SmNd(i,j),x(i),D_Sm_LM(j),D_Nd_LM(j),type_fraction);
        R_solid_SmNd(i,j) = ratio_crys_reservoir...
            (R_previous_SmNd(i,j),x(i),D_Sm_LM(j),D_Nd_LM(j),type_fraction);
        R_liqid_Sm142Nd(i,j) = ratio_residualliquid...
            (R_previous_Sm142Nd(i,j),x(i),D_Sm_LM(j),D_Nd_LM(j),type_fraction);
        R_solid_Sm142Nd(i,j) = ratio_crys_reservoir...
            (R_previous_Sm142Nd(i,j),x(i),D_Sm_LM(j),D_Nd_LM(j),type_fraction);

        % step 2. calculate the daughter isotopes's concentration in the crystalized
        % solid reservoir (LM) and the residue liquid (UM)
        if (type_fraction == 2)% fractional crystallization
            C_Nd_liquid(i,j) = C_previous_Nd(i,j) * (f(i)^(D_Nd_LM(j)-1));
            C_Hf_liquid(i,j) = C_previous_Hf(i,j) * (f(i)^(D_Hf_LM(j)-1));
        else% Equilibrium crystallization
            C_Nd_liquid(i,j) = C_previous_Nd(i,j) / (f(i) + D_Nd_LM(j) - f(i)*D_Nd_LM(j));
            C_Hf_liquid(i,j) = C_previous_Hf(i,j) / (f(i) + D_Hf_LM(j) - f(i)*D_Hf_LM(j));
        end
        C_Nd_solid(i,j) = D_Nd_LM(j) * C_Nd_liquid(i,j);
        C_Hf_solid(i,j) = D_Hf_LM(j) * C_Hf_liquid(i,j);
        
        % step 3. calculate the weighted average isotopic ratio in the LM
        % and UM
        R_average_LuHf_LM(i,j) = Weighted_average_LMratio(i,f_crystal,...
            R_solid_LuHf(2:i,j),C_Hf_solid(2:i,j),R_liqid_LuHf(2:i,j),C_Hf_liquid(2:i,j));
        R_average_SmNd_LM(i,j) = Weighted_average_LMratio(i,f_crystal,...
            R_solid_SmNd(2:i,j),C_Nd_solid(2:i,j),R_liqid_SmNd(2:i,j),C_Nd_liquid(2:i,j));
        R_average_Sm142Nd_LM(i,j) = Weighted_average_LMratio(i,f_crystal,...
            R_solid_Sm142Nd(2:i,j),C_Nd_solid(2:i,j),R_liqid_Sm142Nd(2:i,j),C_Nd_liquid(2:i,j));

        R_average_LuHf_UM(i,j) = R_liqid_LuHf(i,j);
        R_average_SmNd_UM(i,j) = R_liqid_SmNd(i,j);
        R_average_Sm142Nd_UM(i,j) = R_liqid_Sm142Nd(i,j);
    end
    
    % step 4. rewrite the isotopic ratio in the residue liquid as the
    % system solution's inital ratio for next crystallization event
    if (type_crystallization == 2) % fractional crystallization
        R_previous_LuHf(i+1,:) = R_liqid_LuHf(i,:);
        R_previous_SmNd(i+1,:) = R_liqid_SmNd(i,:);
        R_previous_Sm142Nd(i+1,:) = R_liqid_Sm142Nd(i,:);
        C_previous_Nd(i+1,:) = C_Nd_liquid(i,:);
        C_previous_Hf(i+1,:) = C_Hf_liquid(i,:);
    else % equilibrium crystallization
        R_previous_LuHf(i+1,:) = R_LuHf_CHUR_t0;
        R_previous_SmNd(i+1,:) = R_SmNd_CHUR_t0;
        R_previous_Sm142Nd(i+1,:) = R_Sm142Nd_CHUR_t0;
        C_previous_Nd(i+1,:) = C_144Nd_BSE_t0;
        C_previous_Hf(i+1,:) = C_177Hf_BSE_t0;
    end
    
end
    

    
