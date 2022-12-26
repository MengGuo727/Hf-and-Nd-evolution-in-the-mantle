function[Rliq_LuHf_LM_eq,Rliq_SmNd_LM_eq,Rliq_Sm142Nd_LM_eq,Rliq_LuHf_LM_fra,...
    Rliq_SmNd_LM_fra,Rliq_Sm142Nd_LM_fra,R_LuHf_LM_eq,R_SmNd_LM_eq,R_Sm142Nd_LM_eq,...
    R_LuHf_LM_fra,R_SmNd_LM_fra,R_Sm142Nd_LM_fra,...
    C_Hf_LM_eq,C_Nd_LM_eq,C_Hf_LM_fra,C_Nd_LM_fra] = ...
    ratio_concentrration_LM(R_LuHf_CHUR_t0,f,D_Lu_LM,D_Hf_LM,R_SmNd_CHUR_t0,...
    R_Sm142Nd_CHUR_t0,D_Sm_LM,D_Nd_LM,C_Hf_BSE_t0,C_Nd_BSE_t0)

%% Calcualte the ratio of Lu/Hf and Sm/Nd in lower mantle after crystalization
length_CaPv = length(D_Lu_LM);
length_f = length(f);
Rliq_LuHf_LM_eq = nan(length_CaPv,length_f+1);
Rliq_SmNd_LM_eq = nan(length_CaPv,length_f+1); 
Rliq_Sm142Nd_LM_eq = nan(length_CaPv,length_f+1);

Rliq_LuHf_LM_fra = nan(length_CaPv,length_f+1);
Rliq_SmNd_LM_fra = nan(length_CaPv,length_f+1); 
Rliq_Sm142Nd_LM_fra = nan(length_CaPv,length_f+1); 

R_LuHf_LM_eq = nan(length_CaPv,length_f+1);
R_SmNd_LM_eq = nan(length_CaPv,length_f+1);
R_Sm142Nd_LM_eq = nan(length_CaPv,length_f+1);

R_LuHf_LM_fra = nan(length_CaPv,length_f+1);
R_SmNd_LM_fra = nan(length_CaPv,length_f+1);
R_Sm142Nd_LM_fra = nan(length_CaPv,length_f+1);

C_Hf_LM_eq = nan(length_CaPv,length_f); 
C_Nd_LM_eq = nan(length_CaPv,length_f); 
C_Hf_LM_fra = nan(length_CaPv,length_f); 
C_Nd_LM_fra = nan(length_CaPv,length_f); 


for i = 1:length_CaPv
    for j = 1:length_f
        Rliq_LuHf_LM_eq(i,j) = ratio_residualliquid(R_LuHf_CHUR_t0,f(j),D_Lu_LM(i),D_Hf_LM(i),1);
        Rliq_SmNd_LM_eq(i,j) = ratio_residualliquid(R_SmNd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),1);
        Rliq_Sm142Nd_LM_eq(i,j) = ratio_residualliquid(R_Sm142Nd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),1);
        
        Rliq_LuHf_LM_fra(i,j) = ratio_residualliquid(R_LuHf_CHUR_t0,f(j),D_Lu_LM(i),D_Hf_LM(i),2);
        Rliq_SmNd_LM_fra(i,j) = ratio_residualliquid(R_SmNd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),2);
        Rliq_Sm142Nd_LM_fra(i,j) = ratio_residualliquid(R_Sm142Nd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),2);
        
        R_LuHf_LM_eq(i,j) = ratio_crys_reservoir(R_LuHf_CHUR_t0,f(j),D_Lu_LM(i),D_Hf_LM(i),1); 
        R_SmNd_LM_eq(i,j) = ratio_crys_reservoir(R_SmNd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),1);  
        R_Sm142Nd_LM_eq(i,j) = ratio_crys_reservoir(R_Sm142Nd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),1); 
        
        R_LuHf_LM_fra(i,j) = ratio_crys_reservoir(R_LuHf_CHUR_t0,f(j),D_Lu_LM(i),D_Hf_LM(i),2);  
        R_SmNd_LM_fra(i,j) = ratio_crys_reservoir(R_SmNd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),2);            
        R_Sm142Nd_LM_fra(i,j) = ratio_crys_reservoir(R_Sm142Nd_CHUR_t0,f(j),D_Sm_LM(i),D_Nd_LM(i),2);   
        
        C_Hf_LM_eq(i,j) = (C_Hf_BSE_t0 * D_Hf_LM(i)) / ((1-f(j)) + D_Hf_LM(i) - D_Hf_LM(i)*(1-f(j)));
        C_Nd_LM_eq(i,j) = (C_Nd_BSE_t0 * D_Nd_LM(i)) / ((1-f(j)) + D_Nd_LM(i) - D_Nd_LM(i)*(1-f(j)));
        C_Hf_LM_fra(i,j) = C_Hf_BSE_t0 * D_Hf_LM(i) * ((1-f(j))^(D_Hf_LM(i)-1)); %C_Hf_BSE_t0 * (1 - (1-f(j))^(D_Hf_LM(i)-1)) / f(j);
        C_Nd_LM_fra(i,j) = C_Nd_BSE_t0 * D_Nd_LM(i) * ((1-f(j))^(D_Nd_LM(i)-1)); %C_Nd_BSE_t0 * (1 - (1-f(j))^(D_Nd_LM(i)-1)) / f(j);

    end
end

Rliq_LuHf_LM_eq(:,length_f+1) = R_LuHf_CHUR_t0;
Rliq_SmNd_LM_eq(:,length_f+1) = R_SmNd_CHUR_t0; 
Rliq_Sm142Nd_LM_eq(:,length_f+1) = R_Sm142Nd_CHUR_t0;

Rliq_LuHf_LM_fra(:,length_f+1) = R_LuHf_CHUR_t0;
Rliq_SmNd_LM_fra(:,length_f+1) = R_SmNd_CHUR_t0; 
Rliq_Sm142Nd_LM_fra(:,length_f+1) = R_Sm142Nd_CHUR_t0; 

R_LuHf_LM_eq(:,length_f+1) = R_LuHf_CHUR_t0; 
R_SmNd_LM_eq(:,length_f+1) = R_SmNd_CHUR_t0;
R_Sm142Nd_LM_eq(:,length_f+1) = R_Sm142Nd_CHUR_t0; 

R_LuHf_LM_fra(:,length_f+1) = R_LuHf_CHUR_t0;
R_SmNd_LM_fra(:,length_f+1) = R_SmNd_CHUR_t0;
R_Sm142Nd_LM_fra(:,length_f+1) = R_Sm142Nd_CHUR_t0;
