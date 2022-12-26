% This is the function to randomly sample in the mantle and calculate the
% corresponding geochemical signals
function [Signal_sample,M_sample,N_Stable_sample,...
    N_Daughter_sample,index_reservoir,Signal_choose,l_s_choose,Signal_cube,...
    C_Daughter_cube,C_Stable_cube,C_Daughter_reservoir] = ...
    Random_Sampling_inMantle(n_sample,l_samplingbox,rhom_MM,M_WM_sample,...
    R_BSE,l_UM_i,l_LM_i,i,alpha,M_LM,N_Stable_LM,...
    N_Daughter_LM,M_UM,N_Stable_UM,N_Daughter_UM,...
    M_DM,l_DM_i,l_DM_flag,N_Stable_DM,N_Daughter_DM,...
    M_CR,l_CR_i,l_CR_flag,N_Stable_CR,N_Daughter_CR,...
    M_RCC,l_RCC_i,l_RCC_flag,N_Stable_RCC,N_Daughter_RCC,M_Stable,M_Daughter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n_sample: the number of samples we take at each sample time
% l_samplingbox: the length of sampling box
% rhom_MM: the mantle density
% M_reservoir: the matrix of subreservoirs' masses in the mantle
% l_reservoir: the matrix of subreservoirs' length scale in the mantle
% N_Stable_reservoir: the matrix of subreservoirs' # of Stable in the mantle
% N_Parent_reservoir: the matrix of subreservoirs' # of Parent in the mantle
% N_Daughter_reservoir: the matrix of subreservoirs' # of Daughter in the mantle
% M_WM_sample: the mass of the manlte that can be sampled (remeltable mantle)
% R_BSE: the ratio of the daughter to stable isotope in the BSE
% l_UM_t0: the initial length of upper mantle reservoir
% l_LM_t0: the initial length of lower mantle reservoir
% i: the timestep at which we are doing random sampling in the mantle
% alpha: stretching factor
% M_reservoir: the mass of the reservoir in interest
% N_Stable_reservoir: the # of stable isotope of the reservoir in interest
% N_Parent_reservoir: the # of parent isotope of the reservoir in interest
% N_Daughter_reservoir: the # of daughter isotope of the reservoir in interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Streching LM and UM according to stretching factor
min_length = 10;
if (l_UM_i < min_length)
    l_UM_i = min_length;
end
if (l_LM_i < min_length)
    l_LM_i = min_length;
end
range = logical(l_DM_i<min_length);
l_DM_i(range) = min_length;
range = logical(l_CR_i<min_length);
l_CR_i(range) = min_length;
range = logical(l_RCC_i<min_length);
l_RCC_i(range) = min_length;

Na = 6.022E+23; % Avogadro's constant

% Each of the thin laminar has the same concentration of the original LM and UM
M_reservoir = [M_LM;M_UM;M_RCC.*l_RCC_flag; M_DM.*l_DM_flag; M_CR.*l_CR_flag]; M_reservoir = nonzeros(M_reservoir);
l_reservoir = [l_LM_i;l_UM_i;l_RCC_i.*l_RCC_flag; l_DM_i.*l_DM_flag; l_CR_i.*l_CR_flag]; l_reservoir = nonzeros(l_reservoir);
N_Stable_reservoir = [N_Stable_LM; N_Stable_UM;N_Stable_RCC.*l_RCC_flag; N_Stable_DM.*l_DM_flag; N_Stable_CR.*l_CR_flag];N_Stable_reservoir = nonzeros(N_Stable_reservoir);
N_Daughter_reservoir = [N_Daughter_LM; N_Daughter_UM;N_Daughter_RCC.*l_RCC_flag; N_Daughter_DM.*l_DM_flag; N_Daughter_CR.*l_CR_flag;];
N_Daughter_reservoir = nonzeros(N_Daughter_reservoir);
C_Daughter_reservoir = N_Daughter_reservoir ./ M_reservoir;

% Creates empty matrixes to store the results
Signal_sample = nan(n_sample,1); % the geochemical signals in the mantle
l_sample = nan(n_sample,1); % the total mass of the material being sampled
M_sample = nan(n_sample,1); % the total mass of the material being sampled
N_Stable_sample = nan(n_sample,1); % the total # of Stable being sampled
N_Daughter_sample = nan(n_sample,1); % the total # of Daughter being sampled

rng('shuffle');% to avoid have same "rand" results everytime
N_MAX = round(l_samplingbox/min(l_reservoir)) + 5;
index_reservoir = nan(N_MAX,n_sample);

z = 0;
P = M_reservoir / sum(M_reservoir);
population = 1:1:length(M_reservoir); population = population';

Signal_cube = nan(N_MAX,z);
Signal_choose = nan(N_MAX,z);
l_s_choose = nan(N_MAX,z);
C_Daughter_cube = nan(N_MAX,z);
C_Stable_cube = nan(N_MAX,z);

% Sample n_sample times
while (z < n_sample) 
    z = z + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fill the sampling box
    % make sure to start with empty sampling box
    M_sample(z,1) = 0;
    l_sample(z,1) = 0;
    N_Stable_sample(z,1) = 0;
    N_Daughter_sample(z,1) = 0;
    n_sample_repeat = 0;
    
    l_s = l_samplingbox; % the length of the unfilled sampling box
    l_s_origin = l_samplingbox;
    
    index = randsample(population,N_MAX,true,P); % choose one subreservoir according to their probability
    
    index_single = index;
    range = logical(index>2 & index<=2+length(nonzeros(l_RCC_flag)));
    index_single(range) = 3;
    range = logical(index>2+length(nonzeros(l_RCC_flag)) & index<=2+length(nonzeros(l_RCC_flag))+length(nonzeros(l_DM_flag)));
    index_single(range) = 4;
    range = logical(index> 2+length(nonzeros(l_RCC_flag))+length(nonzeros(l_DM_flag)) & index<=2 + length(nonzeros(l_RCC_flag))+length(nonzeros(l_DM_flag))+length(nonzeros(l_CR_flag)));
    index_single(range) = 5;
    index_reservoir(1:N_MAX,z) = index_single;

    
    if mod(z,1) == 0
        disp(['z=' num2str(z)]);% keep track of the calculation
    end % if mod(z,10) == 0
    
    % Random sampling
    % Step 1. fill the depleted proportion of the sampling box
    while l_s > 0 % as long as the sampling box is not completly filled, the random sampling continues
        n_sample_repeat = n_sample_repeat + 1; % the # of sample drawed to fill the sampling box
        % Step 1. Choose a subreservoir to sample
        M_cube = M_reservoir(index(n_sample_repeat)); % the mass of the chosen subreservoir
        N_Stable_cube =  N_Stable_reservoir(index(n_sample_repeat));
        N_Daughter_cube =  N_Daughter_reservoir(index(n_sample_repeat));
        l_cube = l_reservoir(index(n_sample_repeat)); % the length of the chosen subreservoir
        
        % Step 2. Fill the sampling box with the chosen subreservoir
        if (l_cube > l_s)
            f = l_s/l_cube;
        else
            f = 1;
        end
        l_choose = f * l_cube;
        M_choose = f * M_cube;
        N_Stable_choose  = N_Stable_cube/M_cube*l_choose/l_s_origin; % save the concentration of cube I choose
        N_Daughter_choose  = N_Daughter_cube/M_cube*l_choose/l_s_origin;  
        R_cube = N_Daughter_choose / N_Stable_choose;
        Signal_cube(n_sample_repeat,z) =  Epsilon(R_cube, R_BSE,2);
        C_Daughter_cube(n_sample_repeat,z) = N_Daughter_cube* M_Daughter/ Na / 1e3 / M_cube;
        C_Stable_cube(n_sample_repeat,z) = N_Stable_cube* M_Stable/ Na / 1e3 / M_cube;
        
        % fill the sampling box with the chosen material
        l_sample(z,1) = l_sample(z,1) + l_choose;
        M_sample(z,1) = M_sample(z,1) + M_choose;
        N_Stable_sample(z,1) = N_Stable_sample(z,1) + N_Stable_choose;
        N_Daughter_sample(z,1) = N_Daughter_sample(z,1) + N_Daughter_choose;       
        R_choose = N_Daughter_sample(z,1) / N_Stable_sample(z,1);
        Signal_choose(n_sample_repeat,z) =  Epsilon(R_choose, R_BSE,2);

        % Step. 3 Determine if the sampling box is filled
        l_s = l_s_origin - l_sample(z,1);
        l_s_choose(n_sample_repeat,z) = l_s_origin - l_s;
   end % while l_s(z) > 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the total sampling box's signal
    Signal_sample(z,1) =  Signal_choose(n_sample_repeat,z);
end % for z = 1 : n_sample


