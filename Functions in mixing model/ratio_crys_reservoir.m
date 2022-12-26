function [R_solid_i] = ...
    ratio_crys_reservoir(R_liq_previous,xi,D_parent_i,D_daughter_i,type_crystallization)

if (type_crystallization == 1)% equilibrium crystallization
    R_solid_i = R_liq_previous * (D_parent_i/D_daughter_i)...
        *(xi *(D_daughter_i -1) +1) / (xi *(D_parent_i -1) +1);
else % fractional crystallization
    R_solid_i = R_liq_previous * ...
        (1- (1-xi)^ D_parent_i) / (1- (1-xi)^ D_daughter_i);
end