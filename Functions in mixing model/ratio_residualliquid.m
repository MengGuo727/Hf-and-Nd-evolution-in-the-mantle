function [R_liq_i] = ...
    ratio_residualliquid(R_liq_previous,xi,D_parent_i,D_daughter_i,type_crystallization)

if (type_crystallization == 1)% equlibrium crystallization
    R_liq_i = R_liq_previous * ...
        (xi *(D_daughter_i -1) +1) / (xi *(D_parent_i -1) +1);
else % fractional crystallization
    R_liq_i = R_liq_previous * ...
        (1-xi)^(D_parent_i -1) / (1-xi)^(D_daughter_i -1);
end