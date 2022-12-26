function[R_average] = ...
    Weighted_average_LMratio(i,f_crystal,R_solid,C_solid,R_liqid,C_liquid)

% R_average = sum(R_solid .* (f_crystal/i) .* C_solid / ...
%     sum(f_crystal/i * C_solid + C_liquid * (1-f_crystal)/i) ...
%     + R_liqid .* C_liquid * (1-f_crystal)/i / ...
%     sum(f_crystal/i * C_solid + C_liquid * (1-f_crystal)/i));

R_average = ( sum(R_solid .* C_solid * f_crystal/(i-1)) ...
    + sum(R_liqid .* C_liquid * (1-f_crystal)/(i-1)))...
    / (sum(C_solid*f_crystal/(i-1)) + sum(C_liquid)*(1-f_crystal)/(i-1));