function [concentration] = N2C(nt,Numberofatoms,M_atomic,Na,M_reservoir,t,type)
% This function calculates number of atoms for a certain element or isotope
% 
% Meng Guo, Yale University
% Summer 2020

% Input variables and their meanings:
% nt: # of total timesteps
% Numberofatoms: # of atoms of a certain isotope
% M_reservor: mass of the reservor in interest, in unit kg
% Na: Avogadro's constant
% M_atomic: the atomic mass of the element or isotope in interest
% t: time serie
% type: 1 means the output is a series of values; 2 mean the output is a single value
if (type ==1)
    concentration = nan(size(t));
    for i = 1:nt
        if (M_reservoir(i)==0)
            concentration(i)=0;
        else
            concentration(i) = Numberofatoms(i) * M_atomic / Na / 1e3 /M_reservoir(i);
        end
    end
else
    concentration = Numberofatoms * M_atomic / Na / 1e3 /M_reservoir;
end