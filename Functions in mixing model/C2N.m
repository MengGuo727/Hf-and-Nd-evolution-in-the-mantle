 function [N_atoms] = C2N(concentration,fraction,M_reservor,Na,M_atomic)
% This function calculates number of atoms for a certain element or isotope
% 
% Meng Guo, Yale University
% Summer 2020

% Input variables and their meanings:
% concentration: concentration of an element, unitless
% fraction: the fraction of an isotope to its element abundance
% M_reservor: mass of the reservor in interest, in unit kg
% Na: Avogadro's constant
% M_atomic: the atomic mass of the element or isotope in interest

N_atoms = fraction * concentration .* M_reservor * Na *1e3 / M_atomic; % in unit # of atoms
