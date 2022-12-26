%% Simple partial melting calculation
clear all;
close all;
clc;

%% calculation
C_initial = 100;
D = 0.4;
F = [0.1;0.1;0.1;0.1;0.1];


size_F = length(F);
step = 1:1:size_F+1;
C_liquid = nan(size(F));
C_solid = nan(size(F));
C_initial_array = nan(size(F));
C_liquid(1)= C_initial;
C_solid(1) = C_initial;
C_initial_array(1) = C_initial;

for i = 2:size_F+1
%     C_liquid(i) = C_liquid(i-1) / (F(i-1) + D - F(i-1)*D); % equilibrium
    C_liquid(i) = C_liquid(i-1)*F(i-1)^(D-1); % fractional
    C_solid(i) = C_liquid(i) * D;
    C_initial_array(i) = C_initial;
end

%%
figure(1);
plot(step,C_liquid,'r','LineWidth',2);hold on;
plot(step,C_solid,'b','LineWidth',2);
plot(step,C_initial_array,'g','LineWidth',2);