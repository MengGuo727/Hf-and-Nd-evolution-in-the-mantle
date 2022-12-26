function[Qc_model,Qc_backward] = Qc_backward_fun(Qc_tp_model,d_Qc_model,nt,t)
% this function rearrange the Qc to be backward order
%
% Meng Guo, Yale University
% Summer, 2019

% heat flow from core, in unit TW
Qc_backward = nan(size(t));
% present-day heat flow from core, in unit TW
Qc_backward(1) = Qc_tp_model;

% Qc linearly increase from 10TW to 15TW, backward in time
for i = 2:nt
    Qc_backward(i) = (d_Qc_model/4.567)*t(i) + Qc_tp_model;
end
Qc_model = flipud(Qc_backward);