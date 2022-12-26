function [F_Jun_same,S_Jun_same] ...
    = load_FandS_fun(t,nt,data_formationage,data_zircon_surf)
% this function is used for load in data of observed formation age and
% surface age distributions
%
% Meng Guo, Yale University
% Summer, 2019


% data for formation age distribution (Korenaga, 2018)
t_Jun = data_formationage(:,1);
t_Jun = 4.568-(t_Jun/1000);
F_Jun = data_formationage(:,2);

% data for surface age distribution (Roberts & Spencer, 2015)
S_Jun = data_zircon_surf(:,2);

% now need to make them same demonsion as our time series
% formation age distribution (Korenaga, 2018)
n_Jun=length(F_Jun);
Ratio=floor(nt/n_Jun);
F_Jun_same=interp(F_Jun,Ratio);
Diff=nt-length(F_Jun_same);
F_Jun_same=[F_Jun_same;zeros(Diff,1)];
F_Jun_same = flipud(F_Jun_same);

% surface age distribution (Roberts & Spencer, 2015)
n_zircon=length(S_Jun);
Ratio2=floor(nt/n_zircon);
S_Jun_same=interp(S_Jun,Ratio2);
Diff2=nt-length(S_Jun_same);
S_Jun_same=[S_Jun_same;zeros(Diff2,1)];
t_zs_same = 4.567-t;
S_Jun_same = flipud(S_Jun_same);
t_zs_same = flipud(t_zs_same);
