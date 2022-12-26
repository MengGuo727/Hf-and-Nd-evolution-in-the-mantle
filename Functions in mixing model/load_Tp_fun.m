function[Tp_anchorHerz1,Tp_anchorHerz2,Tp_anchorHerz3,Tp_anchorHerz4,...
    t_anchorHerz1,t_anchorHerz2,t_anchorHerz3,t_anchorHerz4,...
    t_Herz,Tp_Herz] = load_Tp_fun(data_Tp)
% this function load in the anchor points for Herzburg's data
%
% Meng Guo, Yale University
% Summer, 2019

% column 1 is the age backward in time, column 2 is the mantle ...
% potential temperature,column 3 is the age forward in time
t_Herz = data_Tp(:,3);
Tp_Herz = data_Tp(:,2);

% chose four anchor points of Herzburg's data to make a trend
Tp_anchorHerz1 = 1470; t_Herz1back=0.8;
Tp_anchorHerz2 = 1554; t_Herz2back=1.87;
Tp_anchorHerz3 = 1578; t_Herz3back=2.75;
Tp_anchorHerz4 = 1567; t_Herz4back=3.39;
t_anchorHerz1 = 4.568-t_Herz1back;
t_anchorHerz2 = 4.568-t_Herz2back;
t_anchorHerz3 = 4.568-t_Herz3back;
t_anchorHerz4 = 4.568-t_Herz4back;
