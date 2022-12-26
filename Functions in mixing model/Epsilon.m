 function [epsilon] = Epsilon(ratio_sample, ratio_reference,mode)
% This function calculates epsilon of the isotope in interest
% 
% Meng Guo, Yale University
% Summer 2020

% Input variables and their meanings:
% ratio_sample: isotope ratio in the sample
% ratio_reference: isotipe ratio in the reference reservoir (CHUR or BSE etc)
% mode: determines how to calculate epsilon
if (mode == 2)
    epsilon = zeros(size(ratio_sample));
    [length,width] = size(ratio_sample);
    for i = 1:length
        if (ratio_sample(i) == 0)
            epsilon(i) = 0;
        else
            epsilon(i) = (ratio_sample(i) - ratio_reference)/ratio_reference * 1e4;     
        end
    end
else
    epsilon = zeros(size(ratio_sample));
    [length,width] = size(ratio_sample);
    for i = 1:length
        if (ratio_sample(i) == 0)
            epsilon(i) = 0;
        else
            epsilon(i) = (ratio_sample(i) - ratio_reference(i))/ratio_reference(i)* 1e4;
        end
    end
end
