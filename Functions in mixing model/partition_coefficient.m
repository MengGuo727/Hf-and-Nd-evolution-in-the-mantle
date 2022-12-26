function[D] = partition_coefficient(a_mineralA,D_mineralA,a_mineralB,D_mineralB)

length_CaPv = length(a_mineralA);
D = nan(length_CaPv,1);
for i = 1:length_CaPv
    D(i) = a_mineralA(i)*D_mineralA + a_mineralB(i)*D_mineralB; 
end