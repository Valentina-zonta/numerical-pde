
clc
clear

nx=[102 202 402 802];

for i=1:4
A=delsq(numgrid('S',nx(i)));
h=1/(nx(i)-2);
i_values = 1:size(A,1);
i_values = i_values(:);
x_exact = 1 ./ sqrt(i_values);
b=A*x_exact;

L_1=ichol(A);

options.type = 'ict';
options.droptol = 1e-2;
L_2 = ichol(A, options);
         
option.type = 'ict';
option.droptol = 1e-3;
L_3=ichol(A,option);

[xa, flaga, relresa, itera(i)] = pcg(A, b, 10^(-8), 5000);
[xb, flagb, relresb, iterb(i)] = pcg(A, b, 10^(-8), 5000,L_1,L_1');
[xc, flagc, relresc, iterc(i)] = pcg(A,b,1e-8,5000,L_2,L_2');
[xd, flagd, relresd, iterd(i)] = pcg(A, b, 10^(-8), 5000,L_3,L_3');

end

%%table on latex