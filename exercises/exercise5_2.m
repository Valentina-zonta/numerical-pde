clc
clear

D = readmatrix('mat13041.rig.txt');
A = spconvert(D);
i_values = 1:size(A,1);
i_values = i_values(:);
x_exact = 1 ./ sqrt(i_values);
b=A*x_exact;

x0=zeros(size(A,1),1);

[x, iter, resvec, flag] = mygmres(x0, A, b, 550, 10^(-10));



