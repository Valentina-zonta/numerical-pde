function [ybd,CPUtimebd,num_stepsbd,infy_norm_error_bd]= BDF3(passo,matrix,y_1,y_2, matrix_width,T_final)

I=speye(matrix_width);
h=passo;
M=matrix;
N=T_final/h;

A=I+(6/11)*h*M;
y_0 = ones(matrix_width,1);
ybd(:,1)= y_0;
ybd(:,2)=y_1;
ybd(:,3)=y_2;

tStart = tic;

for i=1:N-2

b=(18/11)*ybd(:,i+2)-(9/11)*ybd(:,i+1)+(2/11)*ybd(:,i);
ybd(:,i+3)=pcg(A,b,h^3,500);

end

tEnd = toc(tStart);

CPUtimebd = tEnd;
num_stepsbd = N;


y_exact = load('accurate_solution');
error= abs(ybd(:,end)-y_exact);
infy_norm_error_bd=norm(error, Inf);
