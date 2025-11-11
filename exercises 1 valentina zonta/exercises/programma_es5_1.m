clear
%matrix A of the problem
nx = 100;
G = numgrid ('S' , nx ) ;
A = delsq (G) * ( nx-1)^2 ;


%largest eigenvalue of A

lamda = -eigs(A,1,'lm');

% solve inequality for 4 order runge kutta absolute r5egion
min=-2.78;
max=0;

min_lamda=-2.78/lamda;
max_lamda=0;


%solve the problem by ode
count_columns=width(A);
tspan = [0 0.1];
y_0 = ones(count_columns,1);

tic
[t,y] = ode45(@(t,y) -A*y, tspan, y_0);
toc

% find the error in the last time (last row of the solution)
y_exact = load('accurate_solution');
error= y(end,:)'-y_exact;
infy_norm_error=norm(error, Inf);


% Crank Nicolson (CN) method with h ∈ {10−3, 10−4, 10−5}
%{ 

I=speye(count_columns);
h=10^(-3);

N=0.1/h;

M=I+(h/2)*A;
B=I-(h/2)*A;
y_0 = ones(count_columns,1);

ykn(:,1)=y_0;


tStart = tic;

for i=1:N

b=B*ykn(:,i);
ykn(:,i+1)=pcg(M,b,h^3,500);

end
tEnd = toc(tStart);
CPUtime = tEnd;
num_steps = N;

error= abs(ykn(:,end)-y_exact);
infy_norm_error_cn=norm(error, Inf);


[ybd,CPUtimebd,num_stepsbd,infy_norm_error_bd]=BDF3(h,A,ykn(:,2),ykn(:,3), count_columns,0.1);

%}




h=[10^(-3) 10^(-4) 10^(-5)];

for i=1:length(h)
[y_kn,CPUtime,num_steps,norm_error]=crank_nic_pcg(count_columns,h(i),A,0.1);
CPUtime_cn_h(i)=CPUtime;
num_steps_cn_h(i)=num_steps;
infy_norm_error_cn_h(i)=norm_error;



%solve bdf3

[y_bd,CPUtime,num_steps,norm_error]=BDF3(h(i),A,y_kn(:,2),y_kn(:,3), count_columns,0.1);
CPUtime_b3_h(i)=CPUtime;
num_steps_b3_h(i)=num_steps;
infy_norm_error_b3_h(i)=norm_error;
y_b3=y_bd;
end

