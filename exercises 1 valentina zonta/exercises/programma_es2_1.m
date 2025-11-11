clear
format long;
%define the Cauchy problem
a=1;
f = @(y)(-10*(y^2));
 %implement runge kutta 4 order and find y(last t)
lastk= zeros;
for i=5:10
    [outputVector,output2, lasty] = runge_kutta4(f,a,2^(-i), 2);
     lastk(i)=lasty;
end    
 
%evaluate the real y(t=2) 
realy=@(t) (1/(1+10*t));
last_real_y=ones(1,10)*realy(2);

%finding errors
error=abs(last_real_y-lastk);
errors = error(5:10);
N=[2^(6) 2^(7) 2^(8) 2^(9) 2^(10) 2^(11)]


loglog(N, errors, '-o')
title("loglog scale of errors on y(T), T=2 and h=2^k k=5..10")
xlabel("numeber of steps")
ylabel("error")

h=[2^-5 2^-6 2^-7 2^-8 2^-9 2^-10];
errors=errors'
h=h';
T=table(h,errors)


