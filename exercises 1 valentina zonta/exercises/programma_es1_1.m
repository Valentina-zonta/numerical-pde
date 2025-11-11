
format long;
c=exp(-5*0.05);
c
simspon2(c,1);

a=1;
f = @(y)(-5*y);
[outputVector,output2,lasty] = runge_kutta4(f,a,0.05,6);
output2
lasty
simspon2(output2,2);


[outputVectorfe, output2fe] = forward_Euler(a);
output2fe
simspon2(output2fe,3);