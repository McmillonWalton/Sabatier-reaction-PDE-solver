clc;
clear;
m=1;
z = linspace(0,1,5);
t = 0:1:5;
sol=abs(pdepe(m,@pdefun,@pdeic,@pdebc,z,t));