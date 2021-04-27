function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
epsilon=0.5;
vgf=200000;              %All velocity : gas fluid, coolant etc . are 200000 m/s despite their position in the reactor
Cif=1;
rhog=1;
T=(ul(2));                                  %CH4
c1=0.2937e+5;                                %0.333e+5;  
c2=0.3454e+5;                               %0.7993e+5; 
c3=1.428e+3;                                 %2.87e+3; 
c4=0.264e+5;                                %0.416e+5;    
c5=588;                                  %991.96;
Cpg=c1+c2*((c3/T)/sinh(c3/T))^2+c4*((c5/T)/cosh(c5/T))^2;          %at 550K                  %https://gyazo.com/c9012c69258e2cb50730de0503c5d1a8    %https://gyazo.com/3c4f07b55f723d88a9cd04ad4435cc02
Cpg=(Cpg/1000000);
dp=0.003;
mug=2.333e-5;
vg=200000;
Tf=550;
rhoc=2200;
vcf=20000;
Dae= 0.00009; 
kc=5;
lambdac=0.00007281;
Cc=1;
Tcf=415;


pl1=epsilon*vgf*(Cif-(ul(1)));
pl2=epsilon*rhog*vgf*Cpg*(Tf-(ul(2)));
pl3=rhoc*vcf*Cc*(Tcf-(ul(3)));                    % eq 7

pl = [pl1; pl2; pl3];                           %z=0 BC
ql = [1; 1; kc/lambdac];

pr = [0; 0; 0];                                 %z=L BC
qr = [1; 1; 1];
end