function [c,f,s] = pdefun(z,t,u,dudz)
Ci=(u(1));
T=(u(2));
Tc=(u(3));
epsilon=0.5;
GHSV=10000;
L=10;            %m
rhog=1;         %Assumption of molar density 1mol/m3  %CH4
c1=0.2937e+5;                                %0.333e+5;  
c2=0.3454e+5;                               %0.7993e+5; 
c3=1.428e+3;                                 %2.87e+3; 
c4=0.264e+5;                                %0.416e+5;    
c5=588;                                  %991.96;
Cpg=c1+c2*((c3/T)/sinh(c3/T))^2+c4*((c5/T)/cosh(c5/T))^2;          %at 550K                  %https://gyazo.com/c9012c69258e2cb50730de0503c5d1a8    %https://gyazo.com/3c4f07b55f723d88a9cd04ad4435cc02
Cpg=(Cpg/1000000);   %kJ/mol K
rhos=3460;
Cps=0.836;            % at Ni/Al2O3 catalyst      10 %Ni         %kJ/kg K   
rhoc=2200;            %molten salt as coolant   kg/m3
Cpc=1.53;       %molten salt eutectic mixture 0.6 sodium nitrate 0.4 potassium nitrate    %kJ/kg K

rhocp= epsilon*rhog*Cpg+(1-epsilon)*rhos*Cps;       % equation 11 For nickel alumina catalyst
rhoccpc= rhoc* Cpc;                         % heat capacity density
c = [epsilon; rhocp; rhoccpc];          %Solving eqn 4,5,6 simultaneously

Dae= 0.00009;       %Diffusivity at 400K (Dispersion at axial coefficient)
dp=0.003;
sumCi=2;
sumCif=5;
vgf=GHSV*L/epsilon;
%vg=vgf*Ci/sumCif;                %vgf=200000 m/s GHSV=10000
vg=vgf;
mug=1.98                                %2.333e-5;
Rep=vg*rhog*dp/mug;             %gas reynold's number, eq 17/15
lambdag=0.00003757;                %0.00007857;
kae= lambdag*(8+0.05*Rep^1.09);             %eq 15
lambdac=0.000536;                  %0.000536;  
f = [Dae; kae; lambdac]; %.* dudz;

dCidz= (dudz(1));

k1=3.66*10^14 * exp(-240.1/(T*8.314*10^-3));        %rate constants
k2=5.41*10^-3 * exp(-67.1/(T*8.314*10^-3));         
k3=8.83*10^13 * exp(243.9/(T*8.314*10^-3));

delta1 = -164.9;                                     %(kJ/mol) enthalpy change at 25C
delta2 = 41.2;
delta3 = -206.1;

KCO=8.23*10^-10 * exp(70.7/(T*8.314*10^-3));        %Adsorbtion constants
KH2=6.12*10^-14 * exp(82.9/(T*8.314*10^-3));
KCH4=6.65*10^-9 * exp(38.3/(T*8.314*10^-3));
KH2O=1.77*10^5 * exp(-88.7/(T*8.314*10^-3));

K1eq=1.026676*10^10 * exp(-26830/T +30.11);         %equilibrium constant
K2eq=exp(4400/T -4.063);
K3eq=K1eq*K2eq;                 %K1/K2

P0 = 500000;                    %Pa
Ptf=-150*(((1-epsilon)^2 * mug)/(dp^2 * epsilon^3))* vg - 1.75*(((1-epsilon)*rhog)/(dp*epsilon^3))*vg^2*P0; %eq 9

kcap1= k1*rhos*(1-epsilon)/(sqrt(Ptf)*rhog*epsilon);         %eq 13
Dm=(Dae/epsilon -0.0015*vg)/(epsilon^0.5);                  %eq 14
phi1= sqrt(kcap1*dp^2/(4*Dm));
eta1=(3/phi1) * (inv(tanh(phi1))- inv(phi1));

pCO = 100000 ;                                                    %partial pressures
pH2 =400000;
pCH4 =0;
pH2O =0;

den=1+KCO*pCO+KH2*pH2+KCH4*pCH4+(KH2O*pH2O/pH2);

Ri1=(k1/pH2^2.5)*(pCH4*pH2O-(pCO*pH2^3/K1eq))/den^2;     %eq 12a

kcap2= k2*rhos*(1-epsilon)*Ptf/(rhog*epsilon);          %eq 13
phi2= sqrt(kcap2*dp^2/(4*Dm));
eta2=3/phi2 *(inv(tanh(phi2))- inv(phi2));

Ri2= (k2/pH2) *(pCO*pH2O-pH2*pCO/K2eq)/den^2;           %eq 12b

kcap3= k3*rhos*(1-epsilon)/(sqrt(Ptf)*rhog*epsilon);     % eq 13
phi3= sqrt(kcap3*dp^2/(4*Dm));
eta3=3/phi3 * (inv(tanh(phi3))- inv(phi3));

Ri3 =  (k3/pH2^3.5) *(pCH4*pH2O^2-pH2*pCO/K2eq)/den^2;           %eq 12b

sumetajRij= eta1*Ri1 + eta2*Ri2 + eta3*Ri3;
%dCidz=dudz(1);                                                     %eq 4
    s1= -epsilon*vg*dCidz+ rhos*(1-epsilon)*sumetajRij;
%
dTdz=(dudz(2));
   sumdeltaHetajRj= delta1*eta1*Ri1 + delta2*eta2*Ri2 + delta3*eta3*Ri3;
   
Nc=20;                       % minimum number of cooling tubes

hwr=(24+0.34*Rep^0.77)*lambdag/dp;      % eq 17

dw=0.002;
lambdac= 0.00007281 ;                                      %thermal conductivity of coolant
lambdaw=0.045;                                         %Steel
vc=20000;
Dc=0.02;
muc=0.0000936;                              %at 280 C
Rec=vc*rhoc*Dc/muc;                     %coolant reynold's no.
Prc=Cpc*muc/lambdac;                %coolant prandtl no.


if Rec<2030
    hwc=(3.66+(0.065*Rec*Prc*(Dc/L))/(1+0.04*((Rec*Prc*(Dc/L))^(2/3))))*lambdac/Dc;          %eq 18a
elseif 2030 < Rec < 4000
    hwc=(0.012*(Rec^0.87-280)*Prc^0.4*(1+(Dc/L)^(2/3)))*lambdac/Dc;                %eq 18b
elseif Rec>4000
    hwc=(0.027*Rec^0.8*Prc^(1/3))*lambdac/Dc;                                    %eq 18 c
end    
UwHE=inv(1/hwr+dw/lambdaw+1/hwc);                        %eq 16a
arHE=2.463;
UwHL=0.01;
arHL=20;
Te=298;
s2=-epsilon*rhog*Cpg*vg*dTdz+rhos*(1-epsilon)*sumdeltaHetajRj - Nc*UwHE*arHE*(T-Tc)- UwHL*arHL*(T-Te);          %eq 5

dTcdz= (dudz(3));
acHE=200;
    s3= -rhoc*Cpc*vc*dTcdz-UwHE*acHE*(Tc-T);                % eq 6
s = [s1; s2; s3];
end