%symbolic algebra to find optimal kc for full analytical a solution

syms U V W F X a p
syms VPQ VCDE Rc D kmP kmA Rb KPQ KCDE Aout Pout kc jc

X=D/(Rc^2*kc)+1/Rc-1/Rb;
Y=(VCDE*Rc^3)*(D/(kmP*Rb^2)+X)/(3*D*KCDE);
Z=(Pout*(1+jc/kmP))/KCDE;
E=Y-Z+1;
p=(-E+sqrt(E^2+4*Z))/2;
U=(2*VPQ*Rc^3)*(D/(kmP*Rb^2)+X)/(3*D*KPQ);
V=Aout/KPQ;
W=1/2*VCDE/VPQ;
F=1+U-V;

a=(-(F-U*W*p/(1+p))+sqrt((F-U*W*p/(1+p))^2+4*(U*W*p/(1+p)+V)))/2