function [ xi,yi ] = EOSmmp( zi )

%**************************************************************************
% 1.- C1, 2.- C4, 3.- C10 4.- CO2
%Imput data
R = 8.31438;%Pa-m^3/mol-K
P = 2000;%psia
P = P*3892.8571;%Pa
T = 160;%F
T = (T-32)/1.8+273.15;%K


%Compositions of C1, C4, C10 and CO2
xoil = [0.2 0.15 0.65 0];% Oil
ygas = [0.2 0 0 0.8];%Gas
%Properties
Tc = [-116.63 305.65 652.1 87.9];%F
Tc = ((Tc-32)./1.8)+273.15;%K
Pc = [667.8 550.7 305.7 1071];%psia
Pc = Pc.*3892.8571;%Pa
wi = [.0104 .201 .49 .225];%Ascentric factor
Tr = T./Tc;

%VAlues for initialization

E = [1e-6 1e-6 1e-6 1e-6];
dKi = 1000;
absdKi = abs(dKi);
Ki = (Pc./P).*exp(5.37.*(1+wi).*(1-(Tc./T)));



while(absdKi>=E);
% Rachford-Rice Equation 
F0 = sum(zi.*(1-Ki)./Ki);
F1 = 1-sum(zi.*Ki);
if F0>0&&F1<0
    lmin = 0;
    lmax = 1;
elseif F0>0&&F1>0
    lmin = 1;
    lmax = max((Ki.*zi-Ki)./(1-Ki));
elseif F0<0&&F1<0
    lmin = min((zi-Ki)./(1-Ki));
    lmax = 0;
end
if min(Ki)<1
    [Kimin,indexmin] = min(Ki);
    lmin = (zi(indexmin)-Kimin)/(1-Kimin);
end
if max(Ki)>1
    [Kimax,indexmax] = max(Ki);
    lmax = (Kimax*zi(indexmax)-Kimax)/(1-Kimax);
end
l = (lmin+lmax)/2;
F = sum(zi.*(1-Ki)./(Ki+(1-Ki).*l));
while(abs(F)>1e-4)
    if F>0
        lmin = l;
    elseif F<0
        lmax = l;
    end
    af = (1-Ki)./((1-Ki).*l+Ki);
    df = -sum(zi.*af.^2);
    l = l-F/df;
    if (l>=lmax)||(l<=lmin)
        l = (lmin+lmax)/2;
    end
    F = sum(zi.*(1-Ki)./(Ki+(1-Ki).*l));
end

xi = zi./(l+Ki.*(1-l));
yi = Ki.*zi./(l+Ki.*(1-l));

%PENG-ROBINSON (PR)  EOS
u = 2;
w = -1;
kij = zeros(length(wi),length(wi));
%Values from the Table 1 of SPE - 116823 paper.
kij(2,1) = 0.027;
kij(3,1) = 0.042;
kij(4,1) = 0.1;
kij(3,2) = 0.008;
kij(4,2) = 0.1257;
kij(4,3) = 0.0942;
bi = .0778*R.*Tc./Pc;
fwi = .37464+1.54226.*wi-.26992.*(wi.^2);
ai = .45724*(R^2).*(Tc.^2)./Pc;
ai = ai.*(1+fwi.*(1-Tr.^.5)).^2;

% For Liquid
for i=1:length(ai);
aL(i) = sum(xi(i).*xi.*((ai(i).*ai).^.5).*(1-kij(i,:)));
end
aL =sum(aL);
bL = sum(xi.*bi);
AL = aL*P/((R^2)*(T^2));
BL = bL*P/(R*T);
a = 1;
b = -(1+BL-u*BL);
c = AL+w*BL^2-u*BL-u*BL^2;
d = -AL*BL-w*BL^2-w*BL^3;
EqL = [a b c d];
solnPRL = roots(EqL);
irootL = imag(solnPRL);
[indeximL] = find(irootL==0);
ZL = solnPRL(indeximL);
%Volume
VL = ZL*R*T/P;

%Fugacity
for i=1:length(ai);
deltaLi1(i) = sum(xi.*(ai.^.5).*(1-kij(i,:)));
deltaLi(i) = (2.*(ai(i).^.5)./aL)*deltaLi1(i);
end
uu = (u^2-4*w)^.5;
bi_bL = (Tc./Pc)./sum(xi.*Tc./Pc);
coefugLi = bi_bL.*(ZL-1)-log(ZL-BL);
coefugLi1 = AL/(BL*(u^2-4*w)^.5).*(bi_bL-deltaLi);
coefugLi2 = 2*ZL+BL*(u+uu);
coefugLi3 = 2*ZL+BL*(u-uu);
coefugLi = coefugLi + coefugLi1.*log(coefugLi2/coefugLi3);
coefugLi = exp(coefugLi);
fugLi = P.*coefugLi.*xi;

% For Vapor
for i=1:length(ai);
aV(i) = sum(yi(i).*yi.*((ai(i).*ai).^.5).*(1-kij(i,:)));
end
aV =sum(aV);
bV = sum(yi.*bi);
AV = aV*P/((R^2)*(T^2));
BV = bV*P/(R*T);
a = 1;
b = -(1+BV-u*BV);
c = AV+w*BV^2-u*BV-u*BV^2;
d = -AV*BV-w*BV^2-w*BV^3;
EqV = [a b c d];
solnPRV = roots(EqV);
irootV = imag(solnPRV);
[indeximV] = find(irootV==0);
ZV = solnPRV(indeximV);
%Volume
VV = ZV*R*T/P;

%Fugacity
for i=1:length(ai);
deltaVi1(i) = sum(yi.*(ai.^.5).*(1-kij(i,:)));
deltaVi(i) = (2.*(ai(i).^.5)./aV)*deltaVi1(i);
end
bi_bV = (Tc./Pc)./sum(yi.*Tc./Pc);
coefugVi = bi_bV.*(ZV-1)-log(ZV-BV);
coefugVi1 = AV/(BV*(u^2-4*w)^.5).*(bi_bV-deltaVi);
coefugVi2 = 2*ZV+BV*(u+uu);
coefugVi3 = 2*ZV+BV*(u-uu);
coefugVi = coefugVi + coefugVi1.*log(coefugVi2/coefugVi3);
coefugVi = exp(coefugVi);
fugVi = P.*coefugVi.*yi;

if(fugLi~=fugVi)
Ki = (fugLi./fugVi).*Ki;
dKi = (fugLi./fugVi)-1;
absdKi = abs(dKi);
end

end
v = 1-l;
%*************************************************************************

end

