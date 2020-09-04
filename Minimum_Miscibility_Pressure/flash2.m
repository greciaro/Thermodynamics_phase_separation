function [xi,yi,Ki]=flash2(zi) 
  %clear all
  % INPUT: Total composition.
% OUTPUT: Liquid and Vapor composition and equilibrium constants.

% 1.- Calculating initial values and program parameters.
R = 8.31438;%Pa-m^3/mol-K
P = 2000;%psia
P = P*6892.8571;%Pa
T = 160;%F
T = (T-32)/1.8+273.15;%K
Tc = [-116.63 305.65 652.1 87.9];%F
Tc = ((Tc-32)./1.8)+273.15;%K
Pc = [667.8 550.7 305.7 1071];%psia
Pc = Pc.*6892.8571;%Pa
wi = [.0104 .201 .49 .225];%Ascentric factor.
Tr = T./Tc;%Reduced temperatue
E = [1e-6 1e-6 1e-6 1e-6];%Tolerance to reach Ki for equilibrium.
dKi = 1000;%intial value to compare tolerance.
absdKi = abs(dKi);%intial value to compare tolerance.
Ki = (Pc./P).*exp(5.37.*(1+wi).*(1-(Tc./T)));%initial value of Ki with Wilson equation.
alfa = 0.5;%proportion to mix liquid and vapor and get zi.
%%
%****************************F L A S H**********************************
iter=1; % To know how many iterations happen before to rerach equilibrium.
while(absdKi>=E);%Compares if the difference of previous and new Ki is inside of the tolerance.
    
%   2.- Solving Rachford-Rice equation to find “l”, “xi” and “yi”.
    
    %Setting the right and left boundaries of the search for "l".
    F0 = sum(zi.*(1-Ki)./Ki);
    F1 = 1-sum(zi.*Ki);
    if F0>0 && F1<0
        lmin = 0;
        lmax = 1;
    elseif F0>0 && F1>0
        lmin = 1;
        lmax = max((Ki.*zi-Ki)./(1-Ki));
    elseif F0<0 && F1<0
        lmin = min((zi-Ki)./(1-Ki));
        lmax = 0;
    end
    l = (lmin+lmax)/2;%Calculating "l" by bisection.
    
    %Checking for tighter boundaries
    [Kimin,indexmin] = min(Ki);
    if Kimin>1%l<(zi(indexmin)-Kimin)./(1-Kimin)
        lmin = (zi(indexmin)-Kimin)./(1-Kimin);
    end
    [Kimax,indexmax] = max(Ki);
    if Kimax<1%l>((Kimax.*zi(indexmax)-Kimax))./(1-Kimax)
        lmax = ((Kimax.*zi(indexmax)-Kimax))./(1-Kimax);
    end
    
    %Definitive initial values for "l" and "F".
    l = (lmin+lmax)/2;
    F = sum(zi.*(1-Ki)./(Ki+(1-Ki).*l));
    
    epsilon = 1000;%intial value to compare tolerance to find "l".
    
    while(epsilon>1e-6)% Verifies if tolerance for "l" is addressed.
        if F>0
            lmin = l;
        else
            lmax = l;
        end
        af = (1-Ki)./((1-Ki).*l+Ki);
        df = -sum(zi.*af.^2);
        l_old=l;
        l_NR = l-F/df;% Finding "l" with Newton-Raphson.
        
        if (l_NR<lmax) && (l_NR>lmin)% Verifying  if "l" with Newton-Raphson is between boundaries.
            l = l_NR;
        else
            l=(lmin+lmax)/2;% Finding "l" with Bisection.
        end
        epsilon = abs(l-l_old); % Difference between old and new "l".
        F = sum(zi.*(1-Ki)./(Ki+(1-Ki).*l)); %New F.
    end
      
    %Calculating xi and yi for each iteration until equilibrium.
    xi = zi./(l+Ki.*(1-l));
    yi = Ki.*zi./(l+Ki.*(1-l));

%   3.-Selecting root for vapor and liquid of Peng-Robinson equation.
 
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
    
    %Mixing rule for volume and pressure corrections.
    for i=1:length(ai);
        aL(i) = sum(xi(i).*xi.*((ai(i).*ai).^.5).*(1-kij(i,:)));
    end
    aL =sum(aL);
    bL = sum(xi.*bi);
    
    %Fugacity parameters for liquid
    for i=1:length(ai);
        deltaLi1(i) = sum(xi.*(ai.^.5).*(1-kij(i,:)));
        deltaLi(i) = (2.*(ai(i).^.5)./aL)*deltaLi1(i);
    end
    uu = (u^2-4*w)^.5;
    bi_bL = (Tc./Pc)./sum(xi.*Tc./Pc);
    
    %EOS coefficients for liquid
    AL = aL*P/((R^2)*(T^2));
    BL = bL*P/(R*T);
    a = 1;
    b = -(1+BL-u*BL);
    c = AL+w*BL^2-u*BL-u*BL^2;
    d = -AL*BL-w*BL^2-w*BL^3;
    EqL = [a b c d];
    solnPRL = roots(EqL);
    %Choosing right root for liquid.
    solnprlcool=solnPRL([1,3]); %Eliminating middle value.
    irootL = imag(solnprlcool); %Verifying if the possible roots are complex.
    [indeximL] = find(irootL==0); %Identying real roots.
    
    if length(indeximL) == 1 
        ZL = solnprlcool(indeximL); %If there is only one real root, that is the value of Z.
    else
        coefugL1 = bi_bL.*(solnprlcool(1)-1)-log(solnprlcool(1)-BL);
        coefugL11 = AL/(BL*(u^2-4*w)^.5).*(bi_bL-deltaLi);
        coefugL12 = 2*solnprlcool(1)+BL*(u+uu);
        coefugL13 = 2*solnprlcool(1)+BL*(u-uu);
        coefugL111 = coefugL1 + coefugL11.*log(coefugL12/coefugL13);
        coefugL111 = exp(coefugL111);
        
        coefugL3 = bi_bL.*(solnprlcool(2)-1)-log(solnprlcool(3)-BL);
        coefugL31 = AL/(BL*(u^2-4*w)^.5).*(bi_bL-deltaLi);
        coefugL32 = 2*solnprlcool(2)+BL*(u+uu);
        coefugL33 = 2*solnprlcool(2)+BL*(u-uu);
        coefugL333 = coefugL3 + coefugL31.*log(coefugL32/coefugL33);
        coefugL333 = exp(coefugL333);
        
        %If there two real roots, change in Gibbs-free energy is calculated.
        if solnprlcool(1)<solnprlcool(2)
            DeltaG = sum(coefugL111-coefugL333);
            if DeltaG>0
                ZL = solnprlcool(2);
                1;
            else
                ZL = solnprlcool(1);
                2;
            end
        else
            DeltaG = sum(coefugL333-coefugL111);
            if DeltaG>0
                ZL = solnprlcool(2);
                3;
            else
                ZL = solnprlcool(1);
                4;
            end
        end
        
    end
        
    %RIGHT FUGACITY LIQUID(Calculated with right Z)
    coefugLi = bi_bL.*(ZL-1)-log(ZL-BL);
    coefugLi1 = AL/(BL*(u^2-4*w)^.5).*(bi_bL-deltaLi);
    coefugLi2 = 2*ZL+BL*(u+uu);
    coefugLi3 = 2*ZL+BL*(u-uu);
    coefugLi = coefugLi + coefugLi1.*log(coefugLi2/coefugLi3);
    coefugLi = exp(coefugLi);
    fugLi = P.*coefugLi.*xi;
    
    % For Vapor
    
    %Mixing rule for volume and pressure corrections.
    for i=1:length(ai);
        aV(i) = sum(yi(i).*yi.*((ai(i).*ai).^.5).*(1-kij(i,:)));
    end
    aV =sum(aV);
    bV = sum(yi.*bi); 
     
    %Fugacity parameters for vapor
    for i=1:length(ai);
        deltaVi1(i) = sum(yi.*(ai.^.5).*(1-kij(i,:)));
        deltaVi(i) = (2.*(ai(i).^.5)./aV)*deltaVi1(i);
    end
    bi_bV = (Tc./Pc)./sum(yi.*Tc./Pc);
    
    %EOS coefficients for vapor
    AV = aV*P/((R^2)*(T^2));
    BV = bV*P/(R*T);
    a = 1;
    b = -(1+BV-u*BV);
    c = AV+w*BV^2-u*BV-u*BV^2;
    d = -AV*BV-w*BV^2-w*BV^3;
    EqV = [a b c d];
    solnPRV = roots(EqV);
    %Chhosing right root for vapor.
    solnprvcool=solnPRV([1,3]); %Eliminating middle value.
    irootV = imag(solnprvcool); %Verifyng if the posible roots are complex.
    [indeximV] = find(irootV==0); %Identiying real roots.
        
    if length(indeximV) == 1 %If there is only one real root, that is the value of Z.
        ZV = solnprvcool(indeximV);
    else
        coefugV1 = bi_bV.*(solnprvcool(1)-1)-log(solnprvcool(1)-BV);
        coefugV11 = AV/(BV*(u^2-4*w)^.5).*(bi_bV-deltaVi);
        coefugV12 = 2*solnprvcool(1)+BV*(u+uu);
        coefugV13 = 2*solnprvcool(1)+BV*(u-uu);
        coefugV111 = coefugV1 + coefugV11.*log(coefugV12/coefugV13);
        coefugV111 = exp(coefugV111);
        
        coefugV3 = bi_bV.*(solnprvcool(2)-1)-log(solnprvcool(2)-BV);
        coefugV31 = AV/(BV*(u^2-4*w)^.5).*(bi_bV-deltaVi);
        coefugV32 = 2*solnprlcool(2)+BV*(u+uu);
        coefugV33 = 2*solnprlcool(2)+BV*(u-uu);
        coefugV333 = coefugV3 + coefugV31.*log(coefugV32/coefugV33);
        coefugV333 = exp(coefugV333);
        
        %If there two real roots, change in Gibbs-free energy is calculated.
        if solnprvcool(1)<solnprvcool(2)
            DeltaG = sum(coefugV111-coefugV333);
            if DeltaG>0
                ZV = solnprvcool(2);
                1;
            else
                ZV = solnprvcool(1);
                2;
            end
        else
            DeltaG = sum(coefugV333-coefugV111);
            if DeltaG>0
                ZV = solnprvcool(2);
                3;
            elseif DeltaG<0
                ZV = solnprvcool(1);
                4;
            end
        end
        
    end
    
    %RIGHT FUGACITY VAPOR(Calculated with right Z)
    coefugVi = bi_bV.*(ZV-1)-log(ZV-BV);
    coefugVi1 = AV/(BV*(u^2-4*w)^.5).*(bi_bV-deltaVi);
    coefugVi2 = 2*ZV+BV*(u+uu);
    coefugVi3 = 2*ZV+BV*(u-uu);
    coefugVi = coefugVi + coefugVi1.*log(coefugVi2/coefugVi3);
    coefugVi = exp(coefugVi);
    fugVi = P.*coefugVi.*yi;
    
%   4.-Solving Ki to reach equilibrium.
    
    %Updating Ki and parameters to evaluate equilibrium.
    Ki = (fugLi./fugVi).*Ki;
    absdKiold = absdKi;
    dKi = (fugLi./fugVi)-1;
    absdKi = abs(dKi);
    
    %Checking for any change in the direction to reach equilibrium.
    if absdKi > absdKiold
        break
    end
    
    %Updating iteration counter to reach equilibrium.
    iter=iter+1;
end
%************************F L A S H  ENDS**********************************
