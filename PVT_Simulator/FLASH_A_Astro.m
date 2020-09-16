function [xi,yi,Ki,VL,VV,l,v]=FLASH_A_Astro(P,T,zi,Tci,Pci,wi) 
% INPUT: Total composition.
% OUTPUT: Liquid and Vapor composition and equilibrium constants.


% 1.- Calculating initial values and program parameters.
R = 8.31438;%Pa-m^3/mol-K
% P = 2000;%psia
% P = P*6892.8571;%Pa
% T = 162;%F
% T = (T-32)/1.8+273.15;%K
% Tci = [-116.63 305.65 652.1 87.9];%F
% Tci = ((Tci-32)./1.8)+273.15;%K
% Pci = [667.8 550.7 305.7 1071];%psia
% Pci = Pci.*6892.8571;%Pa
% wi = [.0104 .201 .49 .225];%Ascentric factor.
%--------------------------------------------------------------
% T = 162;%F
% T = (T-32)/1.8+273.15;%K
% P = 4.24*10^6;
% Tci = zeros(1,12);
% Pci = zeros(1,12);
% wi = zeros(1,12);
% zi = zeros(1,12);
% 
% Mole_fraction = [1.150 0.208 34.735 11.072 9.169 5.508 3.405 4.490 4.560 3.720 2.540 2.265];%percentage
% zi = (Mole_fraction./sum(Mole_fraction));%fraction
% 
% % 1.- N2
% Tci(1) = 126.2; 
% Pci(1) = 33.9*10^5; 
% wi(1) = 0.039; 
% % 2.- CO2
% Tci(2) = 304.1;
% Pci(2) = 73.8*10^5;
% wi(2) = 0.239;
% % 3.- C1
% Tci(3) = 190.4;
% Pci(3) = 46*10^5;
% wi(3) = 0.011;
% % 4.- C2
% Tci(4) = 305.4;
% Pci(4) = 48.8*10^5;
% wi(4) = 0.099;
% % 5.- C3
% Tci(5) = 369.8;
% Pci(5) = 42.5*10^5;
% wi(5) = 0.153;
% % 6.- C4
% Tci(6) = 425.2;
% Pci(6) = 38*10^5;
% wi(6) = 0.199;%nC4
% % 7.- C5
% Tci(7) = 469.7;
% Pci(7) = 33.7*10^5;
% wi(7) = 0.251;%nC5
% % 8.- C6
% Tci(8) = 507.5;
% Pci(8) = 30.1*10^5;
% wi(8) = 0.299; 
% % 9.- C7
% Tci(9) = 540.3;
% Pci(9) = 27.4*10^5;
% wi(9) = 0.349; 
% % 10.- C8
% Tci(10) = 568.8;
% Pci(10) = 24.9*10^5;
% wi(10) = 0.398;
% % 11.- C9
% Tci(11) = 594.6;
% Pci(11) = 22.9*10^5;
% wi(11) = 0.445;
% % 12.- C10
% Tci(12) = 617.7;
% Pci(12) = 21.2*10^5;
% wi(12) = 0.489;

%-----------------------------------------------------------------

Tr = T./Tci;%Reduced temperatue
E = ones(1,length(Tci))*1e-6;%Tolerance to reach Ki for equilibrium.
dKi = ones(1,length(Tci))*1000;%intial value to compare tolerance.
absdKi = abs(dKi);%intial value to compare tolerance.
Ki = (Pci./P).*exp(5.37.*(1+wi).*(1-(Tci./T)));%initial value of Ki with Wilson equation.
alfa = 0.5;%proportion to mix liquid and vapor and get zi.
%%
% ****************************F L A S H**********************************
iter=1; % To know how many iterations happen before to rerach equilibrium.
while(absdKi>=E);%Compares if the difference of previous and new Ki is inside of the tolerance.
    
%   2.- Solving Rachford-Rice equation to find “l”, “xi” and “yi”.
    
    %Setting the right and left boundaries of the search for "l".
    
    F0 = sum(zi.*(1-Ki)./Ki);
    F1 = 1-sum(zi.*Ki);
    if F0>0 && F1<0
        lmin = 0;
        lmax = 1;
    else
    [Kmin,Index_min] = min(Ki);
    zmin = zi(Index_min);
    [Kmax,Index_max] = max(Ki);
    zmax = zi(Index_max);
    lmin = ((zmin-Kmin)/(1-Kmin));
    lmax = (((Kmax*zmax)-Kmax)/(1-Kmax));
    end
    
    l = (lmin+lmax)/2;%Calculating "l" by bisection.
    F = sum(zi.*(1-Ki)./(Ki+(1-Ki).*l));

%     %Checking for tighter boundaries
%     [Kimin,indexmin] = min(Ki);
%     if Kimin>1%l<(zi(indexmin)-Kimin)./(1-Kimin)
%         lmin = (zi(indexmin)-Kimin)./(1-Kimin);
%     end
%     [Kimax,indexmax] = max(Ki);
%     if Kimax<1%l>((Kimax.*zi(indexmax)-Kimax))./(1-Kimax)
%         lmax = ((Kimax.*zi(indexmax)-Kimax))./(1-Kimax);
%     end
%     
%     %Definitive initial values for "l" and "F".
%     l = (lmin+lmax)/2;
%     F = sum(zi.*(1-Ki)./(Ki+(1-Ki).*l));
%     
%     % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    syms x;
    eqn = sum((zi.*(1-Ki))./(Ki+(1-Ki).*x)) == 0;
    solx = vpasolve(eqn,x,[lmin,lmax]);
    lx = double(solx);
 
    
    if length(lx)==1
        l=lx;
    else
        
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
    end

     

    %Calculating xi and yi for each iteration until equilibrium.
    v = 1-l;
    xi = zi./(l+Ki.*(1-l));
    yi = Ki.*zi./(l+Ki.*(1-l));

    
    %Normalizing Xi and yi
    xi = xi./sum(xi);
    yi = yi./sum(yi);

%   3.-Selecting root for vapor and liquid of Peng-Robinson equation.
 
    u = 2;
    w = -1;
    kij = zeros(length(wi),length(wi));
    %Values from the Table 1 of SPE - 116823 paper.
%     kij(2,1) = 0.027;
%     kij(3,1) = 0.042;
%     kij(4,1) = 0.1;
%     kij(3,2) = 0.008;
%     kij(4,2) = 0.1257;
%     kij(4,3) = 0.0942;
    bi = .0778*R.*Tci./Pci;
    fwi = .37464+1.54226.*wi-.26992.*(wi.^2);
    ai = .45724*(R^2).*(Tci.^2)./Pci;
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
    bi_bL = (Tci./Pci)./sum(xi.*Tci./Pci);
    
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
        
        coefugL3 = bi_bL.*(solnprlcool(2)-1)-log(solnprlcool(2)-BL);
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
    
    %Volume
    VL = ZL*R*T/P;
    
    % For Vapor
    
    %Mixing rule for volume and pressure corrections.
    for i=1:length(ai)
        aV(i) = sum(yi(i).*yi.*((ai(i).*ai).^.5).*(1-kij(i,:)));
    end
    aV =sum(aV);
    bV = sum(yi.*bi); 
     
    %Fugacity parameters for vapor
    for i=1:length(ai);
        deltaVi1(i) = sum(yi.*(ai.^.5).*(1-kij(i,:)));
        deltaVi(i) = (2.*(ai(i).^.5)./aV)*deltaVi1(i);
    end
    bi_bV = (Tci./Pci)./sum(yi.*Tci./Pci);
    
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
    %Volume
    VV = ZV*R*T/P;
    
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
    if iter>=30
        break
    end
    %Updating iteration counter to reach equilibrium.
    iter=iter+1;
 
end
    if l>=1
        l=1;
    elseif l<0
        l=0;
    end
    v = 1-l;
    
    xi = zi./(l+Ki.*(1-l));
    yi = Ki.*zi./(l+Ki.*(1-l));

    
    %Normalizing Xi and yi
    xi = xi./sum(xi);
    yi = yi./sum(yi);
%************************F L A S H  ENDS**********************************
