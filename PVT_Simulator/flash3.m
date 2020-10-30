   3.-Selecting root for vapor and liquid of Peng-Robinson equation.
 
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
    
%