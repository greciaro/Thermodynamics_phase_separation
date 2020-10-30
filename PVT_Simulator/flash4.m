   4.-Solving Ki to reach equilibrium.
    
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