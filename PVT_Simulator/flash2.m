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

%
