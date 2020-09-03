function [yf] = F( l,Ki,zi )

yf = sum(zi.*(1-Ki)./(Ki+(1-Ki).*l));

end

