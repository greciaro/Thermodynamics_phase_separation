function [ A,B,C,D ] = ComputeABCD(Cmax, Cplus, Zplus, Mplus, ...
                                   densityCplus, lastFracDensity)
% This Function Computes the values for A, B, C, D in Pedersen et al.

% Cmax  : maximum carbon number to expand to (e.g. 200)
% Cplus : Carbon Number from Which the expansion starts for example C6+ or
%         C11 + etc..
% Zplus : Mole fraction for the Cplus fraction
% Mplus : Molecular weight for the plus fraction g/mole
% densityCplus : Density of the plus fraction g/cm3
% lastFracDensity : Density of the last measured fraction
%
% by Sergey Klevtsov, 2014 (adapted from Waqas Ali, 2012)

if (Cmax < Cplus),
    error('Max carbon number is less then plus fraction carbon number!');
end;
Crange = Cplus:Cmax;
M = 14 * Crange - 4;
J = zeros(2, 2);
eps = 1e-4;

% initial guess for A and B
v = [-1.0; -0.2];

while true
    z = exp(v(1) + v(2) * Crange);
    % Computing the residual vector
    r = [ Zplus - sum(z) 
          Zplus*Mplus - z * M' ];
    % Convergence criteria
    if all(abs(r) < eps),
        break;
    end;
    % Computing the Jacobian Matrix.
    J(1,1) = - sum(z);
    J(1,2) = - z * Crange';
    J(2,1) = - z * M';
    J(2,2) = - (z .* M) * Crange';
    % Solving the system of equation to find new iteration value
    v = v - J\r;
end;

A = v(1);
B = v(2);
z = exp(A + B * Crange);
sum_zM = z * M';

% initial guess for C and D
v = [0.70; 1.0];

while true
    rho = v(1) + v(2) * log(Crange);
    denominator = z * (M ./ rho)';
    % Computing the residual vector
    r = [ densityCplus - sum_zM / denominator
          lastFracDensity - (v(1) + v(2) * log(Cplus - 1)) ];
    % Convergence criteria
    if all(abs(r) < eps),
        break;
    end;
    % Computing the Jacobian Matrix.
    J(1,1) = - z * (M ./ rho.^2)' * sum_zM / denominator^2;
    J(1,2) = - z * (M .* log(Crange) ./ rho.^2)' * sum_zM / denominator^2;
    J(2,1) = - 1;
    J(2,2) = - log(Cplus - 1);
    % Solving the system of equation to find new iteration value
    v = v - J\r;
end;

C = v(1,1);
D = v(2,1);

end

