
% Cplus = 11;
% Zplus = 0.17178;
% Mplus = 304.993;
% densityCplus = 0.875;
% lastFracDensity = 0.7684;

% This Function Computes the values for A and B
% 
% Cplus : Carbon Number from Which the expansion starts for example C6+ or
%         C11 + etc..
% Zplus : Mole fraction for the Cplus fraction
% Mplus : Molecular weight for the plus fraction g/mole
% densityCplus : Density of the plus fraction g/cm3
% lastFracDensity : Density of the last measured fraction
%
% by Waqas Ali, 2012


%******************************Start**************************************

mass_fracwtplus=Zplus*Mplus;

% Mplus=292.7;

Zplus=mass_fracwtplus/Mplus;

ZplusMplus = Zplus*Mplus;

% Finding the density of the C11+ fraction
% *****************************End of section 2****************************

%******************************Start***************************************
% Section 3
% Initializing the constants A & B to be later used for iteraiton.
A=-1.0;
B=-0.2;
% ****************************End of section 3*****************************

% *****************************Start*************************************
% Section 4
% Formation of the Jacobi matrix and evaluating the constants A and B
% Initializing the the terms of Jacobian Matrix.
Jacobian=zeros(2,2);
% Vector "v" stores the value of A and B at a current iteration step.
v=zeros(2,1);
v(1,1)=A;
v(2,1)=B;
n_eq=0.0;

while (n_eq~=2.0)

% Computing the first term of the Jacobian Matrix.
df1_dA=0.0;

for Ci=Cplus:200% The carbon number starts from C11 and ends at 200.
    df1_dA=df1_dA-exp(v(1,1)+v(2,1)*Ci);
end

Jacobian(1,1)=df1_dA;
% Computing the second term of the Jacobian Matrix.
df1_dB=0.0;


for Ci=Cplus:200
    df1_dB=df1_dB-(exp(v(1,1)+v(2,1)*Ci))*Ci;
end

Jacobian(1,2)=df1_dB;
% Computing the third term of Jacobian Matrix.
df2_dA=0.0;

for Ci=Cplus:200
    df2_dA=df2_dA-exp(v(1,1)+v(2,1)*Ci)*(14*Ci-4);
end

Jacobian(2,1)=df2_dA;
% Computing the fourth term of the Jacobian Matrix.

df2_dB=0.0;
for Ci=Cplus:200
    df2_dB=df2_dB-exp(v(1,1)+v(2,1)*Ci)*(14*Ci-4)*Ci;
end

Jacobian(2,2)=df2_dB;
% computing the value of the of f1
f1_AB=0.0;

for Ci=Cplus:200;
    f1_AB=f1_AB-exp(v(1,1)+v(2,1)*Ci);
end

f1_AB=Zplus+f1_AB;
% computing the value of the of f2
f2_AB=0.0;

for Ci=Cplus:200;
    f2_AB=f2_AB-exp(v(1,1)+v(2,1)*Ci)*(14*Ci-4);
end

f2_AB=ZplusMplus+f2_AB;
% fvec is a vector that stores the function value at each step.
fvec=zeros(2,1);
fvec(1,1)=f1_AB;
fvec(2,1)=f2_AB;
err=zeros(2,1);

% Solving the system of equation to find new iteration value
err=Jacobian\(-fvec);

vk_plus1=v+err;

% evaluating the funciton at next iteration value for checking convergence 
% criteria

f1_AB=0.0;
for Ci=Cplus:200;
    f1_AB=f1_AB-exp(vk_plus1(1,1)+vk_plus1(2,1)*Ci);
end
f1_AB=Zplus+f1_AB;

% computing the value of f2 at k+1 iteration
f2_AB=0.0;
for Ci=Cplus:200;
    f2_AB=f2_AB-exp(vk_plus1(1,1)+vk_plus1(2,1)*Ci)*(14*Ci-4);
end
f2_AB=ZplusMplus+f2_AB;

fvec_plus1=zeros(2,1);
fvec_plus1(1,1)=f1_AB;
fvec_plus1(2,1)=f2_AB;

% Finding whether both functions yeilds function value of zero.
n_eq=0.0;
for i=1:2
    if (abs(fvec_plus1(i,1))<0.0001)
        n_eq=n_eq+1;
    else
    end
end
           
v=vk_plus1;

end

%***************************End of section 4**********************************

%******************************Start**************************************
% Section 5 computing the composition of the expanded fractions
A=v(1,1);
B=v(2,1);
comp_exp=zeros(200-Cplus+1,1);

for Ci=Cplus:200
   comp_exp(Ci-(Cplus-1),1)=exp(A+B*Ci);
end

% ********************************End of section 5*************************

% %*********************************Start************************************
% % Section 6 Plotting the expanded composition 
% subplot(1,2,1),plot(11:200,comp_exp,'go-','MarkerSize',3,'LineWidth',2)...
% ,xlabel('Carbon Number'),ylabel('Molar composition'),title('Compositional Distribution')
% legend('Mole Fraction');
% grid on
% % ********************************End of section 6*************************


%*********************************Start************************************
% Section 7 
% Now we move on to calculate the coefficient 'C' and 'D'
% Firstly computing the numerator of equation 7.29 which is the sum of the
% product of zi and Mi

M=zeros(200-Cplus+1,1);
sum_ziMi200=0.0;
for i=Cplus:200
    M(i-(Cplus-1),1)=14*i-4;
    sum_ziMi200=sum_ziMi200+comp_exp(i-(Cplus-1),1)*M(i-(Cplus-1),1);
end


%***************************End of section 4**********************************



%******************************Start**************************************
% Section 5 computing the composition of the expanded fractions
A=v(1,1);
B=v(2,1);
comp_exp=zeros(200-Cplus+1,1);

for Ci=Cplus:200
   comp_exp(Ci-(Cplus-1),1)=exp(A+B*Ci);
end
% ********************************End of section 5*************************

% %*********************************Start************************************
% % Section 6 Plotting the expanded composition 
% subplot(1,2,1),plot(11:200,comp_exp,'go-','MarkerSize',3,'LineWidth',2)...
% ,xlabel('Carbon Number'),ylabel('Molar composition'),title('Compositional Distribution')
% legend('Mole Fraction');
% grid on
% % ********************************End of section 6*************************


%*********************************Start************************************
% Section 7 
% Now we move on to calculate the coefficient 'C' and 'D'
% Firstly computing the numerator of equation 7.29 which is the sum of the
% product of zi and Mi
M=zeros(200-Cplus+1,1);
sum_ziMi200=0.0;
for i=Cplus:200
    M(i-(Cplus-1),1)=14*i-4;
    sum_ziMi200=sum_ziMi200+comp_exp(i-(Cplus-1),1)*M(i-(Cplus-1),1);
end
% assuming the value of C and D.
C=0.70;
D=1.0;
%********************************End of Section7***************************

%*********************************Start************************************
% Section 8
% Initializing the terms of Jacobian Matrix.
Jacobian=zeros(2,2);
% Vector "v" stores the value of C and D at a current iteration step.
v=zeros(2,1);
v(1,1)=C;
v(2,1)=D;
n_eq=0.0;
% Computing the first term of Jacobian Matrix.
while (n_eq~=2.0)
denominator=0.0;
for Ci=Cplus:200
%Denominator calculates sum of zi*Mi/(C+D*lnCi)
denominator=denominator+comp_exp(Ci-(Cplus-1),1)*M(Ci-(Cplus-1),1)/(v(1,1)+v(2,1)*log(Ci));
end
% Denominator square is utilized in calculating the derivative.
denominatorsq=denominator^2;

numerator_1=0.0;
for Ci=Cplus:200
%numerator_1 is utilized in derivative wrt to C.
numerator_1=numerator_1+comp_exp(Ci-(Cplus-1),1)*M(Ci-(Cplus-1),1)/(v(1,1)+v(2,1)*log(Ci))^2;
end
Jacobian(1,1)=-numerator_1*sum_ziMi200/denominatorsq;


% Evaluating the second term of Jacobian Matrix.
numerator_2=0.0;
for Ci=Cplus:200
%numerator_2 is utilized in derivative wrt to D.    
numerator_2=numerator_2+comp_exp(Ci-(Cplus-1),1)*M(Ci-(Cplus-1),1)*log(Ci)/(v(1,1)+v(2,1)*log(Ci))^2;
end
Jacobian(1,2)=-numerator_2*sum_ziMi200/denominatorsq;

% Evaluating third term of jacobian matrix.
Jacobian(2,1)=-1;

% Evaluating the fourth term of jacobian matrix
Jacobian(2,2)=-log(10);

% Evaluating the first term of function vector
f1_CD=densityCplus-sum_ziMi200/denominator;
% Evaluating the second term of function vector
f2_CD=lastFracDensity-(v(1,1)+v(2,1)*log(10));
% Formation of function vector.
fvec=zeros(2,1);
fvec(1,1)=f1_CD;
fvec(2,1)=f2_CD;

% Solving the matrix vector equation
err=Jacobian\(-1*fvec);
% finding the new iteration step
vk_plus1=v+err;
% Evaluating the function at next iteration step to check convergence.
denominator=0.0;
for Ci=Cplus:200
    denominator=denominator+comp_exp(Ci-(Cplus-1),1)*M(Ci-(Cplus-1),1)/(vk_plus1(1,1)+vk_plus1(2,1)*log(Ci));
end

f1_CD=densityCplus-sum_ziMi200/denominator;
f2_CD=lastFracDensity-(vk_plus1(1,1)+vk_plus1(2,1)*log(10));


fvec_plus1=zeros(2,1);
fvec_plus1(1,1)=f1_CD;
fvec_plus1(2,1)=f2_CD;
% Checking whether the convergence criteria is met for both functions
n_eq=0.0;
for i=1:2
    if (abs(fvec_plus1(i,1))<0.0001)
        n_eq=n_eq+1;
    else
    end
end
           
v=vk_plus1;
end
%********************************End of section 8**************************


%*************************************Start********************************
% Section 9
C=v(1,1);
D=v(2,1);
% Computing density of all expanded fractions.


end

