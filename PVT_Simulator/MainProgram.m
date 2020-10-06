%Main program
%clear
load('Data.mat');
Compo40;
Compo20;
Compo10;
%%
format short g;
char('40 Components EOS parameters');
[(1:length(Tci40))',Tci40',Pci40',zi40',wi40'];
format short g;
char('20 Components EOS parameters');
[(1:length(Tci20))',Tci20',Pci20',zi20',wi20'];
format short g;
char('10 Components EOS parameters');
[(1:length(Tci10))',Tci10',Pci10',zi10',wi10'];
%% Plotting original curves
% CCE  Pressure vs Relative Volume picture
figure('Color','w');set(gcf,'Position',[553   -49   590   354]);
plot(Pexp./(10^6),RelV40,'b','LineWidth',1.5); hold on
plot(Pexp./(10^6),RelV20,'r','LineWidth',1.5);
plot(Pexp./(10^6),RelV10,'g','LineWidth',1.5); 
plot(Pexp./(10^6),RelVexp,'ko','LineWidth',1.5);
title('Pressure VS Relative Volume');xlabel('Pressure (Mpa)');ylabel('Relative Volume');
legend('Sim 40 components','Sim 20 components','Sim 10 components','Experimental Data');
grid on;hold off
% Liquid Fraction
figure('Color','w');set(gcf,'Position',[553   -49   590   354]);
plot(Pexp(14:19)./(10^6),l40(14:19),'b','LineWidth',6); hold on
plot(Pexp(14:19)./(10^6),l20(14:19),'r','LineWidth',3); 
plot(Pexp(14:19)./(10^6),l10(14:19),'g','LineWidth',1.5);
xlabel('Pressure (Mpa)');ylabel('Liquid Fraction');
legend('40 components','20 components','10 components');
grid on;hold off
% Isothermal Oil Compresibility in 1/psia
figure('Color','w');set(gcf,'Position',[553   -49   590   354]);
plot(Pexp(1:14)./(10^6),Compres40,'b','LineWidth',6); hold on;
plot(Pexp(1:14)./(10^6),Compres20,'r','LineWidth',3); 
plot(Pexp(1:14)./(10^6),Compres10,'g','LineWidth',1.5); 
plot(Pexp(1:14)./(10^6),Compresexp,'ko','LineWidth',2);
xlabel('Pressure (Mpa)'); ylabel(' Isothermal Oil Compresibility (1/Psi)');
legend('40 components','20 components','10 components','Experimental Data');
grid on;hold off
%% Tuning
Compo10Tuning;
Compo20Tuning;
Compo40Tuning;
load('Comp10.mat');
load('Comp20.mat');
load('Comp40.mat');
%%
% CCE  Pressure vs Relative Volume picture
figure('Color','w');set(gcf,'Position',[553   -49   590   354]);
plot(Pexp./(10^6),RelV40,'b','LineWidth',1.5); hold on
plot(Pexp./(10^6),RelV20,'r','LineWidth',1.5);
plot(Pexp./(10^6),RelV10,'Color',[0 .7 0],'LineWidth',1.5);
plot(Pexp./(10^6),RelV40t,'.b','MarkerSize',15,'LineWidth',1.5);
plot(Pexp./(10^6),RelV20t,'.','Color',[1 0 0],'MarkerSize',15,'LineWidth',1.5);
plot(Pexp./(10^6),RelV10t,'.','Color',[0 .7 0],'MarkerSize',15,'LineWidth',1.5);
plot(Pexp./(10^6),RelVexp,'ko','LineWidth',1.2);
title('Pressure VS Relative Volume');xlabel('Pressure (Mpa)');ylabel('Relative Volume');
legend('40 Comp', '20 comp', '10 Comp','Tunned 40 Comp','Tunned 20 comp',...
    'Tunned 10 Comp','Experimental Data');
grid on;hold off
%% Liquid Fraction
figure('Color','w');set(gcf,'Position',[553   -49   590   354]);
plot(Pexp(14:19)./(10^6),l40(14:19),'b','LineWidth',1.5); hold on
plot(Pexp(14:19)./(10^6),l20(14:19),'r','LineWidth',1.5);
plot(Pexp(14:19)./(10^6),l10(14:19),'Color',[0 .7 0],'LineWidth',1.5);
plot(Pexp(14:19)./(10^6),l40t(14:19),'.','Color',[0 0 1],'MarkerSize',15,'LineWidth',1.5);
plot(Pexp(14:19)./(10^6),l20t(14:19),'.','Color',[1 0 0],'MarkerSize',15,'LineWidth',1.5);
plot(Pexp(14:19)./(10^6),l10t(14:19),'.','Color',[0 .7 0],'MarkerSize',15,'LineWidth',1.5);
xlabel('Pressure (Mpa)'); ylabel('Liquid Fraction');
legend('40 Comp', '20 comp','10 Comp', 'Tunned 40 Comp','Tunned 20 comp',...
    'Tunned 10 Comp');
grid on;% hold off
%% Isothermal Oil Compresibility in 1/psia
figure('Color','w');set(gcf,'Position',[553   -49   590   354]);
plot(Pexp(1:14)./(10^6),(Compres40),'b','LineWidth',1.5); hold on;
plot(Pexp(1:14)./(10^6),(Compres20),'r','LineWidth',1.5);
plot(Pexp(1:14)./(10^6),(Compres10),'Color',[0 .7 0],'LineWidth',1.5); 
plot(Pexp(1:14)./(10^6),(Compres40t),'.','Color',[0 0 1],'MarkerSize',15,'LineWidth',1.5); 
plot(Pexp(1:14)./(10^6),(Compres20t),'.','Color',[1 0 0],'MarkerSize',15,'LineWidth',1.5);
plot(Pexp(1:14)./(10^6),(Compres10t),'.','Color',[0 .7 0],'MarkerSize',15,'LineWidth',1.5); 
plot(Pexp(1:14)./(10^6),Compresexp,'ko','LineWidth',1.5);
xlabel('Pressure (MPa)'); ylabel(' Isothermal Oil Compresibility (1/Psi)');
legend('40 Comp', '20 comp','10 Comp', 'Tunned 40 Comp','Tunned 20 comp',...
    'Tunned 10 Comp','Experimental');axis([0 60 0 20]);
grid on; hold off;
%%
clear
load('Comp10.mat');
% load('Comp10t.mat');
PTCMEdiagram;
figure('Color','w');set(gcf,'Position',[553   -49   590   354]);
surf(real(lPT)); view(2)
xlabel('Temperature (°C)'); ylabel('Pressure (MPa)'); title('P-T diagram');
set(gca,'XTickLabel',round(linspace(21,490,8).*100)./100);
set(gca,'YTickLabel',round(linspace(3,23.5,9).*10)./10);
axis tight
figure('Color','w');set(gcf,'Position',[553   -49   590   354]);
plot(pointTPT,pointPPT,'bo','LineWidth',1.5); hold on
plot(pointTPTlow,pointPPTlow,'ro','LineWidth',1.5);
plot(linspace(-45,490,100)+273.15,rrtest,'k','LineWidth',1.5);
xlabel('Temperature (K)'); ylabel('Pressure (Pa)'); title('P-T diagram of 10 Components');
grid on
clear
load('Comp10t.mat');
PTCMEdiagramt;
plot(pointTPT,pointPPT,'*'); hold on
plot(pointTPTlow,pointPPTlow,'*');

plot(linspace(-45,550,100)+273.15,rrtest,'r','LineWidth',3);
hold off; axis square
%%
clear
load('Data.mat');
load('Comp10.mat');
load('Comp20.mat');
load('Comp40.mat');
visexp = [0.653,0.635,0.617,0.601,0.585,0.578,0.57,0.556,0.544,0.53,0.514,0.497,0.472,0.46,0.512,0.607,0.698,0.785,0.865];
vis = zeros(length(Pexp));
for i=1:length(Pexp)
vis10(i)=Viscmix(Tci10,Pci10,VL10(i),xi10,Mwi10,wi10,T,Pexp(i));
vis20(i)=Viscmix(Tci20,Pci20,VL20(i),xi20,Mwi20,wi20,T,Pexp(i));
vis40(i)=Viscmix(Tci40,Pci40,VL40(i),xi40,Mwi40,wi40,T,Pexp(i));
end
figure('Color','w');set(gcf,'Position',[553   -49   590   354])
plot(Pexp./(10^6),vis40,'b','LineWidth',2); hold on
plot(Pexp./(10^6),vis20,'r','LineWidth',2);
plot(Pexp./(10^6),vis10,'g','LineWidth',2);
plot(Pexp./(10^6),visexp,'ko','LineWidth',2);
ylabel('Viscosity (cp)'); xlabel('Pressure (MPa)'); title('Oil-phase viscosity');
hold off;

%% 
%mmp
clear all
%Composition
% 1    2   3   4  5   6     7      8       9        10      11    12
% CO2, N2, C1,C2,C3,C2-C3,C4-C6, C7-C9, C10-C13, C14-C18, C19-C26,C27+
xoil = [0.002,0.0106,0,0,0.3208,0.1895,0.1215,0.1428,0.0815,0.0595,0.0445,0.0273];
Mwi = [44,28,16,30.07,44.09,37.778,70.8116,102.5795,157.8634,218.4134,305.5046,510.7578]; %g/mole
Densi = [1.84,1.165,0.668,0,0,1.5484,349.8809,721.0376,782.859,837.701,892.6286,973.4299]; %kg/m^3
Tci = [304.2,126,190.6,305.4,369.8,335.0331,456.6255,572.6343,629.2767,689.7051,773.1023,941.3879];%K
Pci = [74.8792,34.4505,46.6095,48.8*1.01325,42.5*1.01325,46.5093,35.867,29.3561,22.4496,18.4228,15.4334,12.6298];%atm
Pci = Pci.*101325;%To change from atm to Pa
wi = [0.23,0.04,0.01,0.099,0.153,0.123,0.2149,0.3096,0.5626,0.8472,1.1562,1.3681];
ygas(1,:) = [1,0,0,0,0,0,0,0,0,0,0,0];
ygas(2,:) = [0.5,0,0.5,0,0,0,0,0,0,0,0,0];
ygas(3,:) = [0.22,0,0.38,0.20,0.20,0,0,0,0,0,0,0];
ygas(4,:) = [0.0036,0.0198,0.6010,0,0,0.3033,0.0678,0.0042,0,0,0,0]./sum([0.0036,0.0199,0.6010,0,0,0.3033,0.0678,0.0042,0,0,0,0]);
jjj=1;
T = 72.2+273.15;%K
        cellss = 282;
        contacts = (cellss-2)/2;
        alfa = 0.5;
tic
 for jjj = 1:4;
     iii=1;
     for P =[3 6 8 9 10.3 10.5 10.8].*10^6.;%Pa
         Pp=P;
         if jjj==2
             Pp=Pp*2;
         elseif jjj==3
             Pp=Pp*1.20;
         elseif jjj==4
             Pp=Pp*2.5;
         end
         cells = zeros(cellss,length(xoil),contacts);
         ki = zeros(cellss,length(xoil),contacts);
         for ii=1:contacts
             cells(1,:,ii) = ygas(jjj,:);
             cells(2*ii,:,ii) = xoil;
         end
         for ii=2:contacts
             for i=2:2:ii*2-1
                 z=cells(i,:,ii-1)+alfa.*(cells(i-1,:,ii-1)-cells(i,:,ii-1));
                 z = z./sum(z);
                 %                 [x,y,k]=FLASH_mmp(P,T,z,Tci,Pci,wi);
                 [x,y,k,~]=FLASH_A_Astro(Pp,T,z,Tci,Pci,wi);
                 ki(i,:,ii)=k;
                 cells(i,:,ii)=x;
                 cells(i+1,:,ii)=y;
                 %                 ii;jjj;
             end
         end
         for j=2:2:cellss-2
             %         Ln(j-1) = sqrt((cells(j,:,contacts)-cells(j-1,:,contacts)).^2);
             Ln(j-1) =sum(sqrt((cells(j,:,contacts)-cells(j-1,:,contacts)).^2));
         end
         iii;
         jjj;
         Lnp(iii,jjj) = min(Ln(Ln~=0));
         PPP(iii,jjj)=Pp;
         
         iii=iii+1;
         
     end
 end
 Lnp1=Lnp;
 PPP1=PPP;
 toc
 %%
 figure('Color','w');set(gcf,'Position',[553   -49   590   354])
 Lnp(8,1:4)=[0,0,0,0];
 PPP(8,1:4)=[1.24,2.57,1.322,3.07].*10^7;
 for ii=1:4
     if ii==1 | ii==4
         r=polyfit(PPP(4:7,ii)',Lnp(4:7,ii)',4);
     elseif ii==2
         r=polyfit(PPP(4:7,ii)',Lnp(4:7,ii)',4);
     else
         r=polyfit(PPP(3:6,ii)',Lnp(3:6,ii)',4);
     end
     
     rr=polyval(r,[PPP(5,ii)/10^6:.1:37].*10^6);
     plot([PPP(5,ii)/10^6:.1:37].*10^6,rr,'LineWidth',2);
     hold on
     plot(PPP1(:,ii),Lnp1(:,ii),'.','MarkerSize',20);
 end
 axis([5e6 max(max(PPP))+5e6 0 .8]);
 save('MMP.mat')
 
 