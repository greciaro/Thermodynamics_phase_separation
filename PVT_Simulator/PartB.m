
P=[55.2615014,51.8141214,48.3667414,44.9193614,41.4719814,39.989608,38.0246014,34.5772214,31.1298414,27.6824614,24.2350814,...
20.7877014,17.3403214,15.9613694,13.8929414,10.4455614,6.9981814,3.5508014,2.1718494,0.7928974,0.1034214].*10^6;
visexp = [0.653,0.635,0.617,0.601,0.585,0.578,0.57,0.556,0.544,0.53,0.514,0.497,0.472,0.46,0.512,0.607,0.698,0.785,0.865];
vis = zeros(length(Pexp));
for i=1:length(Pexp)
vis10(i)=Viscmix(Tci10,Pci10,VL10(i),xi10,Mwi10,wi10,T,Pexp(i));
vis20(i)=Viscmix(Tci20,Pci20,VL20(i),xi20,Mwi20,wi20,T,Pexp(i));
vis40(i)=Viscmix(Tci40,Pci40,VL40(i),xi40,Mwi40,wi40,T,Pexp(i));
end


close all
plot(Pexp./(10^6),vis40,'b','LineWidth',2);
hold on
plot(Pexp./(10^6),vis20,'r','LineWidth',2);
plot(Pexp./(10^6),vis10,'g','LineWidth',2);
plot(Pexp./(10^6),visexp,'ko','LineWidth',2);
hold off


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
ygas(4,:) = [0.0036,0.0198,0.6010,0,0,0.3033,0.0678,0.0042,0,0,0,0];

iii=1;
jjj=1;
T = 72.2+273.15;%K
tic
for jjj = 1:4;
    figure('Color','w');set(gcf,'Position',[553   -49   590   354])
    for P = [6 9 10 11 12 13 13.05].*10^6;%Pa
%         P = 12*10^6;
        %         for cellss = [152 252 302 402]
        cellss = 282;
        contacts = (cellss-2)/2;
        alfa = 0.5;
        cells = zeros(cellss,length(xoil),contacts);
        ki = zeros(cellss,length(xoil),contacts);
        for ii=1:contacts
            cells(1,:,ii) = ygas(jjj,:);
            cells(2*ii,:,ii) = xoil;
        end
        for ii=2:contacts
            for i=2:2:ii*2-1
                z=cells(i,:,ii-1)+alfa.*(cells(i-1,:,ii-1)-cells(i,:,ii-1));
                [x,y,k]=FLASH_mmp(P,T,z,Tci,Pci,wi);
                ki(i,:,ii)=k;
                cells(i,:,ii)=x;
                cells(i+1,:,ii)=y;
                ii,jjj
            end
        end
        for j=2:2:cellss-2
            %         Ln(j-1) = sqrt((cells(j,:,contacts)-cells(j-1,:,contacts)).^2);
            Ln(j-1) = sqrt(sum((cells(j,:,contacts)-cells(j-1,:,contacts)).^2));
        end
        Lnp(1:length(Ln),jjj) = Ln;
        
        iii=iii+1;
    end
    hold off
end

            plot(min(Ln(1,:)),'.','Markersize',8);
            % ylim([0 .8]);
            hold on

            
            %     [aa,bb]=findpeaks(-Ln(3:2:end-4));
            %     [aam,bbm]=min(Ln(3:2:end-4));
            
            

            
            %
            
            %     plot(P,Ln(23),'.b','Markersize',10);
            %     hold on
            %     plot(P,Ln(123),'.','Color',[0 1 0],'Markersize',10);
            %     plot(P,Ln(223),'.r','Markersize',10);
            
            
            %     plot(P,aa(end),'.r','Markersize',20);
            %     plot(P,aam,'.g','Markersize',20);
            %     plot(Ln(3:2:end-4),'.','Markersize',20);
            %     plot(Ln(:),'.','Markersize',20);
            %     hold on
            %  clear ki cells z x y k
            
            
            % matr(iii)=[Ln(123)];
            
            
            % matr2(iii)=[Ln2(123)];
            % matr3(iii)=[Ln3(123)];
            % matr4(iii)=[Ln4(123)];
           
%         end
        title(['P= ',num2str(P/10^6),'Mpa  Gas inj = ',num2str(jjj)]);
        xlabel('Cell Number');
        ylabel('Tie-line length');
%         legend('150 Cells','250 Cells','300 Cells',' 400 Cells');
        legend('400 Cells');
        grid on
        hold off
%         end
    
        
    % pp=polyfit([500 750 1250 1400 1700 2000 2050 2060 2070 2080 2074 2303],[matr,.2, 0],6);
    % plot( linspace(1700,2303,1000),polyval(pp,linspace(1700,2303,1000)),'-','Color',[0 .6 0],'LineWidth',2)
    % xlim([500 2500]);
    % ylim([0 .65]);
    % xlabel('Pressure (psi)');
    % ylabel('Tie-line length');
    % text(1550,.21,'CO_2')
    % text(1550,.28,'C_1_0')
    % text(1550,.11,'C_1')
    % text(1550,.03,'C_4')
    
    
    
    % pp1=polyfit([500  550 600 750 800 850 900 1000 1250 1400 1700 2000 2030 2040 2050 2055 2060 2065 2070 2075 2080 2303],[matr1,0],6);
    % plot( linspace(1700,2303,1000),polyval(pp1,linspace(1700,2303,1000)),'-','Color',[.8 0 0],'LineWidth',.9)
    % pp2=polyfit([500  550 600 750 800 850 900 1000 1250 1400 1700 2000 2030 2040 2050 2055 2060 2065 2070 2075 2080 2303],[matr2,0],6);
    % plot( linspace(1700,2303,1000),polyval(pp2,linspace(1700,2303,1000)),'-','Color',[0 .8 0],'LineWidth',.9)
    % pp3=polyfit([500  550 600 750 800 850 900 1000 1250 1400 1700 2000 2030 2040 2050 2055 2060 2065 2070 2075 2080 2303],[matr3,0],6);
    % plot( linspace(1700,2303,1000),polyval(pp3,linspace(1700,2303,1000)),'-','Color',[0 0 .8],'LineWidth',.9)
    % pp4=polyfit([500  550 600 750 800 850 900 1000 1250 1400 1700 2000 2030 2040 2050 2055 2060 2065 2070 2075 2080 2303],[matr4,0],6);
    % plot( linspace(1700,2303,1000),polyval(pp4,linspace(1700,2303,1000)),'-','Color',[.8 0 .8],'LineWidth',.9)
    
    
    
% end
% toc



