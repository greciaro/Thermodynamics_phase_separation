
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
kkk=1;
for jjj = 1:4;
%     figure('Color','w');set(gcf,'Position',[553   -49   590   354])
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
%                 ii,jjj
            end
        end
        for j=2:2:cellss-2
            %         Ln(j-1) = sqrt((cells(j,:,contacts)-cells(j-1,:,contacts)).^2);
            Ln(j-1) = sqrt(sum((cells(j,:,contacts)-cells(j-1,:,contacts)).^2));
        end
        Lnp(kkk,iii) = min(Ln);
        kkk=kkk+1;
        iii=iii+1;
    end
end
