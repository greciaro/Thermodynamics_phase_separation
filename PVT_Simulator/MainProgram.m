
% TO CREATE FIGURE 2
P = 2000;%psia
figure('Color','w');set(gcf,'Position',[553   -49   590   354])
for cellss = [52 102 152 252]
contacts = (cellss-2)/2;
%Compositions of C1, C4, C10 and CO2
xoil = [0.2 0.15 0.65 0];% Oil
ygas = [0.2 0 0 0.8];%Gas
alfa = 0.5;
cells = zeros(cellss,length(xoil),contacts);
ki = zeros(cellss,length(xoil),contacts);
for ii=1:contacts
cells(1,1:4,ii) = ygas;
cells(2*ii,1:4,ii) = xoil;
end
for ii=2:contacts
    for i=2:2:ii*2-1
        z=cells(i,1:4,ii-1)+alfa.*(cells(i-1,1:4,ii-1)-cells(i,1:4,ii-1));
        [x,y,k]=flash2(z);
% [x,y]=meshgrid(1:4,1,1);
% x=ii.*[1,2,3,4];
% y=5*i.*[1,2,3,4];
% k=ii.*[1,2,3,4];
        ki(i,1:4,ii)=k;
        cells(i,1:4,ii)=x;
        cells(i+1,1:4,ii)=y;
    end
end

for j=2:2:cellss-2
Ln(j-1) = sqrt(sum((cells(j,1:4,contacts)-cells(j-1,1:4,contacts)).^2));
end

plot(Ln(3:2:end),'.','Markersize',20);
ylim([0 .8]);
hold on
end
xlabel('Cell Number');
ylabel('Tie-line length');
legend('50 Cells','100 Cells','150 Cells',' 250 Cells');
hold off

%% TO CREATE FIGURE 3
% clear all
figure('Color','w');set(gcf,'Position',[553   -49   590   354])
j = 5;
h = plot(ki(2:2*j:end-4,1,contacts)+.21,'--','Markersize',10,'LineWidth',2);
hold on
plot(ki(2:2*j:end-4,2,contacts)./1.6+.07,'--','Markersize',10,'LineWidth',2);
plot(ki(2:2*j:end-4,3,contacts)./2.3,'--','Markersize',10,'LineWidth',2);
plot(ki(2:2*j:end-4,4,contacts)-.1,'--','Markersize',10,'LineWidth',2);
xlabel('Cell Number');
ylabel('K-Value');
ax=gca;
ax.XTickLabel=[0:25:125];
text(22,2.35,'C_1')
text(22,1.3,'CO_2')
text(22,.4,'C_4')
text(22,.1,'C_1_0')
% pcolor(squeeze(cells(:,4,:)));
% title(['CO2 advance with ',num2str(cellss-2),'cells']);
% xlabel('Contacts');
% ylabel('Cells');
% cells(1:end,1:4,end);
hold off

%% TO CREATE FIGURES 4 FOR OIL, GAS AND ALL COMPONENTS
clear all
iii=1;
figure('Color','w');set(gcf,'Position',[553   -49   590   354])
for P = [500  550 600 750 800 850 900 1000 1250 1400 1700 2000 2030 2040 2050 2055 2060 2065 2070 2075 2080 ];%2090 2095 2099 2012 2013]
    cellss = 252;
    contacts = (cellss-2)/2;
    xoil = [0.2 0.15 0.65 0];% Oil
    ygas = [0.2 0 0 0.8];%Gas
    alfa = 0.5;
    cells = zeros(cellss,length(xoil),contacts);
    ki = zeros(cellss,length(xoil),contacts);
    for ii=1:contacts
        cells(1,1:4,ii) = ygas;
        cells(2*ii,1:4,ii) = xoil;
    end
    for ii=2:contacts
        for i=2:2:ii*2-1
            z=cells(i,1:4,ii-1)+alfa.*(cells(i-1,1:4,ii-1)-cells(i,1:4,ii-1));
            [x,y,k]=flash3(z,P); %Flash3 is the same as flash 2 but with Pressure as an input
            %[x,y]=meshgrid(1:4,1,1); Testing Matrix is working
            % x=ii.*[1,2,3,4];
            % y=5*i.*[1,2,3,4];
            % k=ii.*[1,2,3,4];
            ki(i,1:4,ii)=k;
            cells(i,1:4,ii)=x;
            cells(i+1,1:4,ii)=y;
%             pause
            
        end
    end
    for j=2:2:cellss-2
        Ln1(j-1) = sqrt((cells(j,1,contacts)-cells(j-1,1,contacts)).^2);
        Ln2(j-1) = sqrt((cells(j,2,contacts)-cells(j-1,2,contacts)).^2);
        Ln3(j-1) = sqrt((cells(j,3,contacts)-cells(j-1,3,contacts)).^2);
        Ln4(j-1) = sqrt((cells(j,4,contacts)-cells(j-1,4,contacts)).^2);
    end
    
%     [aa,bb]=findpeaks(-Ln(3:2:end-4));
%     [aam,bbm]=min(Ln(3:2:end-4));
    plot(P,Ln1(123),'.','Color',[1 0 0],'Markersize',10);
    hold on
    plot(P,Ln2(123),'.','Color',[0 1 0],'Markersize',10);
    plot(P,Ln3(123),'.','Color',[0 0 1],'Markersize',10);
    plot(P,Ln4(123),'.','Color',[1 0 1],'Markersize',10);
    

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

matr1(iii)=[Ln1(123)];
matr2(iii)=[Ln2(123)];
matr3(iii)=[Ln3(123)];
matr4(iii)=[Ln4(123)];
iii=iii+1;
end  
% pp=polyfit([500 750 1250 1400 1700 2000 2050 2060 2070 2080 2074 2303],[matr,.2, 0],6);
% plot( linspace(1700,2303,1000),polyval(pp,linspace(1700,2303,1000)),'-','Color',[0 .6 0],'LineWidth',2)
xlim([500 2500]);
ylim([0 .65]);
xlabel('Pressure (psi)');
ylabel('Tie-line length');
text(1550,.21,'CO_2')
text(1550,.28,'C_1_0')
text(1550,.11,'C_1')
text(1550,.03,'C_4')

pp1=polyfit([500  550 600 750 800 850 900 1000 1250 1400 1700 2000 2030 2040 2050 2055 2060 2065 2070 2075 2080 2303],[matr1,0],6);
plot( linspace(1700,2303,1000),polyval(pp1,linspace(1700,2303,1000)),'-','Color',[.8 0 0],'LineWidth',.9)
pp2=polyfit([500  550 600 750 800 850 900 1000 1250 1400 1700 2000 2030 2040 2050 2055 2060 2065 2070 2075 2080 2303],[matr2,0],6);
plot( linspace(1700,2303,1000),polyval(pp2,linspace(1700,2303,1000)),'-','Color',[0 .8 0],'LineWidth',.9)
pp3=polyfit([500  550 600 750 800 850 900 1000 1250 1400 1700 2000 2030 2040 2050 2055 2060 2065 2070 2075 2080 2303],[matr3,0],6);
plot( linspace(1700,2303,1000),polyval(pp3,linspace(1700,2303,1000)),'-','Color',[0 0 .8],'LineWidth',.9)
pp4=polyfit([500  550 600 750 800 850 900 1000 1250 1400 1700 2000 2030 2040 2050 2055 2060 2065 2070 2075 2080 2303],[matr4,0],6);
plot( linspace(1700,2303,1000),polyval(pp4,linspace(1700,2303,1000)),'-','Color',[.8 0 .8],'LineWidth',.9)

hold off
