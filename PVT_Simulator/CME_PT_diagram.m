function [pointTPTlow,pointPPTlow]=CME_PTdiagram(zi10,Tci10,Pci10,wi10)


%% *****************************Start**************************************
% Section 
% P-T Diagram
%
toc
TPT=linspace(-45,490,40)+273.15;
PPT=linspace(0.05,23.5,40).*10^6;%
lPT=zeros(length(PPT),length(TPT));
for ii = 1:length(TPT)
       for j=1:length(PPT);
             [xi10,yi10,Ki10,VLii10,VVii10,lii10,vii10]=FLASH_A_Astro(PPT(j),TPT(ii),zi10,Tci10,Pci10,wi10);
        lPT(j,ii) = lii10;
        lii10
       end
end
        toc
figure('Color','w');
surf(real(lPT));
view(2)
xlabel('Temperature (°C)')
ylabel('Pressure (MPa)')
title('P-T diagram')
set(gca,'XTickLabel',round(linspace(21,490,8).*100)./100);
set(gca,'YTickLabel',round(linspace(3,23.5,9).*10)./10);
axis tight
%%
[aa,bb]=ind2sub(size(lPT),find(lPT==1));
jj=1;jjj=1;
for i=1:length(bb)
    if bb(i)==jj
        pointP(jjj)=aa(i);
        jjj=jjj+1;
        jj=jj+1;
    end
end
PPT=linspace(0.05,23.5,40).*10^6;
pointPPT=PPT(pointP(1:31));
TPT=linspace(-45,490,40)+273.15;
pointTPT=TPT(1:31);
[rrr]=polyfit(pointTPT,pointPPT,2);
figure
plot(pointTPT,pointPPT,'o')
hold on
[rr]=polyval(rrr,linspace(-45,490,100)+273.15);
% plot(linspace(-45,490,100)+273.15,rr)
%%
[aa,bb]=ind2sub(size(lPT),find(lPT==0));
jj=34;jjj=1;
for i=1:length(bb)
    if bb(i)==jj
        if aa(i)>12
        pointPlow(jjj)=aa(i);
        jjj=jjj+1;
        jj=jj+1;
        end
    end
end
PPT=linspace(0.05,23.5,40).*10^6;
pointPPTlow=(PPT(pointPlow));
TPT=linspace(-45,490,40)+273.15;
pointTPTlow=TPT(34:40);
[rrrlow]=polyfit(pointTPTlow,pointPPTlow,4);
[rrlow]=polyval(rrrlow,linspace(320,490,100)+273.15);
% plot(linspace(320,490,100)+273.15,rrlow)
plot(pointTPTlow,pointPPTlow,'o')
[test]=polyfit([pointTPT,pointTPTlow],[pointPPT,pointPPTlow],3);
[rrtest]=polyval(test,linspace(-45,490,100)+273.15);
plot(linspace(-45,490,100)+273.15,rrtest,'k','LineWidth',3)
xlabel('Temperature [K]')
ylabel('Pressure [Pa]')
title('P-T diagram of 10')
grid on
toc

end