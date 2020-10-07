% Phase diagram
TPT=linspace(-45,490,40)+273.15;
PPT=linspace(0.05,23.5,40).*10^6;%
lPT=zeros(length(PPT),length(TPT));
for ii = 1:length(TPT)
       for j=1:length(PPT);
             [xi10,yi10,Ki10,VLii10,VVii10,lii10,vii10]=FLASH_A_Astro(PPT(j),TPT(ii),zi10,Tci10,Pci10,wi10);
        lPT(j,ii) = lii10;
        lii10;
       end
end
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
pointPPT=PPT(pointP(1:end));
TPT=linspace(-45,490,40)+273.15;
pointTPT=TPT(1:length(pointPPT));
[rrr]=polyfit(pointTPT,pointPPT,2);
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
[test]=polyfit([pointTPT,pointTPTlow],[pointPPT,pointPPTlow],3);
[rrtest]=polyval(test,linspace(-45,490,100)+273.15);

