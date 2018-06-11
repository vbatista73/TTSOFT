function [a]=sam_plus(b,c)
% Performs addition on single cores from each of the TT tensors being added
d=b.d;
nb=b.n;
nc=c.n;
rb=b.r;
rc=c.r;
r=rb+rc;
posb=b.ps;
posc=c.ps;
crb=b.core;
crc=c.core;
psa=zeros(1,d+1);
psa(1)=1;
if rb(1)==1 && rc(1)==1
    core1b=crb(1:posb(2)-1);
    core1b=reshape(core1b,[nb(1) rb(2)]);
    core1c=crc(1:posc(2)-1);
    core1c=reshape(core1c,[nc(1) rc(2)]);
    if rb(2)==1 && rc(2)==1
        core1=core1b+core1c;
    else
        core1=[core1b,core1c];
    end
    cra=core1(:);
    psa(2)=1+nb(1)*r(2);
    start=2;
    r(1)=1;
else
    start=1;
end
if rb(d+1)==1 && rc(d+1)==1
    fin=d-1;
    r(d+1)=1;
else
    fin=d;
end
for j=start:fin
    corejb=crb(posb(j):posb(j+1)-1);
    corejb=reshape(corejb,[rb(j),nb(j),rb(j+1)]);
    corejc=crc(posc(j):posc(j+1)-1);
    corejc=reshape(corejc,[rc(j),nc(j),rc(j+1)]);
    corej=corejb; % this needs to be fixed for rank-1 cases
    corej(rb(j)+1:r(j),:,rb(j+1)+1:r(j+1))=corejc;
    psa(j+1)=psa(j)+nb(j)*r(j)*r(j+1);
    cra(psa(j):psa(j+1)-1)=corej(:);
end
if fin==d-1 && (d>1 || r(1)>1)
    coredb=crb(posb(d):posb(d+1)-1);
    coredb=reshape(coredb,[rb(d),nb(d)]);
    coredc=crc(posc(d):posc(d+1)-1);
    coredc=reshape(coredc,[rc(d),nc(d)]);
    cored=[coredb;coredc];
    psa(d+1)=psa(d)+nb(d)*r(d);
    cra(psa(d):psa(d+1)-1)=cored(:);
end
a=tt_tensor;
a.d=d;
a.n=nb;
a.r=r;
a.ps=psa;
a.core=cra;
end