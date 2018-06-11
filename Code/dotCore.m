function [p]=dotCore(p,tt1,tt2)
% Takes the conjugate dot product of cores of two tensors as a step in computing dot
% product of whole tensor

cr1=tt1.core;
cr2=tt2.core;
cr2=reshape(cr2,[tt2.r(1),tt2.n(1)*tt2.r(2)]);

p=p*cr2;
p=reshape(p,1,tt1.r(1)*tt1.n(1),tt2.r(2));
p=permute(p,[1,3,2]);
p=reshape(p,tt2.r(2),tt1.r(1)*tt1.n(1));

cr1=reshape(cr1,[tt1.r(1)*tt1.n(1),tt1.r(2)]);

p=p*conj(cr1);
p=reshape(p,1,tt2.r(2),tt1.r(2));
p=permute(p,[1,3,2]);
p=reshape(p,tt1.r(2),tt2.r(2));

end