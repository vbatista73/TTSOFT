function [ m ] = readparameters()
% normal mode frequences
freq(1)=0.0739;
freq(2)=0.1258;
freq(3)=0.1525;
freq(4)=0.1961;
freq(5)=0.3788;
freq(6)=0.1139;
freq(7)=0.0937;
freq(8)=0.1219;
freq(9)=0.0873;
freq(10)=0.1669;
freq(11)=0.1891;
freq(12)=0.3769;
freq(13)=0.0423;
freq(14)=0.1190;
freq(15)=0.1266;
freq(16)=0.1408;
freq(17)=0.1840;
freq(18)=0.3734;
freq(19)=0.1318;
freq(20)=0.1425;
freq(21)=0.1756;
freq(22)=0.3798;
freq(23)=0.0521;
freq(24)=0.0973;
for i=1:24
    m(i)=27.2/freq(i);
end
end

