function [v1,v2,vc]=dimV4D(xVec,d)
% Calculates v1, v2, vc in TT format for dimers
% parameters from paper
% *2 included in diabatic terms

x_tt=tt_tensor(xVec);
nx=size(xVec,2);

ag=5;
b1g=1;
b2g=2;
b3g=4;
au=2;
b1u=4;
b2u=4;
b3u=2;
total=ag+b1g+b2g+b3g+au+b1u+b2u+b3u;

v=zeros(2);
delta=0.4230;

wag=zeros(1,ag); 
wau=zeros(1,au); 
wb1g=zeros(1,b1g); 
wb2g=zeros(1,b2g); 
wb3g=zeros(1,b3g); 
wb1u=zeros(1,b1u); 
wb2u=zeros(1,b2u); 
wb3u=zeros(1,b3u); 

aiag=zeros(1,ag); 
biag=zeros(1,ag); 

aijag=zeros(ag,ag)*NaN; 
bijag=zeros(ag,ag)*NaN;
aijau=zeros(au,au)*NaN;
bijau=zeros(au,au)*NaN;
aijb1g=zeros(b1g,b1g)*NaN; 
bijb1g=zeros(b1g,b1g)*NaN; 
aijb2g=zeros(b2g,b2g)*NaN;
bijb2g=zeros(b2g,b2g)*NaN;
aijb3g=zeros(b3g,b3g)*NaN;
bijb3g=zeros(b3g,b3g)*NaN;
aijb1u=zeros(b1u,b1u)*NaN;
bijb1u=zeros(b1u,b1u)*NaN;
aijb2u=zeros(b2u,b2u)*NaN;
bijb2u=zeros(b2u,b2u)*NaN;
aijb3u=zeros(b3u,b3u)*NaN;
bijb3u=zeros(b3u,b3u)*NaN;

cijb1gag=zeros(b1g,ag)*NaN;
cijb2gb3g=zeros(b2g,b3g)*NaN;
cijaub1u=zeros(au,b1u)*NaN;
cijb3ub2u=zeros(b3u,b2u)*NaN;

cib1g=zeros(1,b1g); 
cib1g(1)= 0.2080; % 10a, lambda in Meyer's

%     frequency;

wag(1)=0.0739;
wag(2)=0.1258;
wag(3)=0.1525;
wag(4)=0.1961;
wag(5)=0.3788;

wb1g(1)=0.1139;

wb2g(1)=0.0937;
wb2g(2)=0.1219;

wb3g(1)=0.0873;
wb3g(2)=0.1669;
wb3g(3)=0.1891;
wb3g(4)=0.3769;

wau(1)=0.0423;
wau(2)=0.1190;

wb1u(1)=0.1266;
wb1u(2)=0.1408;
wb1u(3)=0.1840;
wb1u(4)=0.3734;

wb2u(1)=0.1318;
wb2u(2)=0.1425;
wb2u(3)=0.1756;
wb2u(4)=0.3798;

wb3u(1)=0.0521;
wb3u(2)=0.0973;

%     linear, on-diagonal coupling coefficients mrci;

%      aiag(1)=-0.0964; % 6a
%      aiag(2)=-0.0470; % 1
%      aiag(3)=0.1594; % 9a
%      aiag(4)=-0.0623; % 8a
%      aiag(5)=0.0368; % 2
%      biag(1)=0.1193; % 6a
%      biag(2)=-0.2012; % 1
%      biag(3)=0.0484; % 9a
%      biag(4)=0.0348; % 8a
%      biag(5)=0.0211; % 2

%     linear, on-diagonal coupling coefficients casscf

% aiag(1)=-0.0938; % 6a
% aiag(2)=-0.0808; % 1
% aiag(3)=0.1511; % 9a
% aiag(4)=0.0193; % 8a
% aiag(5)=0.0281; % 2
aiag(1)=-0.09806;           % 6a
aiag(2)=-0.05033;           % 1
aiag(3)=0.14521;
aiag(4)=-0.0445;
aiag(5)=0.0247;

% biag(1)=0.1029; % 6a
% biag(2)=-0.2234; % 1
% biag(3)=0.0509; % 9a
% biag(4)=0.0248; % 8a
% biag(5)=0.0340; % 2
biag(1)=0.13545;          % 6a
biag(2)=-0.17100;           % 1
% biag(2)=-0.1629;
biag(3)=0.03746;           % 9a
biag(4)=0.0168;
biag(5)=0.0162;

%     quadratic, on-diagonal coupling coefficients;

% aijag(1,1)=0.00002; % 6a
% aijag(1,2)=0.00108; % 6ax1
aijag(1,3)=-0.00204; % 6ax9a
 aijag(1,1)=0.0000;          % 6a                                                                                   
      aijag(1,2)=0.00108;        % 6ax1                                                                                  
%       aijag(1,3)=0.00204;        % 6ax9a         

aijag(1,4)=-0.00135; % 6ax8a
aijag(1,5)=-0.00285; % 6ax2

% aijag(2,2)=-0.00810; % 1
aijag(2,3)=0.00474; % 1x9a
      aijag(2,2)=0.0000;          % 1                                                                                    
%       aijag(2,3)=-0.00474;        % 1x9a
      
aijag(2,4)=0.00154; % 1x8a
aijag(2,5)=-0.00163; % 1x2

% aijag(3,3)=-0.00116; % 9a
 aijag(3,3)=0.0000;          % 9a
 
aijag(3,4)=0.00872; % 9ax8a
aijag(3,5)=-0.00474; % 9ax2
% aijag(4,4)=0.00028; % 8a
aijag(4,4)=0;
aijag(4,5)=-0.00143; % 8ax2
% aijag(5,5)=-0.00118; % 2
aijag(5,5)=0;

% bijag(1,1)=-0.00917; % 6a
% bijag(1,2)=-0.00298; % 6ax1
bijag(1,3)=-0.00189; % 6ax9a
      bijag(1,1)=0.00000;           % 6a                                                                                 
      bijag(1,2)=-0.00298;       % 6ax1                                                                                  
%       bijag(1,3)=0.00189;        % 6ax9a     

bijag(1,4)=-0.00203; % 6ax8a
bijag(1,5)=-0.00128; % 6ax2

% bijag(2,2)=0.00488; % 1
bijag(2,3)=0.00155; % 1x9a
      bijag(2,2)=0.00000;            % 1                                                                                 
%       bijag(2,3)=-0.00155; 

bijag(2,4)=0.00311; % 1x8a
bijag(2,5)=-0.00600; % 1x2

% bijag(3,3)=0.00022; % 9a
     bijag(3,3)=0.00000;            % 9a  
     
bijag(3,4)=0.01194; % 9ax8a
bijag(3,5)=-0.00334; % 9ax2
% bijag(4,4)=0.01272; % 8a
bijag(4,4)=0;
bijag(4,5)=-0.00713; % 8ax2
% bijag(5,5)=0.00039; % 2
bijag(5,5)=0;

% aijau(1,1)=0.01626; % 16a
aijau(1,1)=0.01145;
aijau(1,2)=0.00100; % 16ax17a
aijau(2,2)=-0.02040; % 17a

% bijau(1,1)=-0.02026; % 16a
bijau(1,1)=-0.01459;
bijau(1,2)=-0.00091; % 16ax17a
bijau(2,2)=-0.00618; % 17a

aijb1g(1,1)=-0.01159; % 10a

bijb1g(1,1)=-0.01159; % 10a

% aijb2g(1,1)=-0.02427; % 4
aijb2g(1,1)=-0.02252;
aijb2g(1,2)=-0.00049; % 4x5
aijb2g(2,2)=-0.01825; % 5

% bijb2g(1,1)=-0.04439; % 4
bijb2g(1,1)=-0.03445;
bijb2g(1,2)= 0.00911; % 4x5
bijb2g(2,2)=-0.00265; % 5


aijb3g(1,1)=-0.00741; % 6b
aijb3g(1,2)= 0.01321; % 6bx3
aijb3g(1,3)=-0.00717; % 6bx8b
aijb3g(1,4)= 0.00515; % 6bx7b
aijb3g(2,2)= 0.05183; % 3
aijb3g(2,3)=-0.03942; % 3x8b
aijb3g(2,4)= 0.00170; % 3x7b
aijb3g(3,3)=-0.05733; % 8b
aijb3g(3,4)=-0.00204; % 8bx7b
aijb3g(4,4)=-0.00333; % 7b

bijb3g(1,1)=-0.00385; % 6b
bijb3g(1,2)=-0.00661; % 6bx3
bijb3g(1,3)= 0.00429; % 6bx8b
bijb3g(1,4)=-0.00246; % 6bx7b
bijb3g(2,2)= 0.04842; % 3
bijb3g(2,3)=-0.03034; % 3x8b

bijb3g(2,4)=-0.00185; % 3x7b
bijb3g(3,3)=-0.06332; % 8b
bijb3g(3,4)=-0.00388; % 8bx7b
bijb3g(4,4)=-0.00040; % 7b

aijb1u(1,1)=-0.04819; % 12
aijb1u(1,2)= 0.00525; % 12x18a
aijb1u(1,3)=-0.00485; % 12x19a
aijb1u(1,4)=-0.00326; % 12x13
aijb1u(2,2)=-0.00792; % 18a
aijb1u(2,3)= 0.00852; % 18ax19a
aijb1u(2,4)= 0.00888; % 18ax13
aijb1u(3,3)=-0.02429; % 19a
aijb1u(3,4)=-0.00443; % 19ax13
aijb1u(4,4)=-0.00492; % 13

bijb1u(1,1)=-0.00840; % 12
bijb1u(1,2)= 0.00536; % 12x18a
bijb1u(1,3)=-0.00097; % 12x19a
bijb1u(1,4)= 0.00034; % 12x13
bijb1u(2,2)= 0.00429; % 18a
bijb1u(2,3)= 0.00209; % 18ax19a
bijb1u(2,4)=-0.00049; % 18ax13
bijb1u(3,3)=-0.00734; % 19a
bijb1u(3,4)= 0.00346; % 19ax13
bijb1u(4,4)= 0.00062; % 13

% aijb2u(1,1)=-0.02577; % 18b
aijb2u(1,1)=-0.00277;
% aijb2u(1,2)=-0.01115; % 18bx14
aijb2u(1,2)=0.00016;
aijb2u(1,3)=-0.00250; % 18bx19b
aijb2u(1,4)= 0.00357; % 18bx20b
% aijb2u(2,2)= 0.01711; % 14
aijb2u(2,2)=0.03924;
aijb2u(2,3)=-0.00197; % 14x19b
aijb2u(2,4)=-0.00355; % 14x20b
aijb2u(3,3)= 0.00992; % 19b
aijb2u(3,4)= 0.00623; % 19bx20b
aijb2u(4,4)=-0.00110; % 20b

% bijb2u(1,1)= 0.32670; % 18b
bijb2u(1,1)=-0.01179;
% bijb2u(1,2)=-0.26747; % 18bx14
bijb2u(1,2)=-0.00844;
% bijb2u(1,3)= 0.07843; % 18bx19b
bijb2u(1,3)=0.07000;
bijb2u(1,4)=-0.01249; % 18bx20b
% bijb2u(2,2)= 0.17948; % 14
bijb2u(2,2)=0.04000;
% bijb2u(2,3)=-0.06091; % 14x19b
bijb2u(2,3)=-0.05000;
bijb2u(2,4)= 0.00265; % 14x20b
bijb2u(3,3)= 0.01246; % 19b
bijb2u(3,4)=-0.00422; % 19bx20b
bijb2u(4,4)= 0.00069; % 20b

aijb3u(1,1)=-0.02176; % 16b
aijb3u(1,2)=-0.00624; % 16bx11
aijb3u(2,2)= 0.00315; % 11

bijb3u(1,1)=-0.02214; % 16b
bijb3u(1,2)=-0.00261; % 16bx11
bijb3u(2,2)=-0.00496; % 11

% cijb1gag(1,1)= 0.00628; % 10ax6a
      cijb1gag(1,1)= -0.0100;    % 10ax6a 
cijb1gag(1,2)=-0.00551; % 10ax1
%       cijb1gag(1,2)=0.00553;    % 10ax1  
cijb1gag(1,3)= 0.00127; % 10ax9a
%       cijb1gag(1,3)= 0.00126;    % 10ax9a   
cijb1gag(1,4)= 0.00799; % 10ax8a
cijb1gag(1,5)=-0.00512; % 10ax2

cijb2gb3g(1,1)=-0.01372; % 4x6b
cijb2gb3g(1,2)=-0.00466; % 4x3
cijb2gb3g(1,3)= 0.00329; % 4x8b
cijb2gb3g(1,4)=-0.00031; % 4x7b
cijb2gb3g(2,1)= 0.00598; % 5x6b
cijb2gb3g(2,2)=-0.00914; % 5x3
cijb2gb3g(2,3)= 0.00961; % 5x8b
cijb2gb3g(2,4)= 0.00500; % 5x7b

cijaub1u(1,1)=-0.01056; % 16ax12
cijaub1u(1,2)= 0.00559; % 16ax18a
cijaub1u(1,3)= 0.00401; % 16ax19a
cijaub1u(1,4)=-0.00226; % 16ax13
cijaub1u(2,1)=-0.01200; % 17ax12
cijaub1u(2,2)=-0.00213; % 17ax18a
cijaub1u(2,3)= 0.00328; % 17ax19a
cijaub1u(2,4)=-0.00396; % 17ax13

cijb3ub2u(1,1)= 0.00118; % 16bx18b
cijb3ub2u(1,2)=-0.00009; % 16bx14
cijb3ub2u(1,3)=-0.00285; % 16bx19b
cijb3ub2u(1,4)=-0.00095; % 16bx20b
cijb3ub2u(2,1)= 0.01281; % 11x18b
cijb3ub2u(2,2)=-0.01780; % 11x14
cijb3ub2u(2,3)= 0.00134; % 11x19b
cijb3ub2u(2,4)=-0.00481; % 11x20b

%-------------------End of parameters ----------------------------------------;
%     ag b1g b2g b3g gu b1u b2u b3u

%------------------------- v1 v2-------------------------
v1=tt_zeros(nx,d);
v2=v1;

ag=3;
dimsB4=0;
% ag x ag
for i=1:min(ag,d-dimsB4)
    if i-1>0
        poti=tt_ones(nx,i-1);
    else
        poti=[];
    end
    poti=tkron(poti,x_tt);
    for j=i+1:min(ag,d-dimsB4)
        if j-i>1
            potj=tt_ones(nx,j-i-1);
        else
            potj=[];
        end
        potj1=tkron(potj,aijag(i,j)*x_tt);
        potj2=tkron(potj,bijag(i,j)*x_tt);
        if d-j-dimsB4>0
            ott=tt_ones(nx,d-j-dimsB4);
            potj1=tkron(potj1,ott);
            potj2=tkron(potj2,ott);
        end
        v1=v1+tkron(poti,potj1);
        v2=v2+tkron(poti,potj2);
    end
end
dimsB4=dimsB4+ag;

% no b1g x b1g term
dimsB4=dimsB4+b1g;

% one b2g x b2g term
if d>=dimsB4+b2g
    poti=tt_ones(nx,dimsB4);
    poti=tkron(poti,x_tt);
    potj1=aijb2g(1,2)*x_tt;
    potj2=bijb2g(1,2)*x_tt;
    if d>dimsB4+b2g
        ott=tt_ones(nx,d-dimsB4-b2g);
        potj1=tkron(potj1,ott);
        potj2=tkron(potj2,ott);
    end
    v1=v1+tkron(poti,potj1);
    v2=v2+tkron(poti,potj2);
end
dimsB4=dimsB4+b2g;

% b3g x b3g
for i=1:min(b3g,d-dimsB4)
    poti=tt_ones(nx,i-1+dimsB4);
    poti=tkron(poti,x_tt);
    for j=i+1:min(b3g,d-dimsB4)
        if j-i>1
            potj=tt_ones(nx,j-i-1);
        else
            potj=[];
        end
        potj1=tkron(potj,aijb3g(i,j)*x_tt);
        potj2=tkron(potj,bijb3g(i,j)*x_tt);
        if d-j-dimsB4>0
            ott=tt_ones(nx,d-j-dimsB4);
            potj1=tkron(potj1,ott);
            potj2=tkron(potj2,ott);
        end
        v1=v1+tkron(poti,potj1);
        v2=v2+tkron(poti,potj2);
    end
end
dimsB4=dimsB4+b3g;

% one au x au term
if d>=dimsB4+au
    poti=tt_ones(nx,dimsB4);
    poti=tkron(poti,x_tt);
    potj1=aijau(1,2)*x_tt;
    potj2=bijau(1,2)*x_tt;
    if d>dimsB4+au
        ott=tt_ones(nx,d-dimsB4-au);
        potj1=tkron(potj1,ott);
        potj2=tkron(potj2,ott);
    end
    v1=v1+tkron(poti,potj1);
    v2=v2+tkron(poti,potj2);
end
dimsB4=dimsB4+au;

% b1u x b1u
for i=1:min(b1u,d-dimsB4)
    poti=tt_ones(nx,i-1+dimsB4);
    poti=tkron(poti,x_tt);
    for j=i+1:min(b1u,d-dimsB4)
        if j-i>1
            potj=tt_ones(nx,j-i-1);
        else
            potj=[];
        end
        potj1=tkron(potj,aijb1u(i,j)*x_tt);
        potj2=tkron(potj,bijb1u(i,j)*x_tt);
        if d-j-dimsB4>0
            ott=tt_ones(nx,d-j-dimsB4);
            potj1=tkron(potj1,ott);
            potj2=tkron(potj2,ott);
        end
        v1=v1+tkron(poti,potj1);
        v2=v2+tkron(poti,potj2);
    end
end
dimsB4=dimsB4+b1u;

% b2u x b2u
for i=1:min(b2u,d-dimsB4)
    poti=tt_ones(nx,i-1+dimsB4);
    poti=tkron(poti,x_tt);
    for j=i+1:min(b2u,d-dimsB4)
        if j-i>1
            potj=tt_ones(nx,j-i-1);
        else
            potj=[];
        end
        potj1=tkron(potj,aijb2u(i,j)*x_tt);
        potj2=tkron(potj,bijb2u(i,j)*x_tt);
        if d-j-dimsB4>0
            ott=tt_ones(nx,d-j-dimsB4);
            potj1=tkron(potj1,ott);
            potj2=tkron(potj2,ott);
        end
        v1=v1+tkron(poti,potj1);
        v2=v2+tkron(poti,potj2);
    end
end
dimsB4=dimsB4+b2u;

% one b3u x b3u term
if d>=dimsB4+b3u
    poti=tt_ones(nx,dimsB4);
    poti=tkron(poti,x_tt);
    potj1=aijb3u(1,2)*x_tt;
    potj2=bijb3u(1,2)*x_tt;
    if d>dimsB4+b3u
        ott=tt_ones(nx,d-dimsB4-b3u);
        potj1=tkron(potj1,ott);
        potj2=tkron(potj2,ott);
    end
    v1=v1+tkron(poti,potj1);
    v2=v2+tkron(poti,potj2);
end

v1=v1*3.6749d-2;
v2=v2*3.6749d-2;

%----------------------- End of v1,v2 -------------------------

%----------------------- vc --------------------------------------
vc=tt_zeros(nx,d);

% ag(5) x b1g(1)
end1=0;
end2=end1+ag;
for a=1:ag
    if a-1>0
        poti=tt_ones(nx,a-1);
    else
        poti=[];
    end
    poti=tkron(poti,x_tt);
    for b=1:min(b1g,d-end2)
        if a+1<ag+b
            potj=tt_ones(nx,ag+b-a-1);
        else
            potj=[];
        end
        potj=tkron(potj,cijb1gag(b,a)*x_tt);
        if d>end2+b
            ott=tt_ones(nx,d-end2-b);
            potj=tkron(potj,ott);
        end
        vc=vc+tkron(poti,potj);
    end
end

% b2g(2) x b3g(4)
end1=end2+b1g;
end2=end1+b2g;
for a=1:b2g
    poti=tt_ones(nx,a-1+end1);
    poti=tkron(poti,x_tt);
    for b=1:min(b3g,d-end2)
        if a+1<b2g+b
            potj=tt_ones(nx,b2g+b-a-1);
        else
            potj=[];
        end
        potj=tkron(potj,cijb2gb3g(a,b)*x_tt);
        if d>end2+b
            ott=tt_ones(nx,d-end2-b);
            potj=tkron(potj,ott);
        end
        vc=vc+tkron(poti,potj);
    end
end

% au(2) x b1u(4)
end1=end2+b3g;
end2=end1+au;
for a=1:au
    poti=tt_ones(nx,a-1+end1);
    poti=tkron(poti,x_tt);
    for b=1:min(b1u,d-end2)
        if a+1<au+b
            potj=tt_ones(nx,au+b-a-1);
        else
            potj=[];
        end
        potj=tkron(potj,cijb2gb3g(a,b)*x_tt);
        if d>end2+b
            ott=tt_ones(nx,d-end2-b);
            potj=tkron(potj,ott);
        end
        vc=vc+tkron(poti,potj);
    end
end

% b2u(2) x b3u(4)
end1=end2+b1u;
end2=end1+b2u;
for a=1:b2u
    poti=tt_ones(nx,a-1+end1);
    poti=tkron(poti,x_tt);
    for b=1:min(b3u,d-end2)
        if a+1<b2u+b
            potj=tt_ones(nx,b2u+b-a-1);
        else
            potj=[];
        end
        potj=tkron(potj,cijb2gb3g(b,a)*x_tt);
        if d>end2+b
            ott=tt_ones(nx,d-end2-b);
            potj=tkron(potj,ott);
        end
        vc=vc+tkron(poti,potj);
    end
end

vc=vc*3.6749d-2;

end