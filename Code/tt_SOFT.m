% tt-SOFT propagation of non-adiabatic dynamics in pyrazine
% d-dimensional wavepacket on ne=2 electronic surfaces

clear all;
% Uncomment this line if running locally
% toolPath='~/Dropbox/Kenny/TT/newivan/TT-Toolbox';
% Uncomment this line if running on Omega
toolPath='~/scratch/TT-Toolbox';

numCo=11;
d=24; % number of dimensions
ns=620; % # of propagation steps
eps=6e-3; % precision for tt train round off

m = readparameters();
hbar = 1; % hbar
alpha=1; % initial width
x0= 0; % initial position
p0= 0; % initial momentum
q=4;   % quantics
nx=2^q;  % # of grid points in coord representation
np=nx; % # of grid points in momentum representation
L=14; % simulation box
dx= L/nx; % x grid spacing
dp= 2*pi/L/hbar; % spacing in momenta space
tau=10.0; % propagation time step
ne=2; % number of electronic states
Nexp=10; % number of terms in the exponential taylor series

% coordinate and momentum grids
x = ((1:nx)-nx/2)*dx;
p=fftshift(dp*(-np/2:np/2-1));
KE=p.^2/2;

addpath(genpath(toolPath));

% Set up parallel pool
pool=parpool('local',numCo);
addAttachedFiles(pool,toolPath);

dimCo=cell(1,numCo);
dimCo{1}=1:3;
for i=2:9
    dimCo{i}=i+2;
end
dimCo{10}=12:13;
dimCo{11}=14:24;

spmd
    dims=dimCo{labindex};
    nd=length(dims);
    
    % Kinetic propagator
    KP_tt=cell(1,nd);
    % Wave function
    Gx_tt=cell(ne,nd);
    for i=1:nd
        Gx_tt{1,i}=tt_zeros(nx);
        Gx_tt{2,i}=tt_tensor(exp(-alpha/2*(x-x0).^2+1i*p0*(x-x0)/hbar)*(alpha/pi)^0.25);
        KP_tt{i}=tt_tensor(exp(-1i*tau/2*KE/m(dims(i))/hbar));
    end
    G0_tt=Gx_tt(2,:);
end

% Potential
[v1m,v2m,vcm]=monV(x,d);
v1m=round(v1m,eps);
v2m=round(v2m,eps);
vcm=round(vcm,eps);
[v1d,v2d,vcd]=dimV(x,d);
v1d=round(v1d,eps);
v2d=round(v2d,eps);
vcd=round(vcd,eps);

% vmat represents the matrix to be exponentiated
vmat=cell(ne,ne);
vmat{1,1}=-1i*tau*round(v1m+v1d,eps)/hbar;
vmat{2,2}=-1i*tau*round(v2m+v2d,eps)/hbar;
vmat{1,2}=-1i*tau*round(vcm+vcd,eps)/hbar;

[v1max,~]=tt_max_abs(vmat{1,1});
[v2max,~]=tt_max_abs(vmat{2,2});
[vcmax,~]=tt_max_abs(vmat{1,2});
maxEl=max([v1max,v2max,vcmax]);
n0=floor(max(log2(maxEl),0))+1;
scale=1/(2^n0);
vmat{1,1}=vmat{1,1}*scale;
vmat{2,2}=vmat{2,2}*scale;
vmat{1,2}=vmat{1,2}*scale;
vmat{2,1}=vmat{1,2};

% Initialize y to identity matrix
y=cell(ne,ne);
y{1,1}=tt_ones(nx,d);
y{2,2}=y{1,1};
y{1,2}=tt_zeros(nx,d);
y{2,1}=y{1,2};

for l=(Nexp-1):-1:1
    % y=1+y*vmat/l
    res=cell(ne,ne);
    for i=1:ne
        for k=1:ne
            res{i,k}=tt_zeros(nx,d);
            for j=1:ne
                res{i,k}=res{i,k}+y{i,j}.*vmat{j,k};
            end
            res{i,k}=round(res{i,k}/l,eps);
        end
        res{i,i}=res{i,i}+tt_ones(nx,d);
    end
    y=res;
end

for l=1:n0
    % square y
    res=cell(ne,ne);
    for i=1:ne
        for k=1:ne
            res{i,k}=tt_zeros(nx,d);
            for j=1:ne
                res{i,k}=res{i,k}+y{i,j}.*y{j,k};
            end
            res{i,k}=round(res{i,k},eps*(0.5^(n0-l)));
        end
    end
    y=res;
end

PP_tt=cell(1,ne);
PP_tt{1}=y{1,1};
PP_tt{2}=y{2,2};
PA_tt=y{1,2};

PP_split=[cellSplit(PP_tt{1},ones(1,d)); cellSplit(PP_tt{2},ones(1,d))];
PA_tt_s=cellSplit(PA_tt,ones(1,d));

clear PP_tt PA_tt

spmd
    PP_tt=PP_split(:,dims);
    PA_tt=PA_tt_s(dims);
end

clear PP_split PA_tt_s

% Survival amplitude
amp=zeros(1,ns+1);

% Compute dot product
spmd
    if labindex==1
        p=1;
    else
        p=labReceive(labindex-1);
    end
    for i=1:nd
        p=dotCore(p,G0_tt{i},Gx_tt{2,i});
    end
    if labindex<numCo
        labSend(p,labindex+1);
    end
end
amp(1)=p{numCo}*dx^d;
clear p
ranks=zeros(ns,ne,d);

% Propagate psi
disp('starting propagation');
tic;
for k=1:ns
    disp(k);
    % Kinetic propagation
    spmd
        for j=1:ne
            for i=1:nd
                Gp_tt=Sam_ttFT_split(Gx_tt{j,i},1);
                Gp_tt=Gp_tt.*KP_tt{i};
                Gx_tt{j,i}=Sam_ttFT_split(Gp_tt,2);
            end
        end
        Gp_tt=[];
        
        % Potential propagation
        % Multiply
        GA_tt=cell(ne,nd);
        for j=1:ne
            for i=1:nd
                GA_tt{j,i}=Gx_tt{j,i}.*PA_tt{i};
                Gx_tt{j,i}=Gx_tt{j,i}.*PP_tt{j,i};
            end
        end
        
        % Round
        % QR L to R
        nrm=ones(ne,nd,2);
        norm1=zeros(ne,2);
        for j=1:ne
            for i=1:2
                if labindex==1
                    ru=1;
                else
                    ru=labReceive(labindex-1,j+2*i);
                end
                for l=1:nd
                    if i==1
                        currTT=Gx_tt{j,l};
                    else
                        currTT=GA_tt{j,l};
                    end
                    core0=currTT.core;
                    core0=reshape(core0,[currTT.r(1),nx*currTT.r(2)]);
                    core0=ru*core0;
                    newR1=size(core0,1);
                    if dims(l)<d
                        core0=reshape(core0,[newR1*nx,currTT.r(2)]);
                        [core0,ru]=qr(core0,0);
                        newR2=size(core0,2);
                        nrm(j,l,i)=norm(ru,'fro');
                        if(nrm(j,l,i)~=0)
                            ru=ru./nrm(j,l,i);
                        end
                    else
                        newR2=currTT.r(2);
                    end
                    if i==1
                        Gx_tt{j,l}.r(1)=newR1;
                        Gx_tt{j,l}.r(2)=newR2;
                        Gx_tt{j,l}.core=core0(:);
                    else
                        GA_tt{j,l}.r(1)=newR1;
                        GA_tt{j,l}.r(2)=newR2;
                        GA_tt{j,l}.core=core0(:);
                    end
                end
                if dims(nd)<d
                    labSend(ru,labindex+1,j+2*i);
                end
            end
        end
        % SVD R to L
        ep=eps/sqrt(d-1);
        for j=1:ne
            for i=1:2
                if labindex==numCo
                    u=1;
                else
                    u=labReceive(labindex+1,j+2*i);
                end
                for l=nd:-1:1
                    if i==1
                        currTT=Gx_tt{j,l};
                    else
                        currTT=GA_tt{j,l};
                    end
                    core0=currTT.core;
                    core0=reshape(core0,[currTT.r(1)*nx,currTT.r(2)]);
                    core0=core0*u;
                    newR2=size(core0,2);
                    if dims(l)>1
                        core0=reshape(core0,[currTT.r(1),nx*newR2]);
                        [u,s,v]=svd(core0,'econ');
                        evals=diag(s);
                        newR1=sum(evals>ep);
                        newR1=max(1,newR1);
                        u=u(:,1:newR1);
                        v=v(:,1:newR1);
                        s=s(1:newR1,1:newR1);
                        core0=v';
                        u=u*s;
                    else
                        newR1=currTT.r(1);
                        % norm for dim=1
                        norm1(j,i)=norm(core0(:),'fro');
                        if (norm1(j,i) ~= 0)
                            core0=core0./norm1(j,i);
                        end
                    end
                    if i==1
                        Gx_tt{j,l}.r(1)=newR1;
                        Gx_tt{j,l}.r(2)=newR2;
                        Gx_tt{j,l}.core=core0(:);
                        Gx_tt{j,l}.ps(2)=newR1*newR2*nx+1;
                    else
                        GA_tt{j,l}.r(1)=newR1;
                        GA_tt{j,l}.r(2)=newR2;
                        GA_tt{j,l}.core=core0(:);
                        GA_tt{j,l}.ps(2)=newR1*newR2*nx+1;
                    end
                end
                if dims(1)>1
                    labSend(u,labindex-1,j+2*i);
                end
            end
        end
    end
    clear core0 currTT ru u v s
    nrm0=zeros(ne,2);
    for j=1:ne
        for i=1:2
            for l=1:numCo
                tab=nrm{l};
                nrm0(j,i)=nrm0(j,i)+sum(log(abs(tab(j,:,i))));
            end
            tab=norm1{1};
            nrm0(j,i)=nrm0(j,i)+log(abs(tab(j,i)));
            nrm0(j,i)=exp(nrm0(j,i)/d);
        end
    end
    
    spmd
        for j=1:ne
            for i=1:nd
                Gx_tt{j,i}=Gx_tt{j,i}*nrm0(j,1);
                GA_tt{j,i}=GA_tt{j,i}*nrm0(j,2);
            end
        end
        
        % Add
        for i=1:nd
            Gx_tt{1,i}=sam_plus(Gx_tt{1,i},GA_tt{2,i});
            Gx_tt{2,i}=sam_plus(Gx_tt{2,i},GA_tt{1,i});
        end
        GA_tt=[];
        
        % Round
        % QR L to R
        nrm=ones(ne,nd);
        norm1=zeros(ne,1);
        for j=1:ne
            if labindex==1
                ru=1;
            else
                ru=labReceive(labindex-1,j);
            end
            for l=1:nd
                currTT=Gx_tt{j,l};
                core0=currTT.core;
                core0=reshape(core0,[currTT.r(1),nx*currTT.r(2)]);
                core0=ru*core0;
                newR1=size(core0,1);
                if dims(l)<d
                    core0=reshape(core0,[newR1*nx,currTT.r(2)]);
                    [core0,ru]=qr(core0,0);
                    newR2=size(core0,2);
                    nrm(j,l)=norm(ru,'fro');
                    if(nrm(j,l)~=0)
                        ru=ru./nrm(j,l);
                    end
                else
                    newR2=currTT.r(2);
                end
                Gx_tt{j,l}.r(1)=newR1;
                Gx_tt{j,l}.r(2)=newR2;
                Gx_tt{j,l}.core=core0(:);
            end
            if dims(nd)<d
                labSend(ru,labindex+1,j);
            end
        end
        % SVD R to L
        ep=eps/sqrt(d-1);
        for j=1:ne
            if labindex==numCo
                u=1;
            else
                u=labReceive(labindex+1,j);
            end
            for l=nd:-1:1
                currTT=Gx_tt{j,l};
                core0=currTT.core;
                core0=reshape(core0,[currTT.r(1)*nx,currTT.r(2)]);
                core0=core0*u;
                newR2=size(core0,2);
                if dims(l)>1
                    core0=reshape(core0,[currTT.r(1),nx*newR2]);
                    [u,s,v]=svd(core0,'econ');
                    evals=diag(s);
                    newR1=sum(evals>ep);
                    newR1=max(1,newR1);
                    u=u(:,1:newR1);
                    v=v(:,1:newR1);
                    s=s(1:newR1,1:newR1);
                    core0=v';
                    u=u*s;
                else
                    newR1=currTT.r(1);
                    % norm for dim=1
                    norm1(j)=norm(core0(:),'fro');
                    if (norm1(j) ~= 0)
                        core0=core0./norm1(j);
                    end
                end
                Gx_tt{j,l}.r(1)=newR1;
                ra(j,l)=newR1;
                Gx_tt{j,l}.r(2)=newR2;
                Gx_tt{j,l}.core=core0(:);
                Gx_tt{j,l}.ps(2)=newR1*newR2*nx+1;
            end
            if dims(1)>1
                labSend(u,labindex-1,j);
            end
        end
    end
    clear core0 currTT ru u v s
    nrm0=zeros(ne,1);
    for j=1:ne
        for l=1:numCo
            tab=nrm{l};
            nrm0(j)=nrm0(j)+sum(log(abs(tab(j,:))));
            var=ra{l};
            dimen=dims{l};
            for i=1:nd{l}
                ranks(k,j,dimen(i))=var(j,i);
            end
        end
        tab=norm1{1};
        nrm0(j)=nrm0(j)+log(abs(tab(j)));
        nrm0(j)=exp(nrm0(j)/d);
    end
    rat=nrm0(1)/nrm0(2);
    nrm0(2)=dx^(-1/2)*(1+rat^(2*d))^(-1/2/d);
    nrm0(1)=rat*nrm0(2);
    
    
    spmd
        for j=1:ne
            for i=1:nd
                Gx_tt{j,i}=Gx_tt{j,i}*nrm0(j);
            end
        end
        
        % Kinetic propagation
        for j=1:ne
            for i=1:nd
                Gp_tt=Sam_ttFT_split(Gx_tt{j,i},1);
                Gp_tt=Gp_tt.*KP_tt{i};
                Gx_tt{j,i}=Sam_ttFT_split(Gp_tt,2);
            end
        end
        Gp_tt=[];
        
        % Compute overlap with 2nd state
        if labindex==1
            p=1;
        else
            p=labReceive(labindex-1);
        end
        for i=1:nd
            p=dotCore(p,G0_tt{i},Gx_tt{2,i});
        end
        if labindex<numCo
            labSend(p,labindex+1);
        end
    end
    amp(k+1)=p{numCo}*dx^d;
    clear p
    
    save('overlap.mat','amp');
    if mod(k,5)==0
        save('ranks.mat','ranks');
    end
end
toc;
pool.delete;
exit