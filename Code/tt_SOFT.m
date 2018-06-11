% tt-SOFT propagation of non-adiabatic dynamics in pyrazine
% d-dimensional wavepacket on ne=2 electronic surfaces

% Uncomment this line if running locally
toolPath='~/Dropbox/Kenny/TT/newivan/TT-Toolbox';
% Uncomment this line if running on Omega
% toolPath='~/scratch/TT-Toolbox';

d=24; % number of dimensions
ns=620; % # of propagation steps
eps=7e-3; % precision for tt train round off

m = readparameters();
hbar = 1; % hbar
alpha=1; % initial width
x0= 0; % initial position
p0= 0; % initial momentum
nx=16;  % # of grid points in coord representation
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

% Generate rank-1 TT representation of initial wave function
G_tt=cell(1,ne);
G_tt{1}=tt_zeros(nx,d); % no population initially in S1 state
psi_tt=tt_tensor(exp(-alpha/2*(x-x0).^2+1i*p0*(x-x0)/hbar)*(alpha/pi)^0.25);
for i=1:d
    G_tt{2}=tkron(G_tt{2},psi_tt);
end
% Store the initial wave function to calculate autocorrelation function
% later
G0_tt=G_tt{2};

% Generate rank-1 TT representation of kinetic propagator
KP_tt=[];
for i=1:d
    KP_tt=tkron(KP_tt,tt_tensor(exp(-1i*tau/2*KE/m(i)/hbar)));
end

% Construct TT representations of potential operators
[v1m,v2m,vcm]=monV(x,d);
v1m=round(v1m,eps);
v2m=round(v2m,eps);
vcm=round(vcm,eps);
[v1d,v2d,vcd]=dimV(x,d);
v1d=round(v1d,eps);
v2d=round(v2d,eps);
vcd=round(vcd,eps);

%% This section performs the exponentiation of the matrix of tensor trains.
% The matrix is represented as a cell with tensor train elements

% vmat represents the initial matrix to be exponentiated
vmat=cell(ne,ne);
vmat{1,1}=-1i*tau*round(v1m+v1d,eps)/hbar;
vmat{2,2}=-1i*tau*round(v2m+v2d,eps)/hbar;
vmat{1,2}=-1i*tau*round(vcm+vcd,eps)/hbar;

% scale vmat such that it has no elements greater than 0.5
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

% Calculate the Taylor series for Y
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

% square y n0 times
for l=1:n0
    res=cell(ne,ne);
    for i=1:ne
        for k=1:ne
            res{i,k}=tt_zeros(nx,d);
            for j=1:ne
                res{i,k}=res{i,k}+y{i,j}.*y{j,k};
            end
            res{i,k}=round(res{i,k},eps);
        end
    end
    y=res;
end

% Store the diabatic terms in PP_tt{1} and PP_tt{2} and the adiabatic term
% in PA_tt
PP_tt=cell(1,ne);
PP_tt{1}=round(y{1,1},eps);
PP_tt{2}=round(y{2,2},eps);
PA_tt=round(y{1,2},eps);

%% Propagation

% Survival amplitude
amp=zeros(1,ns+1);
% as a check, this should be roughly 1 (i.e., time-0 autocorrelation function)
amp(1)=dot(G0_tt,G_tt{2})*dx^d;

disp('starting propagation');
tic;
for k=1:ns
    disp(k);
    % Kinetic propagation
    for j=1:ne
        Gp_tt=Sam_ttFT(G_tt{j},1,1); % Fourier transform to momentum representation
        Gp_tt=Gp_tt.*KP_tt;
        G_tt{j}=Sam_ttFT(Gp_tt,2,1); % Inverse fourier transform to position representation
    end
    
    % Potential propagation (multiply the S1 and S2 wave functions by the
    % matrix of tensor trains)
    GA_tt=cell(1,ne);
    for j=1:ne
        GA_tt{j}=round(G_tt{j}.*PA_tt,eps);
        G_tt{j}=round(G_tt{j}.*PP_tt{j},eps);
    end
    G_tt{1}=round(G_tt{1}+GA_tt{2},eps);
    G_tt{2}=round(G_tt{2}+GA_tt{1},eps);
    
    % Kinetic propagation
    for j=1:ne
        Gp_tt=Sam_ttFT(G_tt{j},1,1);
        Gp_tt=Gp_tt.*KP_tt;
        G_tt{j}=Sam_ttFT(Gp_tt,2,1);
    end
    
    % Calculate autocorrelation function
    amp(k+1)=dot(G0_tt,G_tt{2})*dx^d;
    
    save('overlap.mat','amp');
    if mod(k,5)==0
        save('ranks.mat','ranks');
    end
end
toc;
exit