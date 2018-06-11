% 4-mode calculation for pyrazine
% comparing tt-SOFT with exact SOFT
% Sam Greene

clear all
addpath(genpath('~/Dropbox/Kenny/TT/newivan/TT-Toolbox'));
ns=680; % # of propagation steps
m = readparameters();
m=m([1,2,3,6]);
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
d=4; % number of dimensions
Nexp=10; % number of terms in the exponential taylor series

% coordinate and momentum grids
x = ((1:nx)-nx/2)*dx;
psi = exp(-alpha/2*(x-x0).^2+1i*p0*(x-x0))*(alpha/pi)^0.25;
p=fftshift(dp*(-np/2:np/2-1));
KE=p.^2/2;

% KE propagator
KP_tt=[];
for i=1:d
    prop=tt_tensor(exp(-1i*tau/2*KE/m(i)/hbar));
    KP_tt=tkron(KP_tt,prop);
end
KP=full(KP_tt,np*ones(1,d));

% Wave function
G_tt=cell(1,ne);
G_tt{1}=tt_zeros(nx,d);
psi_tt=tt_tensor(psi);
G_tt{2}=psi_tt;
for i=2:d
    G_tt{2}=tkron(G_tt{2},psi_tt);
end
G0_tt=G_tt{2};
psiX=cell(1,ne);
psiX{1}=zeros(nx*ones(1,d));
psiX{2}=full(G_tt{2},nx*ones(1,d));
psi0=psiX{2};

% Potential
[v1m,v2m,vcm]=monV4D(x,d);
[v1d,v2d,vcd]=dimV4D(x,d);

PE=cell(1,ne);
PE{1}=zeros(nx*ones(1,d));
PE{2}=zeros(nx*ones(1,d));
PC=zeros(nx*ones(1,d));

PP=cell(1,ne);
PP{1}=zeros(nx*ones(1,1));
PP{2}=zeros(nx*ones(1,d));
PA=zeros(nx*ones(1,d));

tic;
q=zeros(1,24);
for i=1:nx
    q(1)=x(i);
    for j=1:nx
        q(2)=x(j);
        for k=1:nx
            q(3)=x(k);
            for l=1:nx
                q(6)=x(l);
                vm=pot_mono3(q,d);
                vd=pot_dim3(q);
                
                [V,D]=eig(vm+vd);
                D=diag(D);
                expMat=exp(-1i*tau*D/hbar);
                expMat=diag(expMat);
                prop=(V*expMat*V');
                PP{1}(i,j,k,l)=prop(1,1);
                PP{2}(i,j,k,l)=prop(2,2);
                PA(i,j,k,l)=prop(1,2);
            end
        end
    end
end

% Full-rank SOFT propagation
xi=zeros(1,ns+1);
xi(1)=1;
disp('Full-rank');
for k=1:ns
    % Kinetic propagation
    for j=1:ne
        psiP=fftn(psiX{j})* (dx/sqrt(2*pi*hbar))^d;
        psiP=psiP.*KP;
        psiX{j}=ifftn(psiP)*(np*dp/sqrt(2*pi*hbar))^d;
    end
    % Potential propagation
    GA=cell(1,ne);
    for j=1:ne
        GA{j}=psiX{j}.*PA;
        psiX{j}=psiX{j}.*PP{j};
    end
    psiX{1}=psiX{1}+GA{2};
    psiX{2}=psiX{2}+GA{1};
    % Kinetic propagation
    for j=1:ne
        psiP=fftn(psiX{j})* (dx/sqrt(2*pi*hbar))^d;
        psiP=psiP.*KP;
        psiX{j}=ifftn(psiP)*(np*dp/sqrt(2*pi*hbar))^d;
    end
    
    xi(k+1)=dx^d*sum(reshape(psi0.*psiX{2},[1,nx^d]));
end
toc;

% TT propagation
epsVal=[0.02,0.01];
Neps=size(epsVal,2);
xi_tt=zeros(Neps,ns+1);
spec_tt=zeros(Neps,2,2^11);
ra=zeros(Neps,ne,ns);
for a=1:Neps
    tic;
    eps=epsVal(a);
    G_tt{1}=tt_zeros(nx,d);
    G_tt{2}=G0_tt;
    xi_tt(a,1)=dot(G0_tt,G_tt{2})*dx^d;
    
    % vmat represents the matrix to be exponentiated
    vmat=cell(2,2);
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
                res{i,k}=round(res{i,k}/l,eps/10);
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
                res{i,k}=round(res{i,k},eps/10*(0.5^(n0-l)));
            end
        end
        y=res;
    end
    
    PP_tt{1}=y{1,1};
    PP_tt{2}=y{2,2};
    PA_tt=y{1,2};
    
    for k=1:ns
        % Kinetic propagation
        for j=1:ne
            FT_G_tt=Sam_ttFT(G_tt{j},1,1);
            FT_G_tt=FT_G_tt.*KP_tt;
            G_tt{j}=Sam_ttFT(FT_G_tt,2,1);
        end
        % Potential propagation
        GA_tt=cell(1,ne);
        for j=1:ne
            GA_tt{j}=G_tt{j}.*PA_tt;
            GA_tt{j}=round(GA_tt{j},eps);
            G_tt{j}=G_tt{j}.*PP_tt{j};
            G_tt{j}=round(G_tt{j},eps);
        end
        
        G_tt{1}=round(G_tt{1}+GA_tt{2},eps);
        G_tt{2}=round(G_tt{2}+GA_tt{1},eps);
        
        n1old=(dot(G_tt{1},G_tt{1}))^(1/2/d);
        n2old=(dot(G_tt{2},G_tt{2}))^(1/2/d);
        r=n1old/n2old;
        n2=dx^(-1/2)*(1+r^(2*d))^(-1/2/d);
        n1=r*n2;
        G_tt{1}=G_tt{1}*(n1/n1old)^d;
        G_tt{2}=G_tt{2}*(n2/n2old)^d;
        
        % Kinetic propagation
        for j=1:ne
            FT_G_tt=Sam_ttFT(G_tt{j},1,1);
            FT_G_tt=FT_G_tt.*KP_tt;
            G_tt{j}=Sam_ttFT(FT_G_tt,2,1);
        end
        
        xi_tt(a,k+1)=dot(G0_tt,G_tt{2})*dx^d;
        ra(a,1,k)=G_tt{1}.ps(d+1);
        ra(a,2,k)=G_tt{2}.ps(d+1);
    end
    toc;
    spec_tt(a,:,:)=Plotspec(squeeze(xi_tt(a,:)),tau,d);
    spec_tt(a,1,:)=1240./spec_tt(a,1,:);
    spec_tt(a,2,:)=spec_tt(a,2,:)/max(spec_tt(a,2,:));
end

%% Plotting ranks
dtt=tau*2.418e-2;
figSize=[300,200];
fig=figure('Position',[100,400,figSize]);
box on
hold on
plot(dtt*(1:ns),squeeze(ra(1,1,:)),'Color',[0.5 0.5 0.5],'LineWidth',1);
plot(dtt*(1:ns),squeeze(ra(1,2,:)),'--','Color',[0.5 0.5 0.5],'LineWidth',1);
plot(dtt*(1:ns),squeeze(ra(2,1,:)),'-k','LineWidth',1);
plot(dtt*(1:ns),squeeze(ra(2,2,:)),'--k','LineWidth',1);
xlabel('Time (fs)');
ylabel('Tensor-Train Size');
legend([char(949),' = ',num2str(epsVal(1)),', S_1'],[char(949),' = ',num2str(epsVal(1)),', S_2'],[char(949),' = ',num2str(epsVal(2)),', S_1'],[char(949),' = ',num2str(epsVal(2)),', S_2'],'Location','northwest');
axis([0 ns*dtt 0 6500]);
hold off
saveas(fig,'~/Dropbox/Kenny/TT/ttSOFT/Paper/ranks4D.eps','epsc');

%% Plotting spectra
labSize=[15,20];
fig=figure('Position',[100,100,figSize]);
spec=Plotspec(xi,tau,d);
spec(1,:)=1240./spec(1,:); % convert to nm
spec(2,:)=spec(2,:)/max(spec(2,:)); % normalize
plot(spec(1,:),spec(2,:),'-k','LineWidth',1);
axis([220 280 0 1.1]);
hold on
plot(squeeze(spec_tt(1,1,:)),squeeze(spec_tt(1,2,:)),'Color',[0.5 0.5 0.5],'LineWidth',1);
plot(squeeze(spec_tt(2,1,:)),squeeze(spec_tt(2,2,:)),'--k','LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Intensity');
legend('benchmark',[char(949),' = ',num2str(epsVal(1))],[char(949),' = ',num2str(epsVal(2))],'Location','northwest');
uicontrol('Style','text','Position',[0,figSize(2)-labSize(2),labSize],'String','(a)','BackgroundColor','white');
set(gca,'yticklabel',num2str(get(gca,'YTick')','%.1f'))
saveas(fig,'~/Dropbox/Kenny/TT/ttSOFT/Paper/spec.eps','epsc');

fig=figure('Position',[400,100,figSize]);
time=dtt*(0:ns);
plot(time,real(xi),'-k','LineWidth',1);
axis([0 165 -1 1]);
hold on
plot(time,squeeze(real(xi_tt(1,:))),'Color',[0.5 0.5 0.5],'LineWidth',1);
plot(time,squeeze(real(xi_tt(2,:))),'--k','LineWidth',2);
xlabel('Time (fs)');
ylabel('Re[C(t)]');
legend('benchmark',[char(949),' = ',num2str(epsVal(1))],[char(949),' = ',num2str(epsVal(2))],'Location','northeast');
uicontrol('Style','text','Position',[0,figSize(2)-labSize(2),labSize],'String','(b)','BackgroundColor','white');
set(gca,'yticklabel',num2str(get(gca,'YTick')','%.1f'))
saveas(fig,'~/Dropbox/Kenny/TT/ttSOFT/Paper/ReC.eps','epsc');

fig=figure('Position',[700,100,figSize]);
plot(time,imag(xi),'-k','LineWidth',1);
axis([0 165 -1 1]);
hold on
plot(time,squeeze(imag(xi_tt(1,:))),'Color',[0.5 0.5 0.5],'LineWidth',1);
plot(time,squeeze(imag(xi_tt(2,:))),'--k','LineWidth',2);
xlabel('Time (fs)');
ylabel('Im[C(t)]');
legend('benchmark',[char(949),' = ',num2str(epsVal(1))],[char(949),' = ',num2str(epsVal(2))],'Location','northeast');
uicontrol('Style','text','Position',[0,figSize(2)-labSize(2),labSize],'String','(c)','BackgroundColor','white');
set(gca,'yticklabel',num2str(get(gca,'YTick')','%.1f'))
saveas(fig,'~/Dropbox/Kenny/TT/ttSOFT/Paper/ImC.eps','epsc');
