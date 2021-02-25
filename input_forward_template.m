%% input
close all;
clear all;

% load image
vp=ones(4000,4000)*2000;

% dimensions
dt=10^-3;
dx=10;
dz=10;
nt=100;
nx=size(vp,1);
nz=size(vp,2);

X=[0,dx*nx];
Z=[0,-dz*nz];

% PML layers
lp=1;

% PML coefficient, usually 2
nPML=2;

% Theoretical coefficient, more PML layers, less R
% Empirical values
% lp=[10,20,30,40]
% R=[.1,.01,.001,.0001]
R=1;

% generate empty density
rho=ones(nx,nz)*1;
%rho=rho.*A;
% Lame constants for solid
mu=rho.*(vp/sqrt(3)).^2;
lambda=rho.*vp.^2;
% lambda(:,1:29)=0;

C.C11=lambda+2*mu;
C.C13=lambda;
C.C33=lambda+2*mu;
C.C55=mu;
C.rho=rho;
%%
% find surrounding layers
[a1,b1]=meshgrid(1:lp+2,1:nz);
[a2,b2]=meshgrid(nx-(lp+2)+1:nx,1:nz);
[a3,b3]=meshgrid(1:nx,1:lp+2);
[a4,b4]=meshgrid(1:nx,nz-(lp+2)+1:nz);
IND_air_layer=sub2ind([nx,nz],[reshape(a1,[1,(lp+2)*nz]),reshape(a2,[1,(lp+2)*nz]),reshape(a3,[1,(lp+2)*nx]),reshape(a4,[1,(lp+2)*nx])],...
    [reshape(b1,[1,(lp+2)*nz]),reshape(b2,[1,(lp+2)*nz]),reshape(b3,[1,(lp+2)*nx]),reshape(b4,[1,(lp+2)*nx])]);
clear a1 a2 a3 a4 b1 b2 b3 b4
% assign surrounding layers with air
%{
C11(IND_air_layer)=1145*340^2;
C13(IND_air_layer)=1145*340^2;
C15(IND_air_layer)=0;
C33(IND_air_layer)=1145*340^2;
C35(IND_air_layer)=0;
C55(IND_air_layer)=0;
rho(IND_air_layer)=1145;
%}


% source
% magnitude
M=2.7;

s_s1=[fix(nx/2)];
s_s3=ones(size(s_s1))*fix(nz/2);
%%
%{
figure;
set(gcf,'position',[80,80,1800,900]);
imagesc(X,Z,C.C55');
colorbar;
hold on
ax=plot(dx*s_s1,-dz*s_s3,'v','color',[1,0,0]);
set(gca,'ydir','normal');
hold on;
ax2=plot(dx*r1,-dz*r3,'^','color',[0,1,1]);
hold on;
colorbar;
xlabel('x [m]');
ylabel('z [m]');
title('C55 [Pa]');
legend([ax,ax2],...
    'source','receiver', ...
    'Location',[0.5,0.02,0.005,0.002],'orientation','horizontal');
print(gcf,['C55.png'],'-dpng','-r200');
%}
%%
% point interval in time steps
plot_interval=0;

p2=mfilename('fullpath');
if ~exist(p2,'dir')
    mkdir(p2)
end
tic;
for source_code=1:length(s_s1)
    % source locations
    s1=s_s1(source_code);
    s3=s_s3(source_code);
    
    % source frequency [Hz]
    freq=5;
    
    % source signal
    singles=rickerWave(freq,dt,nt,M);
    
    % give source signal to x direction
    src1=zeros(nt,1);
    src1=0*repmat(singles,[1,length(s3)]);
    
    % give source signal to z direction
    src3=src1;
    src3=1*repmat(singles,[1,length(s3)]);
    
    % receiver locations [m]
    %r1=[22:5:nx-22,22:5:nx-22,ones(size(22:5:nz-22))*25,ones(size(22:5:nz-22))*120];
    %r3=[ones(size(22:5:nx-22))*25,ones(size(22:5:nx-22))*80,22:5:nz-22,22:5:nz-22];
    r1=[4];
    r3=[4];
    
    s1t=dx*s1;
    s3t=-dz*s3;
    r1t=dx*r1;
    r3t=-dz*r3;
    
    % source type. 'D' for directional source. 'P' for P-source.
    source_type='D';
    
    
    % plot source
    plot_source=1;
    
    % figure path
    path=[p2 '/source_code_' num2str(source_code) '/'];
    
    % save wavefield
    save_wavefield=0;
    %% pass parameters to solver
    [v1,v3,R1,R3,P]=VTI_2D(dt,dx,dz,nt,nx,nz,...
        X,Z,...
        r1,r3,...
        s1,s3,src1,src3,source_type,...
        r1t,r3t, ...
        s1t,s3t, ...
        lp,nPML,R,...
        C, ...
        plot_interval,plot_source,...
        path,...
        save_wavefield);
    %{
    %% write to gif
    sources=[path 'pic/'];
    delaytime=.2;
    filename=['animation'];
    gifmaker(filename,delaytime,sources);
        %}
        %% write recordings
        %% recording saving
        if ~exist([path '/rec/'],'dir')
            mkdir([path '/rec/'])
        end
        
        DATA=rec_conversion(nt,nx,nz,dt,dx,dz,s1,s3,r1,r3);
        parsave([path '/rec/simu_info.mat'],DATA);
        
        DATA=R1;
        parsave([path '/rec/tR1.mat'],DATA);
        DATA=R3;
        parsave([path '/rec/tR3.mat'],DATA);
        DATA=P;
        parsave([path '/rec/tP.mat'],DATA);
        %% source saving
        %% recording saving
        if ~exist([path '/source/'],'dir')
            mkdir([path '/source/'])
        end
        DATA=src1;
        parsave([path '/source/src1.mat'],DATA);
        DATA=src3;
        parsave([path '/source/src3.mat'],DATA);
        DATA=source_type;
        parsave([path '/source/source_type.mat'],DATA);
end
toc;
%%