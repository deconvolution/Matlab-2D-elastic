function [v1,v3,R1,R3,P]=VTI_2D(dt,dx,dz,nt,nx,nz,...
    X,Z,...
    r1,r3,...
    s1,s3,src1,src3,source_type,...
    r1t,r3t, ...
    s1t,s3t, ...
    lp,nPML,R,...
    C, ...
    plot_interval,plot_source,...
    path,...
    save_wavefield)
%% create folder for figures

if ~exist(path,'dir')
    mkdir(path)
end

n_picture=1;
if save_wavefield==1
    if ~exist([path '/forward_wavefield/'],'dir')
        mkdir([path '/forward_wavefield/'])
    end
end

if plot_interval~=0
    if ~exist([path '/forward_pic/'],'dir')
        mkdir([path '/forward_pic/'])
    end
end

%% PML
vmax=sqrt((C.C33)./C.rho);
beta0=(ones(nx,nz,'single').*vmax*(nPML+1)*log(1/R)/2/lp/dx);
beta1=(zeros(nx,nz,'single'));
beta3=beta1;
tt=(1:lp)/lp;
tt2=repmat(reshape(tt,[lp,1]),[1,nz]);
plane_grad1=zeros(nx,nz);
plane_grad3=plane_grad1;

plane_grad1(2:lp+1,:)=flip(tt2,1);
plane_grad1(nx-lp:end-1,:)=tt2;
plane_grad1(1,:)=plane_grad1(2,:);
plane_grad1(end,:)=plane_grad1(end-1,:);

tt2=repmat(reshape(tt,[1,lp]),[nx,1]);
plane_grad3(:,2:lp+1)=flip(tt2,2);
plane_grad3(:,nz-lp:end-1)=tt2;
plane_grad3(:,1)=plane_grad3(:,2);
plane_grad3(:,end)=plane_grad3(:,end-1);

beta1=beta0.*plane_grad1.^nPML;
beta3=beta0.*plane_grad3.^nPML;

IND=unique(find(beta1.*beta3~=0));
beta=beta1+beta3;
beta(IND)=beta(IND)/2;
clear beta1 beta3 plane_grad1 plane grad3 vmax
%% source
ind_sor=sub2ind([nx,nz],s1,s3);
%% model configuration
R1=(zeros(nt,length(r1),'single'));
R3=(zeros(nt,length(r3),'single'));
P=(zeros(nt,length(r3),'single'));
ind_rec=sub2ind([nx,nz,3],r1,r3,ones(1,length(r1))*3);

%% monoclinic 2D solver xz plane (symmetric plane)
v1=(zeros(nx,nz,3,'single'));
v3=v1;
lim2=(zeros(4,2,'single'));

sigmas11=(zeros(nx,nz,'single'));
sigmas13=sigmas11;
sigmas33=sigmas11;
p=sigmas11;
ts=sigmas11;
ts2=sigmas11;

l=1;

switch source_type
    case 'D'
        v1(ind_sor)=v1(ind_sor)+.5./C.rho(ind_sor).*src1(l,:);
        v3(ind_sor)=v3(ind_sor)+.5./C.rho(ind_sor).*src3(l,:);
    case 'P'
        p(ind_sor)=p(ind_sor)+.5./C.rho(ind_sor).*src3(l,:);
end
%% save wavefield
if save_wavefield==1
    data=v1(:,:,2);
    save([path '/forward_wavefield/v1_' num2str(l) '.mat'],'data');
    data=v3(:,:,2);
    save([path '/forward_wavefield/v3_' num2str(l) '.mat'],'data');
    data=sigmas11;
    save([path '/forward_wavefield/sigma11_' num2str(l) '.mat'],'data');
    data=sigmas33;
    save([path '/forward_wavefield/sigma33_' num2str(l) '.mat'],'data');
    data=sigmas13;
    save([path '/forward_wavefield/sigma13_' num2str(l) '.mat'],'data');
    data=p;
    save([path '/forward_wavefield/p_' num2str(l) '.mat'],'data');
    
    data=v1(:,:,3);
    save([path '/forward_wavefield/v1_' num2str(l+1) '.mat'],'data');
    data=v3(:,:,3);
    save([path '/forward_wavefield/v3_' num2str(l+1) '.mat'],'data');
    data=sigmas11;
    save([path '/forward_wavefield/sigma11_' num2str(l+1) '.mat'],'data');
    data=sigmas33;
    save([path '/forward_wavefield/sigma33_' num2str(l+1) '.mat'],'data');
    data=sigmas13;
    save([path '/forward_wavefield/sigma13_' num2str(l+1) '.mat'],'data');
    data=p;
    save([path '/forward_wavefield/p_' num2str(l+1) '.mat'],'data');
end
%%
tic;
for l=2:nt-1
    for l2=1:2
        v1(:,:,l2)=v1(:,:,l2+1);
        v3(:,:,l2)=v3(:,:,l2+1);
    end
    %% compute sigma
    v1_x=D(v1(:,:,2),1)/dx;
    v3_x=D(v3(:,:,2),-1)/dx;
    v1_z=D(v1(:,:,2),-3)/dz;
    v3_z=D(v3(:,:,2),3)/dz;  
    
    sigmas11=dt*.5*((C.C11-C.C13).*v1_x...
        +(C.C13-C.C33).*v3_z)...
        +sigmas11...
        -dt*beta.*sigmas11;
    
    sigmas33=dt*.5*((C.C33-C.C13).*v3_z...
        +(C.C13-C.C11).*v1_x)...
        +sigmas33...
        -dt*beta.*sigmas33;
    
    sigmas13=dt*(C.C55.*(v1_z+v3_x))...
        +sigmas13...
        -dt*beta.*sigmas13;
    %% p
    p=-dt*((C.C11+C.C33)*.5.*v1_x+(C.C13+C.C33)*.5.*v3_z) ...
        +p ...
        -dt*beta.*p;
    %% compute v
    v1(:,:,3)=dt./C.rho.*(D(sigmas11-p,-1)/dx...
        +D(sigmas13,3)/dz)...
        +v1(:,:,3)...
        -dt*beta.*v1(:,:,2);
    v3(:,:,3)=dt./C.rho.*(D(sigmas13,1)/dx...
        +D(sigmas33-p,-3)/dz)...
        +v3(:,:,3)...
        -dt*beta.*v3(:,:,2);
    

    switch source_type
        case 'D'
            v1(ind_sor)=v1(ind_sor)+1./C.rho(ind_sor).*src1(l,:);
            v3(ind_sor)=v3(ind_sor)+1./C.rho(ind_sor).*src3(l,:);
        case 'P'
            p(ind_sor)=p(ind_sor)+1./C.rho(ind_sor).*src3(l,:);
    end
    %% fixed boundary condition
    %{
    v1(1,:,3)=0;
    v1(end,:,3)=0;
    v1(:,1,3)=0;
    v1(:,end,3)=0;
    
    v3(1,:,3)=0;
    v3(end,:,3)=0;
    v3(:,1,3)=0;
    v3(:,end,3)=0;
    %}
    %% assign recordings
    R1(l+1,:)=v1(ind_rec);
    R3(l+1,:)=v3(ind_rec);
    P(l+1,:)=p(ind_rec-nx*nz*2);
    %% save wavefield
    if save_wavefield==1
        
        data=v1(:,:,3);
        save([path '/forward_wavefield/v1_' num2str(l+1) '.mat'],'data');
        data=v3(:,:,3);
        save([path '/forward_wavefield/v3_' num2str(l+1) '.mat'],'data');
        data=sigmas11;
        save([path '/forward_wavefield/sigma11_' num2str(l+1) '.mat'],'data');
        data=sigmas33;
        save([path '/forward_wavefield/sigma33_' num2str(l+1) '.mat'],'data');
        data=sigmas13;
        save([path '/forward_wavefield/sigma13_' num2str(l+1) '.mat'],'data');
        data=p;
        save([path '/forward_wavefield/p_' num2str(l+1) '.mat'],'data');
    end
    %% plot
    if plot_interval~=0
        if mod(l,plot_interval)==0 || l==nt-1
            lim2(1,1)=min(min(v3,[],'all'),lim2(1,1));
            lim2(1,2)=max(max(v3,[],'all'),lim2(1,2));
            lim2(2,1)=min(min(sigmas33,[],'all'),lim2(2,1));
            lim2(2,2)=max(max(sigmas33,[],'all'),lim2(2,2));
            lim2(3,1)=min(min(sigmas13,[],'all'),lim2(3,1));
            lim2(3,2)=max(max(sigmas13,[],'all'),lim2(3,2));
            lim2(4,1)=min(min(p,[],'all'),lim2(4,1));
            lim2(4,2)=max(max(p,[],'all'),lim2(4,2));
            
            
            hfig=figure('Visible','off');
            set(gcf,'position',[80,80,1300,600]);
            subplot(2,3,1)
            imagesc(X,Z,v3(:,:,3)');
            set(gca,'ydir','normal');
            colorbar;
            xlabel({['x [m]']});
            ylabel({['z [m]']});
            title({['t=' num2str(l*dt) 's'],['v_3 [m/s]']});
            xlabel('x [m]');
            ylabel('z [m]');
            colorbar;
            hold on;
            if plot_source==1
                for i=1:length(s1)
                    ax2=plot(s1t(i),s3t(i),'v','color',[1,0,0]);
                    hold on;
                end
            end
            for i=1:length(r1)
                ax4=plot(r1t(i),r3t(i),'^','color',[0,1,1]);
                hold on;
            end
            ax3=plot(X(1)+[lp+1,lp+1]*dx,Z(1)-[lp+1,nz-lp-1]*dz,'color','blue');
            hold on;
            ax3=plot(X(1)+[nx-lp-1,nx-lp-1]*dx,Z(1)-[lp+1,nz-lp-1]*dz,'color','blue');
            hold on;
            hold on;
            ax3=plot(X(1)+[lp+1,nx-lp-1]*dx,Z(1)-[nz-lp-1,nz-lp-1]*dz,'color','blue');
            axis on;
            ax3=plot(X(1)+[lp+1,nx-lp-1]*dx,Z(1)-[lp+1,lp+1]*dz,'color','blue');
            axis on;
            
            subplot(2,3,2)
            imagesc(X,Z,sigmas33');
            set(gca,'ydir','normal');
            xlabel('x [m]');
            ylabel('z [m]');
            title('\sigma_{s33} [Pa]');
            colorbar;
            
            
            subplot(2,3,3)
            imagesc(X,Z,sigmas13');
            set(gca,'ydir','normal');
            xlabel('x [m]');
            ylabel('z [m]');
            title('\sigma_{s13} [Pa]');
            colorbar;
            
            subplot(2,3,4)
            imagesc(X,Z,p');
            set(gca,'ydir','normal');
            xlabel('x [m]');
            ylabel('z [m]');
            title('p [Pa]');
            colorbar;
            
            subplot(2,3,5)
            imagesc([1,length(r1)],[1,(l+1)]*dt,R3(1:l+1,:));
            colorbar;
            xlabel('Nr');
            ylabel('t [s]');
            title('R3 [m/s]');
            ylim([1,nt]*dt);
            
            subplot(2,3,6)
            imagesc(X,Z,C.C33');
            set(gca,'ydir','normal');
            xlabel('x [m]');
            ylabel('z [m]');
            title('C33 [Pa]');
            colorbar;
            hold on;
            if plot_source==1
                for i=1:length(s1)
                    ax2=plot(s1t(i),s3t(i),'v','color',[1,0,0]);
                    hold on;
                end
            end
            for i=1:length(r1)
                ax4=plot(r1t(i),r3t(i),'^','color',[0,1,1]);
                hold on;
            end
            ax3=plot(X(1)+[lp+1,lp+1]*dx,Z(1)-[lp+1,nz-lp-1]*dz,'color','blue');
            hold on;
            ax3=plot(X(1)+[nx-lp-1,nx-lp-1]*dx,Z(1)-[lp+1,nz-lp-1]*dz,'color','blue');
            hold on;
            hold on;
            ax3=plot(X(1)+[lp+1,nx-lp-1]*dx,Z(1)-[nz-lp-1,nz-lp-1]*dz,'color','blue');
            axis on;
            ax3=plot(X(1)+[lp+1,nx-lp-1]*dx,Z(1)-[lp+1,lp+1]*dz,'color','blue');
            axis on;
            if plot_source==1
                legend([ax2,ax3,ax4],...
                    'source','PML boundary','receiver',...
                    'Location',[0.5,0.02,0.005,0.002],'orientation','horizontal');
            else
                legend([ax3,ax4],...
                    'PML boundary','receiver',...
                    'Location',[0.5,0.02,0.005,0.002],'orientation','horizontal');
            end
            
            print(gcf,[path './forward_pic/' num2str(n_picture) '.png'],'-dpng','-r200');
            n_picture=n_picture+1;
            
        end
    end
    
    fprintf('\n time step=%d/%d',l+1,nt);
    fprintf('\n    epalsed time=%.2fs',toc);
    fprintf('\n    n_picture=%d',n_picture);
    d=clock;
    fprintf('\n    current time=%d %d %d %d %d %.0d',d(1),d(2),d(3),d(4),d(5),d(6));
    
end
%%
