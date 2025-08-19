function Hs_ustar_C35

clear; close all


%  10-m neutral wind array
U10N=[0.5:0.1:20];
% U10N=[7:1:10];
nU10=length(U10N);
% %  10-m neutral wind array
% U10N=[0.5:0.1:20];

%  U10 value for figure(9) wind+wave drift solution
%    (jUj_WD={2,4} for ms. Figs. 7,8)
% jU10j: U10N(jU10j)=[4 7 11 16 20] m/s
jU10j=[1+35 1+35+30 1+35+30+40 1+35+30+40+50 nU10];
% jUj_WD=1
jUj_WD=2
% jUj_WD=3
% jUj_WD=4
% jUj_WD=5
%
U10N_WD=U10N(jU10j(jUj_WD))

%  factor to reduce b0 by setting integration limit
%    b0 is not so sensitive to b0_fac 
b0_fac=1.0

% %  smallest wavenumber controls b0
% %
lambda_k0=0.15*U10N;
% lambda_k0=0.05*(1+0.001*floor(1000*(U10N/3).^2));
% lambda_k0=0.05*max(ones(1,nU10),0.001*floor(1000*(U10N/3).^2));
% whos lambda_k0
% lambda_k=[lambda_k0:0.05:50 51:0.25:100 101:2000];
% % lambda_k=[0.1:0.05:50 51:0.25:100 101:2000];
% kk=2*pi./lambda_k;

%%  SDRO
Omega=7.29*10^-5;
% % phi_lat=40;
% % phi_lat=20;
% phi_lat=10;
% fcor=2*Omega*sind(phi_lat)
% 
%  Set Coriolis parameter for re-scaling dimensionless solution
%   NH only; use complex conjugate for SH conversion
% phi_lat_array=[phi_lat];
% phi_lat_array=[40];
phi_lat_array=[0:5:90];
%
% phi_lat_array=[-90:2.5:-5 5:2.5:90];
% phi_lat_array=[-40 40];
nphi=length(phi_lat_array);


%  Phillips parameters
%  Phillips constant (1958)
alpha_P=0.0083;
%  equilibrium (1985)
gamma_beta3_I3pp1=10^-3
%
% %  Phillips 1985 p. 526; interpolating between two beta values
% %   for p=2 or so.
% betaIp=1.2*10^-2
%  Phillips 1985 p. 526; p=1/2
betaIp=3*10^-2
% %
% %  Value small enough to keep b0 < 1 for U10N = 20 m/s
% betaIp=0.7*10^-2



CDN0=1.3*10^-3;
ustar=sqrt(CDN0)*U10N;

% z0x=0.15*10^-3;
z0=zeros(1,nU10);

kappa=0.40;
gamma=0.11;
nu=1.5*10^-5;
g=9.81;
m_alpha= 0.0017;
b_alpha=-0.005;
alpha_c=min(m_alpha*U10N+b_alpha,0.028);
%     alpha_c=m_alpha*U10N+b_alpha;

rho_o=1025;
rho_a=1.25;
alpha_ao=sqrt(rho_a/rho_o);

rho0=rho_o;
mu=1.3*10^-2;
mu_MKS=10^-3*100*mu
nu_MKS=mu_MKS/rho0

%  figure variables
fig_a=0;
fig_b=0;
fsa=14;
fsa2=18;
msa=16;
Av_axis_limit=0;
z_grid_fac_j=[80 20 5 3];
% z_grid_fac_j=[20 20 5 2];


addpath ~/Documents/MATLAB/cbrewer
cmap128=cbrewer('div','RdBu',128);
cmap_b=cmap128(64:128,:);
cmap_bi=flipud(cmap_b);
cmap128=flipud(cmap128);
cmap_r=cmap128(64:128,:);

%%
%   H_s from u_star via COARE 3.5, derived from Edson et al 2013.

%  get COARE 3.5 z0

% first guess


for jU10=1:nU10

    U10Nx=U10N(jU10);
    alphax=alpha_c(jU10);
    ustara=ustar(jU10);
    ustar(jU10)=fzero(@get_ustar,ustara);
%         ustar(jU10)=fzero(@get_ustar,U10Nx);
    if(jU10<nU10)
        ustar(jU10+1)=ustar(jU10);
    end

    z0(jU10)=gamma*nu/ustar(jU10)+alphax*ustar(jU10)^2/g;

end

%  COARE 3.5 CDN
CDN=(kappa./log(10./z0)).^2;
%  COARE 3.5 z0_rough
z0_rough=max(0.,z0-gamma*nu./ustar);
% z0_rough=z0-gamma*nu./ustar;
%%   inferred Hs
Dcap=0.09;
ustar_cp=max(10^-9,0.03*(U10N-2)/8);
Hs=z0_rough./(Dcap*ustar_cp.^2);
% %
% ustar(end-20:end)
% 10^3*ustar_cp(end-20:end).^2
% 10^3*z0_rough(end-20:end)
% Hs(end-20:end)
% figure
% plot(U10N,10^3*z0,'-r','LineWidth',[1])
% hold on
% plot(U10N,ustar,'-g','LineWidth',[1])
% plot(U10N,10^2*ustar_cp,'-k','LineWidth',[1])
% plot(U10N,10^3*z0_rough,'-b','LineWidth',[1])
% plot(U10N,Hs,'-b','LineWidth',[2])
% axis([0 20 0 7])

%
cp_s0=ustar./ustar_cp;
k_s0=g./cp_s0.^2;
sigma_s0=sqrt(g*k_s0);
lambda_s0=2*pi./k_s0;
% %   inferred wavelength
% lambda_s0=2*pi*Hs*Dcap./alpha;
% %  wavenumber and frequency
% k_s0=2*pi./lambda_s0;
% sigma_s0=sqrt(g*k_s0);
% cp_s0=sigma_s0/k_s0;
%  Stokes wave mean
Us0=0.5*cp_s0.*(0.5*alpha_c/Dcap).^2;
% Us0=0.5*(sigma_s0./k_s0).*(0.5*alpha/Dcap).^2;
% Us0=0.5*(sigma_s0./k_s0).*(Hs*pi./lambda_s0).^2;
%
%   inferred wavelength
%    arbitrary factor of 1/2
lambda_s=0.5*lambda_s0;
%  wavenumber and frequency
k_s=2*pi./lambda_s;
sigma_s=sqrt(g*k_s);
cp_s=sigma_s/k_s;
%  Stokes wave mean
Us=0.5*cp_s.*(alpha_c/Dcap).^2;
% Us=0.5*(sigma_s./k_s).*(alpha_c/Dcap).^2;
% Us=0.5*(sigma_s./k_s).*(Hs*pi./lambda_s).^2;
% % steepnesses
% ka_s0=Hs*pi./lambda_s0;
% ka_s=Hs*pi./lambda_s;

%  Weber: Pierson etc
%  U_19.5: Weber (6.6)
C10=zeros(size(CDN));
C10(U10N < 15)=1.8*10^-3;
C10(nU10)=2.7*10^-3;
iC10a=nU10;
iC10b=1;
for iC10=1:nU10
    if(iC10a==nU10 && U10N(iC10)>=15)
        iC10a=iC10;
    end
end
iC10a;
iC10b=nU10-1;
for iC10=iC10a:iC10b
    C10(iC10)=C10(iC10a-1) ...
        +(C10(iC10b+1)-C10(iC10a-1))*(iC10-iC10a+1)/(iC10b+1-iC10a+1);
end
U19p5=U10N.*(1+sqrt(C10)*log(19.5/10)/kappa);
% U19p5=U10N.*(1+sqrt(CDN)*log(19.5/10)/kappa);
%  lambda: Weber (6.5)
lambda_W=(10^4/10^2)*2.803*10^-3*U19p5.^2;
% lambda_W=10^4*2.803*10^-3*U19p5.^2;
k_W=2*pi./lambda_W;
%  k_W a = 0.055 => a_W = 0.055/k_W = 0.055 x lambda_W/2pi
Hs_W=2*0.055./k_W;
% %  factor to make results agree with Weber ???
% Hs_W=3*Hs_W;
%  frequency and drift
sigma_W=sqrt(g*k_W);
cp_W=sigma_W./k_W;
Us_W=0.5*cp_W.*(0.055).^2;
% Us_W=0.5*(sigma_W./k_W).*(0.055).^2;
% Us_W=(sigma_W./k_W).*(Hs_W*pi./lambda_W).^2;
% Us_W=0.5*(sigma_W./k_W).*(Hs_W*pi./lambda_W).^2;
Us_W(end)

% Resio et al 1999
Ur=0.516*U10N.^(1.244);
Hs_fd=0.21*Ur.^2/g;
% Hs_fd=0.056*U10N.^(2.488)/g;
%  c_P=u_r
L_p=2*pi.*Ur.^2/g;
ustar_r=Ur/24.18;
z0_r=0.015*ustar_r.^2/g;
CDN_r=(kappa./log(10./z0_r)).^2;
% Us_W=(sigma_W./k_W).*(0.055).^2;
%  Stokes drift mean
k_r=2*pi./L_p;
sigma_r=sqrt(g*k_r);
cp_r=sigma_r./k_r;
Us_r=0.5*cp_r.*(Hs_fd*pi./L_p).^2;
% Us_r=0.5*(sigma_r./k_r).*(Hs_fd*pi./L_p).^2;


%  T_S estimates
% b0=0.5
b2=5
b0_bf=min(1,b2.*0.5*k_s0.*Hs);
b0_fd=min(1,b2.*0.5*k_r.*Hs_fd);
% bulk-flux
TS_bf=0.125*rho0*g*Hs.^2./(b0_bf.*sqrt(g./k_s0).*rho_a.*ustar.^2);
% fully-developed sea
TS_fd=0.125*rho0*g*Hs_fd.^2./(b0_fd.*sqrt(g./k_r).*rho_a.*ustar.^2);
%
P_a=b0_bf.*rho_a.*ustar.^2./(0.5*k_s0.*Hs);
% P_a=b0*rho_a*ustar.^2./(0.5*k_s0.*Hs);
P_a(Hs<=0)=NaN;
T_S0=rho0*sqrt(g)*0.5*Hs./(sqrt(k_s0).*P_a);
T_S0(Hs<=0)=NaN;
% fully-developed sea
P_a_fd=b0_fd.*rho_a.*ustar.^2./(0.5*k_r.*Hs_fd);
% P_a_fd=b0*rho_a*ustar.^2./(0.5*k_r.*Hs_fd);
P_a_fd(Hs_fd<=0)=NaN;
T_Sfd=rho0*sqrt(g)*0.5*Hs_fd./(sqrt(k_r).*P_a_fd);
T_Sfd(Hs_fd<=0)=NaN;



%  smallest wavenumber controls b0
%
lambda_k=[0.05:0.01:4.99 5:0.05:50 51:0.25:100 101:2000];
% lambda_k=[0.1:0.05:50 51:0.25:100 101:2000];
kk=2*pi./lambda_k;

% 
% % kk=(2*pi./2000)*[10000:-1:1];
% % lambda_k=2*pi./kk;
% lambda_k=[0.2:0.1:10 11:2000];
% kk=2*pi./lambda_k;

% jU10k=[16 28 44 64];
% U10N(jU10k)
% ustar(jU10k)

Tk_j=nan(length(kk),nU10);
Tk_j_plt=nan(length(kk),nU10);
% Tk_j_wt=nan(length(kk),nU10);

for jU10=1:nU10
    for jkk=1:length(kk)
        if(lambda_k(jkk) < b0_fac*L_p(jU10))
            Tk_j_plt(jkk,jU10)= 25*sqrt(g)*ustar(jU10)^-2*kk(jkk)^(-3/2);
            if(lambda_k(jkk) > lambda_k0(jU10))
              Tk_j(jkk,jU10)= 25*sqrt(g)*ustar(jU10)^-2*kk(jkk)^(-3/2);
            % % Tk_j(jkk,jU10)= alpha_P*g ...
            % %     /(2.*gamma_beta3_I3pp1*ustar(jU10)^3*kk(jkk)^2);
            % Tk_j_wt(jkk,jU10)=Tk_j(jkk,jU10) ...
            %     *sqrt(L_p(jU10)/(3*lambda_k(jkk)));
            end
        end
    end    
end

% Tk_j=zeros(length(kk),length(jU10k));
% 
% for jj=1:length(jU10k)
%     jU10=jU10k(jj);
%     Tk_j(:,jj)=(alpha_P*sqrt(g)/(gamma_beta3_I3pp1*ustar(jU10)^3))./sqrt(kk);
% end
% 

% zzk=[0:-0.25:-15];
% zzk=[0:-0.25:-20];
zzk=[0:-0.25:-25];
% zzk=[0:-0.25:-100];
% zzk(1)=-0.001;
zzk(1)=-0.01;
nzzk=length(zzk);

%  surface Us vs. latitude
Us_j_sfc=complex(zeros(nphi,nU10),zeros(nphi,nU10));

for jphi=1:nphi

    phi_lat=phi_lat_array(jphi)
    fcor=2*Omega*sind(phi_lat)

    Us_j=complex(zeros(nzzk,nU10),zeros(nzzk,nU10));
    GfUs_j=complex(zeros(nzzk,nU10),zeros(nzzk,nU10));
    Up_j=zeros(nzzk,nU10);
    Px_zint=zeros(1,nU10);
    b0_U10N=zeros(1,nU10);
    
    for jU10=1:nU10
        for jkk=2:length(kk)  
            %  0.5 factor to reduce b0
            if(  lambda_k(jkk) > lambda_k0(jU10) ...
              && lambda_k(jkk) < b0_fac*L_p(jU10))
                gammak_j=atan(fcor*Tk_j(jkk,jU10));
                Gammakf_j=1./(fcor*Tk_j(jkk,jU10));
                kk_j=0.5*(kk(jkk-1)+kk(jkk));
                for jz=1:nzzk
                    Us_j(jz,jU10)=Us_j(jz,jU10) ...
                        + ( cos(gammak_j)*exp(2*kk_j*zzk(jz)) ...
                            *exp(-i*gammak_j)/kk_j ) ...
                          *(kk(jkk-1)-kk(jkk));
                    GfUs_j(jz,jU10)=GfUs_j(jz,jU10) ...
                        + Gammakf_j*( cos(gammak_j)*exp(2*kk_j*zzk(jz)) ...
                            *exp(-i*gammak_j)/kk_j ) ...
                          *(kk(jkk-1)-kk(jkk));
                    Up_j(jz,jU10)=Up_j(jz,jU10) ...
                        + ( exp(2*kk(jkk)*zzk(jz))/kk_j ) ...
                          *(kk(jkk-1)-kk(jkk));
                end
                Px_zint(jU10)=Px_zint(jU10)+ ...
                         (kk(jkk-1)-kk(jkk))/(Tk_j(jkk,jU10)*kk_j^2);
           end
        end
    
        Us_j(:,jU10)=2*betaIp*ustar(jU10)*Us_j(:,jU10);
        GfUs_j(:,jU10)=2*betaIp*ustar(jU10)*GfUs_j(:,jU10);
        Up_j(:,jU10)=2*betaIp*ustar(jU10)*Up_j(:,jU10);
        Px_zint(jU10)=betaIp*ustar(jU10)*Px_zint(jU10);
    
        Us_j_sfc(jphi,jU10)=Us_j(1,jU10);

    end
    
    Us_j(Us_j==0)=NaN;
    Up_j(Up_j==0)=NaN;
end

Us_j_sfc(Us_j_sfc==0)=NaN;

% Us_j_sfc

% 
% b0_U10N=(rho_o/rho_a)*Px_zint./ustar.^2;
% 
% 
% rSTzeta=[0:0.01:5];
% ZZn_ap2_2=2.*(sinh(0.5*rSTzeta).^2)./sinh(rSTzeta);
% 
% ZZ_ap2_2=1.-(1./(2.*rSTzeta)) ...
%             .*( 4.*(1.-exp(-rSTzeta)) ...
%                 -(1.+ZZn_ap2_2).*(1.-exp(-2.*rSTzeta)) );
% %%

fsa=14


%%  plots

figure(1)
%
subplot(2,1,1)
%
pcolor(U10N,phi_lat_array,abs(Us_j_sfc))
shading interp
cbar=colorbar;
caxis([0 0.3])
cbar.Label.String='|U_S_0| (m s^-^1)';
cbar.Label.FontSize=fsa;
axis([0 20 0 90])
text(0.25,85,'a','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('NH latitude (^o)','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(2,1,2)
%
pcolor(U10N,phi_lat_array,(180/pi)*abs(angle(Us_j_sfc)))
% pcolor(U10N,phi_lat_array,(180/pi)*angle(Us_j_sfc))
shading interp
cbar=colorbar;
caxis([0 15])
% caxis([-15 0])
cbar.Label.String='Clockwise wind-relative angle(U_S_0) (^o)';
cbar.Label.FontSize=fsa;
axis([0 20 0 90])
text(0.25,85,'b','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('NH latitude (^o)','FontSize',fsa)
set(gca,'FontSize',fsa)


figure(2)
%
subplot(1,2,1)
%
pcolor(U10N,phi_lat_array,abs(Us_j_sfc))
shading interp
cbar=colorbar;
caxis([0 0.3])
cbar.Label.String='|U_S_0| (m s^-^1)';
cbar.Label.FontSize=fsa;
axis([0 20 0 90])
text(0.25,85,'a','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('NH latitude (^o)','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(1,2,2)
%
pcolor(U10N,phi_lat_array,(180/pi)*abs(angle(Us_j_sfc)))
% pcolor(U10N,phi_lat_array,(180/pi)*angle(Us_j_sfc))
shading interp
cbar=colorbar;
caxis([0 15])
% caxis([-15 0])
cbar.Label.String='Clockwise wind-relative angle(U_S_0) (^o)';
cbar.Label.FontSize=fsa;
axis([0 20 0 90])
text(0.25,85,'b','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('NH latitude (^o)','FontSize',fsa)
set(gca,'FontSize',fsa)


%%

%  nested function - COARE 3.5 iteration
function dustar = get_ustar(ustara)

    z0x=gamma*nu/ustara+alphax*ustara^2/g;
    CDN10sqrt=kappa/log(10/z0x);
    ustarx=CDN10sqrt*U10Nx;

    dustar=ustarx-ustara;

end




end
