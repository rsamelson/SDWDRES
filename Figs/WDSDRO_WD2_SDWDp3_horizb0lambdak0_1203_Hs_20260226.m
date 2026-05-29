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

% %  smallest wavelength controls b0
% %
lambda_k1=0.15*U10N;
% lambda_k1=0.05*(1+0.001*floor(1000*(U10N/3).^2));
% lambda_k1=0.05*max(ones(1,nU10),0.001*floor(1000*(U10N/3).^2));
% whos lambda_k1
% lambda_k=[lambda_k1:0.05:50 51:0.25:100 101:2000];
% % lambda_k=[0.1:0.05:50 51:0.25:100 101:2000];
% kk=2*pi./lambda_k;

%%  SDRO
Omega=7.29*10^-5;
phi_lat=40;
% phi_lat=14;
fcor=2*Omega*sind(phi_lat)

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

% z0x=0.152*10^-3;
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



%  smallest wavelength controls b0
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
            if(lambda_k(jkk) > lambda_k1(jU10))
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

Us_j=complex(zeros(nzzk,nU10),zeros(nzzk,nU10));
GfUs_j=complex(zeros(nzzk,nU10),zeros(nzzk,nU10));
Up_j=zeros(nzzk,nU10);
Px_zint=zeros(1,nU10);
b0_U10N=zeros(1,nU10);
Hs_Psi=zeros(1,nU10);

whos Us_j

for jU10=1:nU10
    for jkk=2:length(kk)  
        %  0.5 factor to reduce b0
        if(  lambda_k(jkk) > lambda_k1(jU10) ...
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
            Hs_Psi(jU10)=Hs_Psi(jU10)+ ...
                     (kk(jkk-1)-kk(jkk))/(kk_j^2.5);
            
       end
    end

    Us_j(:,jU10)=2*betaIp*ustar(jU10)*Us_j(:,jU10);
    GfUs_j(:,jU10)=2*betaIp*ustar(jU10)*GfUs_j(:,jU10);
    Up_j(:,jU10)=2*betaIp*ustar(jU10)*Up_j(:,jU10);
    Px_zint(jU10)=betaIp*ustar(jU10)*Px_zint(jU10);
    Hs_Psi(jU10)=4*sqrt(betaIp*ustar(jU10)*sqrt(1/g)*Hs_Psi(jU10));

end

Us_j(Us_j==0)=NaN;
Up_j(Up_j==0)=NaN;

b0_U10N=(rho_o/rho_a)*Px_zint./ustar.^2;


rSTzeta=[0:0.01:5];
ZZn_ap2_2=2.*(sinh(0.5*rSTzeta).^2)./sinh(rSTzeta);

ZZ_ap2_2=1.-(1./(2.*rSTzeta)) ...
            .*( 4.*(1.-exp(-rSTzeta)) ...
                -(1.+ZZn_ap2_2).*(1.-exp(-2.*rSTzeta)) );
% %%

fsa=14


%%  plots

figure(1)
%
subplot(3,3,1)
plot(U10N,rho_a*ustar.^2,'-b','LineWidth',[2])
% hold on
% plot(U10N,b0.*ustar.^2,'-b','LineWidth',[1])
% plot(U10N,0.5*k_s0.*Hs.*P_a,'-b')
% hold on
% plot(U10N,P_a_fd,'-b','LineWidth',[2])
text(1,1.1,'a','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('\tau_0 (N m^-^2)','FontSize',fsa)
% ylabel(strcat('P_a (N m^-^2), b_0=',num2str(b0)),'FontSize',fsa)
set(gca,'FontSize',fsa)
%
%
subplot(3,3,2)
plot(U10N,rho_a*b0_bf.*ustar.^2,'-b','LineWidth',[1])
hold on
plot(U10N,rho_a*b0_fd.*ustar.^2,'-b','LineWidth',[2])
% plot(U10N,0.5*k_s0.*Hs.*P_a,'-b')
% hold on
% plot(U10N,P_a_fd,'-b','LineWidth',[2])
text(1,0.9,'b','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('b_0 \tau_0 (N m^-^2)','FontSize',fsa)
% ylabel(strcat('P_a (N m^-^2), b_0=',num2str(b0)),'FontSize',fsa)
set(gca,'FontSize',fsa)
%
%
subplot(3,3,3)
plot(U10N,b0_bf,'-b','LineWidth',[1])
hold on
plot(U10N,b0_fd,'-b','LineWidth',[2])
% plot(U10N,0.5*k_s0.*Hs.*P_a,'-b')
% hold on
% plot(U10N,P_a_fd,'-b','LineWidth',[2])
axis([0 20 0 0.8])
text(1,0.72,'c','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('b_2ka','FontSize',fsa)
% ylabel(strcat('b_0 = b_2ka; b_2=',num2str(b2)),'FontSize',fsa)
% ylabel(strcat('P_a (N m^-^2), b_0=',num2str(b0)),'FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(3,3,4)
plot(U10N,0.5*Hs,'-b')
hold on
% plot(U10N,Hs_W,'--r')
plot(U10N,0.5*Hs_fd,'-b','LineWidth',[2])
axis([0 20 0 5])
text(1,4.5,'d','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('a=H_s/2 (m)','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(3,3,5)
plot(U10N,lambda_s0,'-b')
hold on
plot(U10N,L_p,'-b','LineWidth',[2])
text(1,270,'e','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('\lambda (m)','FontSize',fsa)
axis([0 20 0 300])
set(gca,'FontSize',fsa)
%
subplot(3,3,6)
plot(U10N,k_s0.*0.5.*Hs,'-b','LineWidth',[1])
hold on
plot(U10N,k_r.*0.5.*Hs_fd,'-b','LineWidth',[2])
% plot(U10N,10*0.5*(k_s0.*0.5.*Hs).^2,'--b','LineWidth',[1])
% hold on
% plot(U10N,10*0.5*(k_r.*0.5.*Hs_fd).^2,'--b','LineWidth',[2])
%
% plot(U10N,0.5*k_s0.*Hs.*P_a,'-b')
% hold on
% plot(U10N,P_a_fd,'-b','LineWidth',[2])
text(1,0.135,'f','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('ka','FontSize',fsa)
% ylabel('ka, 10 \times k^2a^2/2','FontSize',fsa)
% ylabel(strcat('P_a (N m^-^2), b_0=',num2str(b0)),'FontSize',fsa)
set(gca,'FontSize',fsa)
%
%
subplot(3,3,7)
plot(U10N,k_s0,'-b')
hold on
plot(U10N,k_r,'-b','LineWidth',[2])
text(1,0.9,'g','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('k=2\pi/\lambda (m^-^1)','FontSize',fsa)
axis([0 20 0 1])
set(gca,'FontSize',fsa)
%
subplot(3,3,8)
plot(U10N,TS_bf/(3600),'-b')
hold on
plot(U10N,TS_fd/(3600),'-b','LineWidth',[2])
% plot(U10N,T_S0/(3600),'-b')
% hold on
% plot(U10N,T_Sfd/(3600),'-b','LineWidth',[2])
axis([0 20 0 8])
text(1,7.2,'h','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('T_e_q (hr)','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(3,3,9)
plot(U10N,(180/pi)*atan(fcor*TS_bf),'-b')
hold on
plot(U10N,(180/pi)*atan(fcor*TS_fd),'-b','LineWidth',[2])
% plot(U10N,(180/pi)*atan(fcor*T_S0),'-b')
% hold on
% plot(U10N,(180/pi)*atan(fcor*T_Sfd),'-b','LineWidth',[2])
axis([0 20 0 70])
text(1,063,'i','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('arctan(fT_e_q) (^o)','FontSize',fsa)
set(gca,'FontSize',fsa)


figure
%
subplot(3,2,1)
%
pcolor(U10N,kk,log10(Tk_j_plt/3600))
shading interp
colorbar
caxis([-2 1.5])
axis([0 20 0 2*pi])
text(0.5,0.9*2*pi,'a','FontSize',fsa)
ylabel('k (m^-^1)','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
% title('log_1_0[T_k (hr)]','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(3,2,2+1)
%
pcolor(U10N,lambda_k,log10(Tk_j_plt/3600))
shading interp
colorbar
caxis([-2 1.5])
axis([0 20 0 200])
% set('XDir','reverse')
text(0.5,180,'c','FontSize',fsa)
ylabel('\lambda_k (m)','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
% title('log_1_0[T_k (hr)]','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(3,2,2*2+1)
%
pcolor(U10N,lambda_k/(4*pi),log10(Tk_j_plt/3600))
shading interp
colorbar
caxis([-2 1.5])
axis([0 20 0 16])
% set('XDir','reverse')
text(0.5,14.4,'e','FontSize',fsa)
ylabel('\delta_s (m)','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
% title('log_1_0[T_k (hr)]','FontSize',fsa)
set(gca,'FontSize',fsa)

% figure
%
subplot(3,2,2)
%
pcolor(U10N,kk,(180/pi)*atan(fcor*Tk_j_plt))
shading interp
colorbar
caxis(70*[0 1])
axis([0 20 0 2*pi])
text(0.5,0.9*2*pi,'b','FontSize',fsa)
ylabel('k (m^-^1)','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
% title('Stokes-drift angle (^o) at 40^o N ','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(3,2,2*2)
%
pcolor(U10N,lambda_k,(180/pi)*atan(fcor*Tk_j_plt))
shading interp
colorbar
caxis(70*[0 1])
axis([0 20 0 200])
% set('XDir','reverse')
text(0.5,180,'d','FontSize',fsa)
ylabel('\lambda_k (m)','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
% title('Stokes-drift angle (^o) at 40^o N','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(3,2,3*2)
%
pcolor(U10N,lambda_k/(4*pi),(180/pi)*atan(fcor*Tk_j_plt))
shading interp
colorbar
caxis(70*[0 1])
axis([0 20 0 16])
% set('XDir','reverse')
text(0.5,14.4,'f','FontSize',fsa)
ylabel('\delta_s (m)','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
% title('Stokes-drift angle (^o) at 40^o N','FontSize',fsa)
set(gca,'FontSize',fsa)


figure

colormap(cmap_r)

subplot(1,3,1)
%
pcolor(U10N,zzk,abs(Us_j))
shading interp
axis([0 20 -10 0])
colorbar
caxis([0 0.20])
text(0.5,-0.5,'a','FontSize',fsa)
ylabel('z (m)','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
title('|U_S_d| at 40^o N','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(1,3,2)
%
pcolor(U10N,zzk,-(180/pi)*angle(Us_j))
shading interp
axis([0 20 -10 0])
text(0.5,-0.5,'b','FontSize',fsa)
ylabel('z (m)','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
title('U_S_d angle (^o) at 40^o N','FontSize',fsa)
set(gca,'FontSize',fsa)
colorbar
%
subplot(1,3,3)
%
pcolor(U10N,zzk,Up_j)
shading interp
axis([0 20 -10 0])
colorbar
caxis([0 0.20])
text(0.5,-0.5,'c','FontSize',fsa)
ylabel('z (m)','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
title('U_S at 40^o N','FontSize',fsa)
set(gca,'FontSize',fsa)


figure
%
% %  10-m neutral wind array
% U10N=[0.5:0.1:20];
jU10j=[1+35 1+35+30 1+35+30+40 1+35+30+40+50 nU10];
jpanelj=['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j'];
xpanelj=[0.05 0.1 0.2 0.25 0.4];
%
for jUj=1:length(jU10j)

    subplot(2,5,jUj)
    plot(abs(Us_j(:,jU10j(jUj))),zzk,'-b','Linewidth',2)
    hold on
    plot(Up_j(:,jU10j(jUj)),zzk,'--b','Linewidth',1)    
    axis([0 xpanelj(jUj) -15 0])
    text(0.9*xpanelj(jUj),-14,jpanelj(jUj),'FontSize',fsa)
    ylabel('z (m)','FontSize',fsa)
    xlabel('|U_S_d|','FontSize',fsa)
    % xlabel('|U_S_d| at 40^o N','FontSize',fsa)
    plt_name=strcat('U10N=',num2str(U10N(jU10j(jUj))),'m/s')
    title(plt_name,'FontSize',fsa)
    set(gca,'FontSize',fsa)
    %
    %  Breivik - Phillips spectrum profile function
    % k_r=omega_P^2/g;
    k_p=2*pi/L_p(jU10j(jUj));
    omega_p=sqrt(g*k_p);
    ubar_P=(2*alpha_P*g/omega_p) ...
            .*(exp(2*k_p.*zzk) ...
                  -sqrt(-2*pi*k_p.*zzk).*erfc(sqrt(-2*k_p.*zzk)));
    % hold on
    plot(ubar_P,zzk,'-g','Linewidth',1)

    %
    subplot(2,5,5+jUj)
    plot(-(180/pi)*angle(Us_j(:,jU10j(jUj))),zzk,'-b','Linewidth',2)
    axis([0 50 -15 0])
    text(4,-14,jpanelj(5+jUj),'FontSize',fsa)
    ylabel('z (m)','FontSize',fsa)
    xlabel('U_S_d angle (^o)','FontSize',fsa)
    % xlabel('U_S_d angle (^o) at 40^o N','FontSize',fsa)
    set(gca,'FontSize',fsa)

end



figure

plot(U10N,-(180/pi)*angle(Us_j(1,:)),'-b','Linewidth',2)
axis([0 20 0 12])
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('U_S_d(z=0) angle (^o)','FontSize',fsa)
title('Surface dynamic Stokes drift angle','FontSize',fsa)
% title('Surface Stokes angle right of downwind in Northern Hemisphere','FontSize',fsa)
set(gca,'FontSize',fsa)

figure
%
plot(2.*rSTzeta,ZZ_ap2_2,'-b','LineWidth',2)
hold on
plot(2.*rSTzeta,ZZn_ap2_2,'--b','LineWidth',1)
xlabel('\Gamma_e_q T_\zeta = 2 r_S T_\zeta','FontSize',fsa)
ylabel('2|ZZ|/a_p^2 (-), 2|ZZ(t_n)|/a_p^2 (--)','FontSize',fsa)
set(gca,'FontSize',fsa)

% 
% figure
% %
% subplot(2,1,1)
% %
% plot(U10N,b0_U10N,'-b','LineWidth',2)
% axis([0 20 0 1])
% text(1,0.9,'a','FontSize',fsa)
% xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
% ylabel('b_0','FontSize',fsa)
% set(gca,'FontSize',fsa)
% %
% subplot(2,1,2)
% %
% plot(U10N,lambda_k1,'-b','LineWidth',2)
% axis([0 20 0 3])
% text(1,2.7,'b','FontSize',fsa)
% xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
% ylabel('\lambda_k_0 (m)','FontSize',fsa)
% set(gca,'FontSize',fsa)


figure
%
subplot(1,3,1)
%
plot(U10N,b0_U10N,'-b','LineWidth',2)
axis([0 20 0 1])
text(1,0.9,'a','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('b_0','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(1,3,2)
%
plot(U10N,lambda_k1,'-b','LineWidth',2)
axis([0 20 0 3])
text(1,2.7,'b','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('\lambda_k_0 (m)','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(1,3,3)
%
plot(U10N,Hs_Psi,'-b','LineWidth',2)
hold on
plot(U10N,Hs,'--k','LineWidth',1)
plot(U10N,Hs_fd,'--g','LineWidth',1)
axis([0 20 0 10])
text(1,9,'c','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('H_s (m)','FontSize',fsa)
set(gca,'FontSize',fsa)

figure
%
subplot(1,2,1)
%
plot(U10N,b0_U10N,'-b','LineWidth',2)
axis([0 20 0 1])
text(1,0.9,'a','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('b_0','FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(1,2,2)
%
plot(U10N,lambda_k1,'-b','LineWidth',2)
axis([0 20 0 3])
text(1,2.7,'b','FontSize',fsa)
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('\lambda_k_1 (m)','FontSize',fsa)
% ylabel('\lambda_k_0 (m)','FontSize',fsa)
set(gca,'FontSize',fsa)
% 
% 
% figure(9)
% 
% U10N(jU10j(jUj_WD))
% 
% for jplot=1:4
%     %
%     subplot(2,4,jplot)
%     %
%     plot(abs(Us_j(:,jU10j(jUj_WD))),log10(abs(zzk)),'-b','Linewidth',2)
%     hold on
%     plot(real(Us_j(:,jU10j(jUj_WD))),log10(abs(zzk)),'-b','Linewidth',1)
%     plot(imag(Us_j(:,jU10j(jUj_WD))),log10(abs(zzk)),'--b','Linewidth',1)
%     % axis([max(abs(Us_j(:,jU10j(jUj_WD))))*[-1 1] -3 2]);
%     % plot(Up_j(:,jU10j(jUj_WD)),log10(abs(zzk)),'--b','Linewidth',1)
%     % ylabel('z (m)','FontSize',fsa)
%     % xlabel('|U_S_d|,u_p,u_B at 40^o N','FontSize',fsa)
%     % plt_name=strcat('U10N=',num2str(U10N(jU10j(jUj_WD))),'m/s')
%     % title(plt_name,'FontSize',fsa)
%     % set(gca,'FontSize',fsa)
%     %
%     %  Breivik - Phillips spectrum profile function
%     % k_r=omega_P^2/g;
%     k_p=2*pi/L_p(jU10j(jUj_WD));
%     omega_p=sqrt(g*k_p);
%     ubar_P=(2*alpha_P*g/omega_p) ...
%             .*(exp(2*k_p.*zzk) ...
%                   -sqrt(-2*pi*k_p.*zzk).*erfc(sqrt(-2*k_p.*zzk)));
%     hold on
%     plot(ubar_P,log10(abs(zzk)),':r','Linewidth',1)
%     % axis([max(abs(ubar_P))*[-1 1] -3 2]);
%     %
%     subplot(2,4,4+jplot)
%     %
%     plot(abs(Us_j(:,jU10j(jUj_WD))),zzk,'-b','Linewidth',2)
%     hold on
%     plot(real(Us_j(:,jU10j(jUj_WD))),zzk,'-b','Linewidth',1)
%     plot(imag(Us_j(:,jU10j(jUj_WD))),zzk,'--b','Linewidth',1)
%     % axis([max(abs(Us_j(:,jU10j(jUj_WD))))*[-1 1] -15 0]);
%     % plot(Up_j(:,jU10j(jUj_WD)),zzk,'--b','Linewidth',1)
%     % ylabel('z (m)','FontSize',fsa)
%     % xlabel('|U_S_d|,u_p,u_B at 40^o N','FontSize',fsa)
%     % plt_name=strcat('U10N=',num2str(U10N(jU10j(jUj_WD))),'m/s')
%     % title(plt_name,'FontSize',fsa)
%     % set(gca,'FontSize',fsa)
%     %
%     %  Breivik - Phillips spectrum profile function
%     % k_r=omega_P^2/g;
%     k_p=2*pi/L_p(jU10j(jUj_WD));
%     omega_p=sqrt(g*k_p);
%     ubar_P=(2*alpha_P*g/omega_p) ...
%             .*(exp(2*k_p.*zzk) ...
%                   -sqrt(-2*pi*k_p.*zzk).*erfc(sqrt(-2*k_p.*zzk)));
%     hold on
%     plot(ubar_P,zzk,':r','Linewidth',1)
%     % axis([max(abs(ubar_P))*[-1 1] -15 0]);
% 
% end

            % ustar_o=ustar_o_U10N(1)
            % tau0_o=rho_o*ustar_o^2
            % % ustar_o=ustar_o_U10N(jU10)
            % %  remove pressure stress from wave-drift model
            % b0tau0=0.01
            % b0_ustar_o=sqrt(b0tau0/rho_o)
            % ustar_o=ustar_o-b0_ustar_o
            % 
            % 
tau0_U10N=rho_a*ustar(jU10j(jUj_WD))^2
tauw_U10N=rho_o*Px_zint(jU10j(jUj_WD))
b0=tauw_U10N/tau0_U10N


%%  WDHES

%
%  Set physical and COARE3.5 constants
%
Omega=7.29*10^-5;
%
kappa=0.40;
gamma=0.11;
nu=1.5*10^-5;
g=9.81;
%
rho_o=1025;
rho_a=1.25;
alpha_ao=sqrt(rho_a/rho_o);
%
rho0=rho_o;
% Ocean viscosities for viscous sublayer
% value in g cm^-1 s^-1 from Phillips (1977)
mu=1.3*10^-2;
mu_MKS=10^-3*100*mu;
nu_MKS=mu_MKS/rho0;


%  Set 10-m neutral wind array for profiles

%  10-m neutral wind array
% U10N_max=7
% U10N_max=20
% U10N_max=30
% U10N_max=50
% U10N=[7];
U10N=U10N_WD
% U10N=[0.25:0.25:U10N_max];
nU10=length(U10N)

%  get COARE 3.5 neutral z0 and friction velocities
%
m_alpha= 0.0017;
b_alpha=-0.005;
alpha_c=m_alpha*U10N+b_alpha;
% alpha_c=min(m_alpha*U10N+b_alpha,0.028);
% first guess array for air friction velocity
%
CDN0=1.3*10^-3;
ustar_a=sqrt(CDN0)*U10N;
%
% z0x=0.152*10^-3;
z0=zeros(1,nU10);
%
%  loop over U10N and use get_ustar function to solve neutral COARE 3.5
%  
for jU10=1:nU10    
    U10Nx=U10N(jU10);
    alphax=alpha_c(jU10);
    ustara=ustar_a(jU10);
    ustar_a(jU10)=fzero(@get_ustar,ustara);
    if(jU10<nU10)
        ustar_a(jU10+1)=ustar_a(jU10);
    end    
    z0(jU10)=gamma*nu/ustar_a(jU10)+alphax*ustar_a(jU10)^2/g;    
end
%
%  COARE 3.5 CDN
CDN=(kappa./log(10./z0)).^2;
%  COARE 3.5 z0_rough
z0_rough=z0-gamma*nu./ustar_a;
%  COARE 3.5 water friction velocity
ustar_o_U10N=alpha_ao*ustar_a;

% %  index for cut-off at U10N = 20 m/s
% % U10N_20=[0.25:0.25:20];
% U10N_cutoff=20;
% jU10N_20=length([0.25:0.25:U10N_cutoff]);

%  Set Coriolis parameter for re-scaling dimensionless solution
%   NH only; use complex conjugate for SH conversion
phi_lat_array=[40];
% phi_lat_array=[5:2.5:90];
%
% phi_lat_array=[-90:2.5:-5 5:2.5:90];
% phi_lat_array=[-40 40];
nphi=length(phi_lat_array);

for j_mono_spec=1:4

    %  log-Ekman layer parameters
    if(j_mono_spec==1)
        % param 1
        z0o_U10N=10.^(-5+7*(1-exp(-U10N/10)));
        phi_U10N=exp(-U10N/10);
        %
        lab_ac='a';
        lab_bd='b';
    elseif(j_mono_spec==2)
        % param 2
        z0o_U10N=10.^(-6+7.5*(1-exp(-U10N/10)));
        phi_U10N=exp(-U10N/15);
        %
        lab_ac='c';
        lab_bd='d';
   elseif(j_mono_spec==3)
        % param 3
        z0o_U10N=10.^(-5+7*(1-exp(-U10N/5)));
        phi_U10N=0.05+0.95*exp(-U10N/5);
        %
        lab_ac='a';
        lab_bd='b';
    elseif(j_mono_spec==4)
        % param 4
        z0o_U10N=10.^(-6+7.5*(1-exp(-U10N/5)));
        phi_U10N=0.05+0.95*exp(-U10N/7.5);
        %
        lab_ac='c';
        lab_bd='d';
    end
    % 
    % U10N_WD=zeros(nU10+1,nphi);
    % U0_WD=zeros(nU10+1,nphi);
    % U_DSm=zeros(nU10+1,nphi);
    % U_DSz=zeros(nU10+1,nphi);
    % U4_WD=zeros(nU10+1,nphi);
    % U6_WD=zeros(nU10+1,nphi);
    % U8_WD=zeros(nU10+1,nphi);
    % U_Drm=zeros(nU10+1,nphi);
    % U_Drz=zeros(nU10+1,nphi);
    % U_Mim=zeros(nU10+1,nphi);
    % U_Miz=zeros(nU10+1,nphi);


    for jphi=1:nphi

        phi_lat=phi_lat_array(jphi)
        % %
        % % For Southern Hemisphere, set flag for complex conjugate profiles
        % if(phi_lat < 0)
        %     SH_flag=1;
        %     phi_lat=abs(phi_lat)
        % else
        %     SH_flag=0;
        % end
        fcor=2*Omega*sind(phi_lat)
        
        
        
        %%  Generate WDHES model wind-drift profiles
        %
        %
        %  loop over mono vs. spec parameterizations
        %    (WDHES JPO 2022, eqs. (50)-(51) and (52)-(53)
        %
    
        % %  For U10N > 20 m/s, set z0o and phi_w to values for U10N = 20 m/s
        % z0o_U10N(U10N>20)=z0o_U10N(jU10N_20);
        % phi_U10N(U10N>20)=phi_U10N(jU10N_20);
    
        %  arrays for profile plots and output
        nzprof=5000;
        UprofU10N=complex(nan(nzprof,nU10),nan(nzprof,nU10));
        U10profU10N=nan(nzprof,nU10);
        zprofU10N=nan(nzprof,nU10);
        % epsprofU10N=nan(nzprof,nU10);
        
        %%
        
        
        %%  Solve wind-drift model
        %    Eqs.  (12), (31), (32), (36)-(38) with
        %     j_mono_spec=1 => monochromatic calibration parameters:  (50)-(51)
        %     j_mono_spec=2 => spectral calibration parameters:  (52)-(53)
       
        %  U10N loop
        for jU10=1:nU10
            %
            phi=phi_U10N(jU10);
            z0_o_dim=z0o_U10N(jU10);
            ustar_o=ustar_o_U10N(jU10)
            % % ---
            % %  Do not
            % %  remove pressure stress from wave-drift model
            % %  scaling here.  Instead, adjust forcing amplitude only
            % %  in equation for Ah.
            % ustar_o=sqrt(1-b0)*ustar_o
            % % ---
            %
            D_o=ustar_o/fcor;
            %
 
            % b0tau0=0.01
            % b0_ustar_o=sqrt(b0tau0/rho_o)
            % ustar_o=ustar_o-b0_ustar_o
        %     ustar_o=alpha_ao*ustar_a(jU10);
        
            %%  log-Ekman layer
            % %  dimensionless depth z/D, D=ustar_o/f
        
            pkappa=phi*kappa;
            
            %  check for viscous sublayer
            z0_o_dim=max(z0_o_dim,nu_MKS/(pkappa*ustar_o));
        
            % ocean - dimensionless
            z0_o=z0_o_dim*fcor/ustar_o;
            %
            %
        
            %  Set profile of dimensionless eddy viscosity Av/Av0,
            %      where Av0=ustar_o*D=ustar_o^2/f
        
            Av_max=0.03;
            Av_min=0.1*Av_max;
            z1=z0_o-Av_max/pkappa;
            z2=2*z1;
            pkappa2=pkappa;
            z3=z2-(Av_max-Av_min)/pkappa2;
             %     % z1 z2 z3
            zza=z1*(exp(0.01*[0:1000])-1)/(exp(10)-1);
            xia=(1+i)*sqrt(2*(z0_o-zza)/pkappa);
            dz1=0.001*(z1-z2);
            zz1=[z1-dz1:-dz1:z2];
            dz2=0.001*(z2-z3);
            zz2=[z2-dz2:-dz2:z3];
            xi2=(1+i)*sqrt(2*(z0_o-z1-z2+zz2)/pkappa2);
            dz3=0.01*(z3+20);
            zz3=[z3-dz3:-dz3:-20];
            zz=[zza zz1 zz2 zz3];

    
%  wave-breaking forcing
%    particular solutions Ap*K0+Bp*I0
            %
            xi_K0=(1+i).*sqrt(2.*(z0_o-zz)./pkappa);
            %  integrate to get Ap, Bp
            nzprof1=min(length(zz),nzprof);
            nzz=nzprof1
            zprof_K0=D_o*zz(1:nzprof1);
            zprof_K0(end)
            zzk(end)
            zprof_K0(1)
            zzk(1)
            % nzz=length(zz);
            Ap=complex(zeros(1,nzz),zeros(1,nzz));
            Bp=complex(zeros(1,nzz),zeros(1,nzz));
            GfWs=complex(zeros(1,nzz),zeros(1,nzz));
            zzk_int=[zprof_K0(1) zzk zzk(end)-1 zprof_K0(end)-1];
            zzk_int(end-2:end)
            zzk_int(1:3)
            GfUs_int=[GfUs_j(1,jU10j(jUj_WD)) ...
                      GfUs_j(:,jU10j(jUj_WD))' 0+0*i 0+0*i];
            GfUs_int(end-2:end)
            GfUs_int(1:3)
            GfWs=interp1(zzk_int,GfUs_int,zprof_K0);
            % GfWs=interp1(zzk,GfUs_j(:,jU10j(jUj_WD)),zprof_K0);
            whos GfWs
            max(abs(GfWs))
            GfWs(1:10)
            for jz=nzz-1:-1:1
                Ap(jz)=Ap(jz+1) ...
                  + 0.5*(xi_K0(jz)*besseli(0,xi_K0(jz))*GfWs(jz) ...
                        +xi_K0(jz+1)*besseli(0,xi_K0(jz+1))*GfWs(jz+1))...
                      *(xi_K0(jz)-xi_K0(jz+1));
                Bp(jz)=Bp(jz+1) ...
                  + 0.5*(xi_K0(jz)*besselk(0,xi_K0(jz))*GfWs(jz) ...
                        +xi_K0(jz+1)*besselk(0,xi_K0(jz+1))*GfWs(jz+1))...
                      *(xi_K0(jz)-xi_K0(jz+1));
            end
            whos Ap
            max(abs(Ap))
            Ap=-i*Ap;
            Bp=i*Bp;
            Wp=Ap.*besselk(0,xi_K0(1:nzprof1)) ...
              +Bp.*besseli(0,xi_K0(1:nzprof1));

            whos Wp
            max(abs(Wp))
            min(abs(Wp))



            xi0=(1+i)*sqrt(2*z0_o/pkappa);
            Kcap0z0=besselk(0,xi0);
            Kcap1z0=besselk(1,xi0);
            Icap1z0=besseli(1,xi0);

            %  Remove pressure-stress fraction here.
            Ah=(exp(-i*0.25*pi) ...
                *sqrt(ustar_o^2*D_o/(pkappa*z0_o_dim))*(1-b0) ...
                -Ap(1)*Kcap1z0+Bp(1)*Icap1z0)/Kcap1z0;

            %  wind-drift profile
            Wcap=Ah*besselk(0,xi_K0(1:nzprof1))+Wp;

            %  Us wave-drift profile
            Ws=complex(zeros(1,nzz),zeros(1,nzz));
            % zzk_int=[zprof_K0(1) zzk zzk(end)-1 zprof_K0(end)-1];
            % zzk_int(end-2:end)
            % zzk_int(1:3)
            Us_int=[Us_j(1,jU10j(jUj_WD)) ...
                      transpose(Us_j(:,jU10j(jUj_WD))) 0+0*i 0+0*i];
            Ws=interp1(zzk_int,Us_int,zprof_K0);

            %  wind-drift + wave-drift
            Wtot=Wcap+Ws;

            % %  K0 only - no wave-breaking forcing
            %  Get K0 solution for z1,z2,z3 -> -\infty
            xi0=(1+i)*sqrt(2*z0_o/pkappa);
            Kcap0z0=besselk(0,xi0);
            Kcap1z0=besselk(1,xi0);
            A0_K0=1./(-pkappa*z0_o*Kcap1z0*(-2*i/(pkappa*xi0)));
            % A0_K0=0.5*A0_K0;
            xi_K0=(1+i).*sqrt(2.*(z0_o-zz)./pkappa);
            Bcap_K0=A0_K0.*besselk(0,xi_K0);
            % %  Pressure-stress has not been removed from scaling
            % %   so do not add it again here.
            % % Uprof_K0=(ustar_o/(1-b0))*Bcap_K0(1:nzprof1);
            Uprof_K0=ustar_o*Bcap_K0(1:nzprof1);

            % % A0_K0=1./(-pkappa*z0_o*Kcap1z0*(-2*i/(pkappa*xi0)));
            % % A0_K0=(1-b0)*A0_K0;
            % % Bcap_K0=A0_K0.*besselk(0,xi_K0);
            % % 
            % % whos Bcap_K0
            % 
            % %  get dimensional ocean velocity profiles
            % %
            % % D_o=ustar_o/fcor;
            % %
            % % nzprof1=min(length(zz),nzprof);
            % Uprof_K0=ustar_o*Bcap_K0(1:nzprof1);
            % zprof_K0=D_o*zz(1:nzprof1);
            log10_zprof_K0=log10(-zprof_K0);

            fsa=14;
            %
            % % U10Nj = 7 m/s
            % jU10p=7*4;
            %
            % U10Nj = max(U10N)
            jU10p=nU10;
            %
            U10N(jU10p);



            % figure(9)
            % 
            % %
            % subplot(2,4,j_mono_spec)
            % plot(abs(Wcap),log10_zprof_K0, ...
            %     '-g','LineWidth',2)
            % hold on
            % plot(real(Wcap),log10_zprof_K0, ...
            %     '-g','LineWidth',1)
            % plot(imag(Wcap),log10_zprof_K0, ...
            %     '--g','LineWidth',1)
            % plot([0 0],[-6 4],':k')
            % axis([max([abs(Wcap) ...
            %            abs(Us_j(1,jU10j(jUj_WD))) ...
            %                ubar_P(1)])*[-1 1] -3 2]);
            % xlabel('Velocity (m s^-^1)','FontSize',fsa)
            % ylabel('log_1_0[|z| (m)]','FontSize',fsa)
            % set(gca,'YDir','reverse')
            % pltnam=strcat('|U|,(U,V) for U_1_0_N =',num2str(U10N(jU10p)), ...
            %     ' m s^-^1',': p',num2str(j_mono_spec));
            % % pltnam=strcat('|U|,(U,V) for U_1_0_N = 7 m s^-^1', ...
            % %    ': Monochromatic calibration parameters');
            % title(pltnam,'FontSize',fsa)
            % set(gca,'FontSize',fsa)
            % % 
            %  %
            % subplot(2,4,4+j_mono_spec)
            % %
            % plot(abs(Wcap),zprof_K0, ...
            %     '-g','LineWidth',2)
            % hold on
            % plot(real(Wcap),zprof_K0, ...
            %     '-g','LineWidth',1)
            % plot(imag(Wcap),zprof_K0, ...
            %     '--g','LineWidth',1)
            % plot([0 0],[-15 0],':k')
            % axis([max([abs(Wcap) ...
            %            abs(Us_j(1,jU10j(jUj_WD))) ...
            %                ubar_P(1)])*[-1 1] -15 0]);
            % xlabel('Velocity (m s^-^1)','FontSize',fsa)
            % ylabel('z (m)]','FontSize',fsa)
            % % set(gca,'YDir','reverse')
            % pltnam=strcat('|U|,(U,V) for U_1_0_N =',num2str(U10N(jU10p)), ...
            %     ' m s^-^1',': p',num2str(j_mono_spec));
            % % pltnam=strcat('|U|,(U,V) for U_1_0_N = 7 m s^-^1', ...
            % %    ': Monochromatic calibration parameters');
            % title(pltnam,'FontSize',fsa)
            % set(gca,'FontSize',fsa)
            % % 
            % 
            % figure(10)
            % %
            % subplot(2,4,j_mono_spec)
            % %
            % plot(abs(Wp),log10_zprof_K0,'-b','LineWidth',2);
            % hold on
            % plot(real(Wp),log10_zprof_K0,'-b','LineWidth',1);
            % plot(imag(Wp),log10_zprof_K0,'--b','LineWidth',1);
            % plot([0 0],[-6 4],':k')
            % axis([max(abs(Wp))*[-1 1] -3 2]);
            % xlabel('W_p (m s^-^1)','FontSize',fsa)
            % ylabel('log_1_0[|z| (m)]','FontSize',fsa)
            % set(gca,'YDir','reverse')
            % pltnam=strcat('W_p for U_1_0_N =',num2str(U10N(jU10p)), ...
            %     ' m s^-^1',': p',num2str(j_mono_spec));
            % title(pltnam,'FontSize',fsa)
            % set(gca,'FontSize',fsa)
            % %
            % subplot(2,4,4+j_mono_spec)
            % %
            % plot(abs(Wp),zprof_K0,'-b','LineWidth',2);
            % hold on
            % plot(real(Wp),zprof_K0,'-b','LineWidth',1);
            % plot(imag(Wp),zprof_K0,'--b','LineWidth',1);
            % plot([0 0],[-15 0],':k')
            % axis([max(abs(Wp))*[-1 1] -15 0]);
            % xlabel('W_p (m s^-^1)','FontSize',fsa)
            % ylabel('z (m)','FontSize',fsa)
            % % set(gca,'YDir','reverse')
            % pltnam=strcat('W_p for U_1_0_N =',num2str(U10N(jU10p)), ...
            %     ' m s^-^1',': p',num2str(j_mono_spec));
            % title(pltnam,'FontSize',fsa)
            % set(gca,'FontSize',fsa)
            % 
            figure(9)
            %
            subplot(2,4,j_mono_spec)
            %
            plot(abs(Uprof_K0),log10_zprof_K0,'-g','LineWidth',2);
            hold on
            % plot(real(Uprof_K0),log10_zprof_K0,'-g','LineWidth',1);
            plot(imag(Uprof_K0),log10_zprof_K0,'--g','LineWidth',1);
            plot(abs(Wtot),log10_zprof_K0,'-b','LineWidth',2);
            % plot(real(Wtot),log10_zprof_K0,'-b','LineWidth',1);
            plot(imag(Wtot),log10_zprof_K0,'--b','LineWidth',1);
            plot([0 0],[-6 4],':k')
            axis([max([abs(Wtot) abs(Uprof_K0)])*[-1 1] -3 2]);
            xlabel('W_p (m s^-^1)','FontSize',fsa)
            ylabel('log_1_0[|z| (m)]','FontSize',fsa)
            set(gca,'YDir','reverse')
            pltnam=strcat('W for U_1_0_N =',num2str(U10N(jU10p)), ...
                ' m s^-^1',': p',num2str(j_mono_spec));
            title(pltnam,'FontSize',fsa)
            set(gca,'FontSize',fsa)
            %
            subplot(2,4,4+j_mono_spec)
            %
            plot(abs(Uprof_K0),zprof_K0,'-g','LineWidth',2);
            hold on
            % plot(real(Uprof_K0),zprof_K0,'-g','LineWidth',1);
            plot(imag(Uprof_K0),zprof_K0,'--g','LineWidth',1);
            plot(abs(Wtot),zprof_K0,'-b','LineWidth',2);
            % plot(real(Wtot),zprof_K0,'-b','LineWidth',1);
            plot(imag(Wtot),zprof_K0,'--b','LineWidth',1);
            plot([0 0],[-15 0],':k')
            axis([max(abs([abs(Wtot) abs(Uprof_K0)]))*[-1 1] -100 0]);
            xlabel('W_p (m s^-^1)','FontSize',fsa)
            ylabel('z (m)','FontSize',fsa)
            % set(gca,'YDir','reverse')
            pltnam=strcat('W for U_1_0_N =',num2str(U10N(jU10p)), ...
                ' m s^-^1',': p',num2str(j_mono_spec));
            title(pltnam,'FontSize',fsa)
            set(gca,'FontSize',fsa)

            Wtot(log10_zprof_K0<-2.1)=NaN;
            Wcap(log10_zprof_K0<-2.1)=NaN;
            Ws(log10_zprof_K0<-2.1)=NaN;
            Ws(zprof_K0<-25)=NaN;
            jms_lab=['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'];
            deg_min=[-120 -120 -90 -90 -90];

            if(mod(j_mono_spec,2)==1)

                figure(10)
                %
                subplot(2,4,j_mono_spec)
                %
                plot(abs(Uprof_K0),log10_zprof_K0,'-g','LineWidth',2);
                hold on
                plot(abs(Wcap),log10_zprof_K0,'-r','LineWidth',2);
                plot(abs(Ws),log10_zprof_K0,'--r','LineWidth',2);
                plot(abs(Wtot),log10_zprof_K0,'-b','LineWidth',2);
                % plot(real(Wtot),log10_zprof_K0,'-b','LineWidth',1);
                % plot(imag(Wtot),log10_zprof_K0,'--b','LineWidth',1);
                % plot(abs(Uprof_K0),log10_zprof_K0,'-g','LineWidth',2);
                % plot(real(Uprof_K0),log10_zprof_K0,'-g','LineWidth',1);
                % plot(imag(Uprof_K0),log10_zprof_K0,'--g','LineWidth',1);
                plot([0 0],[-6 4],':k')
                axis([max([abs(Wtot) abs(Uprof_K0)])*[0 1] -2 1.5]);
                text(0.9*max(abs(Wtot)),1.2, ...
                    jms_lab(j_mono_spec),'FontSize',fsa)
                xlabel('|U_m|,|W|,|U_S_d| (m s^-^1)','FontSize',fsa)
                ylabel('log_1_0[|z| (m)]','FontSize',fsa)
                set(gca,'YDir','reverse')
                pltnam=strcat('U_1_0_N =',num2str(U10N(jU10p)), ...
                    ' m s^-^1',': p',num2str(j_mono_spec));
                title(pltnam,'FontSize',fsa)
                set(gca,'FontSize',fsa)
                %
                subplot(2,4,4+j_mono_spec)
                %
                plot(abs(Uprof_K0),zprof_K0,'-g','LineWidth',2);
                hold on
                plot(abs(Wcap),zprof_K0,'-r','LineWidth',2);
                plot(abs(Ws),zprof_K0,'--r','LineWidth',2);
                plot(abs(Wtot),zprof_K0,'-b','LineWidth',2);
                % plot(real(Wtot),zprof_K0,'-b','LineWidth',1);
                % plot(imag(Wtot),zprof_K0,'--b','LineWidth',1);
                % plot(abs(Uprof_K0),zprof_K0,'-g','LineWidth',2);
                % plot(real(Uprof_K0),zprof_K0,'-g','LineWidth',1);
                % plot(imag(Uprof_K0),zprof_K0,'--g','LineWidth',1);
                plot([0 0],[-15 0],':k')
                axis([max(abs([abs(Wtot) abs(Uprof_K0)]))*[0 1] -15 0]);
                text(0.9*max(abs(Wtot)),-13, ...
                    jms_lab(4+j_mono_spec),'FontSize',fsa)
                xlabel('|U_m|,|W|,|U_S_d| (m s^-^1)','FontSize',fsa)
                ylabel('z (m)','FontSize',fsa)
                % set(gca,'YDir','reverse')
                pltnam=strcat('U_1_0_N =',num2str(U10N(jU10p)), ...
                    ' m s^-^1',': p',num2str(j_mono_spec));
                title(pltnam,'FontSize',fsa)
                set(gca,'FontSize',fsa)
        %
                subplot(2,4,j_mono_spec+1)
                %
                plot((180/pi)*abs(angle(Uprof_K0)),log10_zprof_K0,'-g','LineWidth',2);
                hold on
                plot((180/pi)*abs(angle(Wcap)),log10_zprof_K0,'-r','LineWidth',2);
                plot((180/pi)*abs(angle(Ws)),log10_zprof_K0,'--r','LineWidth',2);
                plot((180/pi)*abs(angle(Wtot)),log10_zprof_K0,'-b','LineWidth',2);
                % plot(real(Wtot),log10_zprof_K0,'-b','LineWidth',1);
                % plot(imag(Wtot),log10_zprof_K0,'--b','LineWidth',1);
                % plot((180/pi)*angle(Uprof_K0),log10_zprof_K0,'-g','LineWidth',2);
                % plot(real(Uprof_K0),log10_zprof_K0,'-g','LineWidth',1);
                % plot(imag(Uprof_K0),log10_zprof_K0,'--g','LineWidth',1);
                % plot([0 0],[-6 4],':k')
                axis([0 abs(deg_min(jUj_WD)) -2 1.5]);
                text(0.1*abs(deg_min(jUj_WD)),1.2, ...
                    jms_lab(j_mono_spec+1),'FontSize',fsa)
                xlabel('angle(U_m,W,U_S_d) (^o)','FontSize',fsa)
                ylabel('log_1_0[|z| (m)]','FontSize',fsa)
                set(gca,'YDir','reverse')
                pltnam=strcat('U_1_0_N =',num2str(U10N(jU10p)), ...
                    ' m s^-^1',': p',num2str(j_mono_spec));
                title(pltnam,'FontSize',fsa)
                set(gca,'FontSize',fsa)
                %
                subplot(2,4,4+j_mono_spec+1)
                %
                plot((180/pi)*abs(angle(Uprof_K0)),zprof_K0,'-g','LineWidth',2);
                hold on
                plot((180/pi)*abs(angle(Wcap)),zprof_K0,'-r','LineWidth',2);
                plot((180/pi)*abs(angle(Ws)),zprof_K0,'--r','LineWidth',2);
                plot((180/pi)*abs(angle(Wtot)),zprof_K0,'-b','LineWidth',2);
                % plot(real(Wtot),zprof_K0,'-b','LineWidth',1);
                % plot(imag(Wtot),zprof_K0,'--b','LineWidth',1);
                % plot((180/pi)*angle(Uprof_K0),zprof_K0,'-g','LineWidth',2);
                % plot(real(Uprof_K0),zprof_K0,'-g','LineWidth',1);
                % plot(imag(Uprof_K0),zprof_K0,'--g','LineWidth',1);
                plot([0 0],[-15 0],':k')
                axis([0 abs(deg_min(jUj_WD)) -15 0]);
                text(0.1*abs(deg_min(jUj_WD)),-13, ...
                    jms_lab(4+j_mono_spec+1),'FontSize',fsa)
                xlabel('angle(U_m,W,U_S_d) (^o)','FontSize',fsa)
                ylabel('z (m)','FontSize',fsa)
                % set(gca,'YDir','reverse')
                pltnam=strcat('U_1_0_N =',num2str(U10N(jU10p)), ...
                    ' m s^-^1',': p',num2str(j_mono_spec));
                title(pltnam,'FontSize',fsa)
                set(gca,'FontSize',fsa)

                % sfc_angles=(180/pi)*angle(Wtot(1:25))

            end
                     % end loop U10N
        end

        % end loop phi_lat
    end
           
    % end loop mono vs. spec. parameterizations
end

figure(10)

ustar_o=ustar_o_U10N(1)
tau0_o=rho_o*ustar_o^2
% ustar_o=ustar_o_U10N(jU10)
%  remove pressure stress from wave-drift model
% b0tau0=0.01
% b0_ustar_o=sqrt(b0tau0/rho_o)
% ustar_o=ustar_o-b0_ustar_o

max(lambda_k1)
min(lambda_k1)

b0_tauw=tauw_U10N/tau0_o
b0_tau0=b0*tau0_U10N/tau0_o
b0
% b0p=tau0_U10N/tau0_o

%%

%  nested function - COARE 3.5 iteration
function dustar = get_ustar(ustara)

    z0x=gamma*nu/ustara+alphax*ustara^2/g;
    CDN10sqrt=kappa/log(10/z0x);
    ustarx=CDN10sqrt*U10Nx;

    dustar=ustarx-ustara;

end




end
