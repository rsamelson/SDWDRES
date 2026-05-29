function SDWDRES_fn_z_U10_lat

% Compute solutions to Samelson and Zippel (2026)
%  wave-forced wind + wind-wave drift model
%  and save output in Matlab .mat file

clear; close all

%%
%  specify output depths (m)
zz_out=[0:-0.1:-1 -1.25:-0.25:-5 -5.5:-0.5:-10 -11:-1:-25];
nzz_out=length(zz_out);

%  10-m neutral wind array (m/s)
U10N=[0.25:0.25:20];
nU10=length(U10N);

%  Set latitudes for Coriolis parameter
%   (for re-scaling dimensionless solution)
%   NH latitudes (degrees) only; for SH latitudes,
%   use complex conjugates of NH velocities.
%   Note:  phi_w is used for wave-correction factor in wind-drift model
%           so use phi_lat for latitude
phi_lat_array=[5:2.5:90];
nphi_lat=length(phi_lat_array);

%  Set up output arrays
%   Arrays are indexed by (depth,10-m wind,north latitude)
%   Velocities W (m/s) are in complex form, W = U + iV,
%      with x-axis aligned positive downwind:
%       U = Re{W} = downwind (x) velocity
%       V = Im{W} = cross-wind (y) velocity
%   For southern hemisphere (S) latitudes, use complex conjugates
%    of computed velocities:
%       W_{Southern hemisphere} = conj(W) = U - iV
%
%  Um_SZ:   Total wave-forced wind drift + dynamic wave drift
%             from eq. (70) of Samelson and Zippel (2026).
Um_SZ=zeros(nzz_out,nU10,nphi_lat);
%
%  Wcap_SZ: Wave-forced wind drift
%             from eq. (63) of Samelson and Zippel (2026).
Wcap_SZ=zeros(nzz_out,nU10,nphi_lat);
%
%  USd_SZ:   Dynamic wave drift
%             from eq. (58) of Samelson and Zippel (2026).
USd_SZ=zeros(nzz_out,nU10,nphi_lat);
%
%  Uwd_SZ:   Wind-drift-only wind drift
%             (b0=0; no wave-correlated pressure forcing)
%             obtained by setting b_0=W_F=0 in eqs. (A22)-(A26)
%             of Samelson and Zippel (2026).
Uwd_SZ=zeros(nzz_out,nU10,nphi_lat);
%
%  US_SZ:   Kinematic wave-drift
%             from the spectral integral of eq. (56)
%             of Samelson and Zippel (2026).
US_SZ=zeros(nzz_out,nU10,nphi_lat);
%
%  Also saved on output:
%    ustar_o_U10N (m/s): Ocean (water) friction velocity
%    L_p (m):      Resio et al. (1999) fully-developed sea peak wavelength
%    Psi_int (m):  Integrated Phillips (1985) mean-square displacement
%                    spectrum with Resio et al. (1999) peak wavelength L_p
%                    as lower integration limit and derived upper
%                    integration limit k_1 as in Samelson and Zippel (2026)

%%
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

%  Phillips parameters
%  Phillips constant (1958)
alpha_P=0.0083;
%  equilibrium (1985)
gamma_beta3_I3pp1=10^-3
%
% betaIp=1.2*10^-2
%  Phillips 1985 p. 526; p=1/2
betaIp=3*10^-2
% %
%  factor to reduce b0 by setting integration limit:
%   smallest wavelength (upper limit k1 of wavenumber integral) controls b0
lambda_k1=0.15*U10N;

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

%%  SDRO - Solve Samelson and Zippel (2026) model equations

ustar=ustar_a;

% Resio et al 1999
Ur=0.516*U10N.^(1.244);
Hs_fd=0.21*Ur.^2/g;
% Hs_fd=0.056*U10N.^(2.488)/g;
%  c_P=u_r
L_p=2*pi.*Ur.^2/g;
ustar_r=Ur/24.18;
z0_r=0.015*ustar_r.^2/g;
CDN_r=(kappa./log(10./z0_r)).^2;
%  Stokes drift mean
k_r=2*pi./L_p;
sigma_r=sqrt(g*k_r);
cp_r=sigma_r./k_r;
Us_r=0.5*cp_r.*(Hs_fd*pi./L_p).^2;

% %  "bulk-flux" wave state (Edson et al. 2013 wave-age relation)
% ustar_cp=max(10^-9,0.03*(U10N-2)/8);
% %   inferred Hs
% Dcap=0.09;
% Hs=z0_rough./(Dcap*ustar_cp.^2);
% %
% cp_s0=ustar./ustar_cp;
% k_s0=g./cp_s0.^2;
% sigma_s0=sqrt(g*k_s0);
% lambda_s0=2*pi./k_s0;
% %  Stokes wave mean
% Us0=0.5*cp_s0.*(0.5*alpha_c/Dcap).^2;
% %
% %   inferred wavelength
% %    arbitrary factor of 1/2
% lambda_s=0.5*lambda_s0;
% %  wavenumber and frequency
% k_s=2*pi./lambda_s;
% sigma_s=sqrt(g*k_s);
% cp_s=sigma_s/k_s;
% %  Stokes wave mean
% Us=0.5*cp_s.*(alpha_c/Dcap).^2;
% 


%  T_eq (T_S) estimates
% 
% fully-developed sea
b2=5
b0_fd=min(1,b2.*0.5*k_r.*Hs_fd);
TS_fd=0.125*rho0*g*Hs_fd.^2./(b0_fd.*sqrt(g./k_r).*rho_a.*ustar.^2);
P_a_fd=b0_fd.*rho_a.*ustar.^2./(0.5*k_r.*Hs_fd);
% P_a_fd=b0*rho_a*ustar.^2./(0.5*k_r.*Hs_fd);
P_a_fd(Hs_fd<=0)=NaN;
T_Sfd=rho0*sqrt(g)*0.5*Hs_fd./(sqrt(k_r).*P_a_fd);
T_Sfd(Hs_fd<=0)=NaN;

% % bulk-flux
% b0_bf=min(1,b2.*0.5*k_s0.*Hs);
% TS_bf=0.125*rho0*g*Hs.^2./(b0_bf.*sqrt(g./k_s0).*rho_a.*ustar.^2);
% P_a=b0_bf.*rho_a.*ustar.^2./(0.5*k_s0.*Hs);
% % P_a=b0*rho_a*ustar.^2./(0.5*k_s0.*Hs);
% P_a(Hs<=0)=NaN;
% T_S0=rho0*sqrt(g)*0.5*Hs./(sqrt(k_s0).*P_a);
% T_S0(Hs<=0)=NaN;


%  spectral wavenumber kk
%   lower integration limit k0 set equal to Resio et al. 1999
%   upper integration limit k1 (corresponding to 
%   smallest wavelength lambda_k1) set to control b0
%
lambda_k=[0.05:0.01:4.99 5:0.05:50 51:0.25:100 101:2000];
% lambda_k=[0.1:0.05:50 51:0.25:100 101:2000];
kk=2*pi./lambda_k;

Tk_j=nan(length(kk),nU10);
Tk_j_plt=nan(length(kk),nU10);

for jU10=1:nU10
    for jkk=1:length(kk)
        if(lambda_k(jkk) < L_p(jU10))
            Tk_j_plt(jkk,jU10)= 25*sqrt(g)*ustar(jU10)^-2*kk(jkk)^(-3/2);
            if(lambda_k(jkk) > lambda_k1(jU10))
              Tk_j(jkk,jU10)= 25*sqrt(g)*ustar(jU10)^-2*kk(jkk)^(-3/2);
            end
        end
    end    
end

%  get profile of spectral terms for wave-drift solution
%  and for particular solution of wave-forced wind-drift model
zzk=[0:-0.25:-25];
zzk(1)=-0.01;
nzzk=length(zzk);
% 

%  use P3 in wave-forced wind-drift model
for j_mono_spec=3:3
% for j_mono_spec=1:4

    j_mono_spec

    %  log-Ekman layer parameters
    if(j_mono_spec==1)
        % param 1
        z0o_U10N=10.^(-5+7*(1-exp(-U10N/10)));
        phi_w_U10N=exp(-U10N/10);
    elseif(j_mono_spec==2)
        % param 2
        z0o_U10N=10.^(-6+7.5*(1-exp(-U10N/10)));
        phi_w_U10N=exp(-U10N/15);
   elseif(j_mono_spec==3)
        % param 3
        z0o_U10N=10.^(-5+7*(1-exp(-U10N/5)));
        phi_w_U10N=0.05+0.95*exp(-U10N/5);
    elseif(j_mono_spec==4)
        % param 4
        z0o_U10N=10.^(-6+7.5*(1-exp(-U10N/5)));
        phi_w_U10N=0.05+0.95*exp(-U10N/7.5);
    end

    %  latitude loop
    for jphi_lat=1:nphi_lat

        phi_lat=phi_lat_array(jphi_lat)
        % %
        fcor=2*Omega*sind(phi_lat)
        
        %  SDRO
        
        Us_j=complex(zeros(nzzk,nU10),zeros(nzzk,nU10));
        GfUs_j=complex(zeros(nzzk,nU10),zeros(nzzk,nU10));
        Up_j=zeros(nzzk,nU10);
        Px_zint=zeros(1,nU10);
        b0_U10N=zeros(1,nU10);
        Psi_int=zeros(1,nU10);

        % whos Us_j
        
        for jU10=1:nU10
            for jkk=2:length(kk)  
                if(  lambda_k(jkk) > lambda_k1(jU10) ...
                  && lambda_k(jkk) < L_p(jU10))
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
                    Psi_int(jU10)=Psi_int(jU10) ...
                        + (kk_j^(-3.5))*(kk(jkk-1)-kk(jkk));
                end
            end
        
            Us_j(:,jU10)=2*betaIp*ustar(jU10)*Us_j(:,jU10);
            GfUs_j(:,jU10)=2*betaIp*ustar(jU10)*GfUs_j(:,jU10);
            Up_j(:,jU10)=2*betaIp*ustar(jU10)*Up_j(:,jU10);
            Px_zint(jU10)=betaIp*ustar(jU10)*Px_zint(jU10);
            Psi_int(jU10)=(betaIp*ustar(jU10)/sqrt(g))*Psi_int(jU10);

        end
        
        % Us_j(Us_j==0)=NaN;
        % Up_j(Up_j==0)=NaN;
        
        %  Compute total wave-correlated pressure fraction of stress as
        %  check
        b0_U10N=(rho_o/rho_a)*Px_zint./ustar.^2;
        
        
        
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
        % % UprofU10N=complex(nan(nzprof,nU10),nan(nzprof,nU10));
        % U10profU10N=nan(nzprof,nU10);
        zprofU10N=nan(nzprof,nU10);
        % epsprofU10N=nan(nzprof,nU10);
        Wtot_profU10N=complex(nan(nzprof,nU10),nan(nzprof,nU10));
        Wcap_profU10N=complex(nan(nzprof,nU10),nan(nzprof,nU10));
        WSd_profU10N=complex(nan(nzprof,nU10),nan(nzprof,nU10));
        W0_profU10N=complex(nan(nzprof,nU10),nan(nzprof,nU10));
        WS_profU10N=complex(nan(nzprof,nU10),nan(nzprof,nU10));

        %%
        
        
        %%  Solve wind-drift model
        %    Eqs.  (12), (31), (32), (36)-(38) with
        %     j_mono_spec=1 => monochromatic calibration parameters:  (50)-(51)
        %     j_mono_spec=2 => spectral calibration parameters:  (52)-(53)
       
        %  U10N loop
        for jU10=1:nU10

            %
            phi_w=phi_w_U10N(jU10);
            z0_o_dim=z0o_U10N(jU10);
            ustar_o=ustar_o_U10N(jU10);
            %     ustar_o=alpha_ao*ustar_a(jU10);

            D_o=ustar_o/fcor;

            %%  log-Ekman layer
            % %  dimensionless depth z/D, D=ustar_o/f
            pkappa=phi_w*kappa;
            %  check for viscous sublayer
            z0_o_dim=max(z0_o_dim,nu_MKS/(pkappa*ustar_o));
            % ocean - dimensionless
            z0_o=z0_o_dim*fcor/ustar_o;
            %
            %
            %  Set depth-coordinate array based on WDHES
            %   profile of dimensionless eddy viscosity Av/Av0,
            %      where Av0=ustar_o*D=ustar_o^2/f
            %   but use constant Av to get solution
            %   in terms of single modified Bessel function K_0,
            %   for which particular solution is simple to compute.
            % 
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
            dz3=0.002*(z3+25);
            zz3=[z3-dz3:-dz3:-25];
            zz=[zza zz1 zz2 zz3];
    
%  wave-breaking forcing
            % particular solutions Ap*K0+Bp*I0
            %
            xi_K0=(1+i).*sqrt(2.*(z0_o-zz)./pkappa);
            %  integrate to get Ap, Bp
            nzprof1=min(length(zz),nzprof);
            % 
            Ap=complex(zeros(1,nzprof1),zeros(1,nzprof1));
            Bp=complex(zeros(1,nzprof1),zeros(1,nzprof1));
            GfWSd=complex(zeros(1,nzprof1),zeros(1,nzprof1));
            % 
            zprof_K0=D_o*zz(1:nzprof1);
            zzk_int=[zprof_K0(1) zzk zzk(end)-1 zprof_K0(end)-1];
            GfUs_int=[GfUs_j(1,jU10) ...
                      GfUs_j(:,jU10)' 0+0*i 0+0*i];
            GfWSd=interp1(zzk_int,GfUs_int,zprof_K0);
            for jz=nzprof1-1:-1:1
                Ap(jz)=Ap(jz+1) ...
                  + 0.5*(xi_K0(jz)*besseli(0,xi_K0(jz))*GfWSd(jz) ...
                        +xi_K0(jz+1)*besseli(0,xi_K0(jz+1))*GfWSd(jz+1))...
                      *(xi_K0(jz)-xi_K0(jz+1));
                Bp(jz)=Bp(jz+1) ...
                  + 0.5*(xi_K0(jz)*besselk(0,xi_K0(jz))*GfWSd(jz) ...
                        +xi_K0(jz+1)*besselk(0,xi_K0(jz+1))*GfWSd(jz+1))...
                      *(xi_K0(jz)-xi_K0(jz+1));
            end
            % whos Ap
            % max(abs(Ap))
            Ap=-i*Ap;
            Bp=i*Bp;
            Wp=Ap.*besselk(0,xi_K0(1:nzprof1)) ...
              +Bp.*besseli(0,xi_K0(1:nzprof1));

            % whos Wp
            % max(abs(Wp))
            % min(abs(Wp))

            %  homogeneous solution
            xi0=(1+i)*sqrt(2*z0_o/pkappa);
            Kcap0z0=besselk(0,xi0);
            Kcap1z0=besselk(1,xi0);
            Icap1z0=besseli(1,xi0);
            %
            %  Remove pressure-stress fraction here.
            % b0_U10N=(rho_o/rho_a)*Px_zint./ustar.^2;
            b0=b0_U10N(jU10);
            %
            Ah=(exp(-i*0.25*pi) ...
                *sqrt(ustar_o^2*D_o/(pkappa*z0_o_dim))*(1-b0) ...
                -Ap(1)*Kcap1z0+Bp(1)*Icap1z0)/Kcap1z0;
            %
            %  wave-forced wind-drift profile
            Wcap=Ah*besselk(0,xi_K0(1:nzprof1))+Wp;

            %  dynamic wave-drift profile
            WSd=complex(zeros(1,nzprof1),zeros(1,nzprof1));
            Us_int=[Us_j(1,jU10) ...
                      transpose(Us_j(:,jU10)) 0+0*i 0+0*i];
            WSd=interp1(zzk_int,Us_int,zprof_K0);

            %  total wave-forced wind-drift + wave-drift
            Wtot=Wcap+WSd;

            %  wind-drift model with b0=0 (no pressure forcing)
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
            W0=ustar_o*Bcap_K0(1:nzprof1);
            % Uprof_K0=ustar_o*Bcap_K0(1:nzprof1);

            %  kinematic wave-drift profile
            WS=complex(zeros(1,nzprof1),zeros(1,nzprof1));
            Up_int=[Up_j(1,jU10) ...
                transpose(Up_j(:,jU10)) 0+0*i 0+0*i];
            WS=interp1(zzk_int,Up_int,zprof_K0);

            %  fill arrays for output grid interpolation
            %  Wtot, Wcap, WSd, W0 are dimensional
            Wtot_profU10N(1:nzprof1,jU10)=Wtot(1:nzprof1);
            Wcap_profU10N(1:nzprof1,jU10)=Wcap(1:nzprof1);
            WSd_profU10N(1:nzprof1,jU10)=WSd(1:nzprof1);
            W0_profU10N(1:nzprof1,jU10)=W0(1:nzprof1);
            WS_profU10N(1:nzprof1,jU10)=WS(1:nzprof1);
            %  dimensional vertical grid (wind-speed dependent)
            zprofU10N(1:nzprof1,jU10)=D_o*zz(1:nzprof1);
                
            % end U10N loop
        end

        %%
        %  Interpolate to output vertical grid
        %  loop over U10N 
        for jU10=1:nU10
            Um_SZ(:,jU10,jphi_lat)= ...
                interp1(zprofU10N(1:nzprof1,jU10) ...
                       ,Wtot_profU10N(1:nzprof1,jU10),zz_out);
            Wcap_SZ(:,jU10,jphi_lat)= ...
                interp1(zprofU10N(1:nzprof1,jU10) ...
                       ,Wcap_profU10N(1:nzprof1,jU10),zz_out);
            USd_SZ(:,jU10,jphi_lat)= ...
                interp1(zprofU10N(1:nzprof1,jU10) ...
                       ,WSd_profU10N(1:nzprof1,jU10),zz_out);
            Uwd_SZ(:,jU10,jphi_lat)= ...
                interp1(zprofU10N(1:nzprof1,jU10) ...
                       ,W0_profU10N(1:nzprof1,jU10),zz_out);
            US_SZ(:,jU10,jphi_lat)= ...
                interp1(zprofU10N(1:nzprof1,jU10) ...
                ,WS_profU10N(1:nzprof1,jU10),zz_out);
            % 
            % %  reset surface velocity
            % Um_SZ(1,jU10,jphi_lat)=Wtot_profU10N(1,:);
            % Wcap_SZ(1,jU10,jphi_lat)i=Wcap_profU10N(1,:);
            % USd_SZ(1,jU10,jphi_lat)=WSd_profU10N(1,:);
            % Uwd_SZ(1,jU10,jphi_lat)=W0_profU10N(1,:);
        end

        % end loop over phi_lat_array
        % pause
    end

    %%  save output file
    %
    U10N_WD=U10N';
    fname=strcat('SDWDRES_z_U10_lat_',num2str(j_mono_spec),'.mat');
    save(fname,'zz_out','U10N_WD','phi_lat_array' ...
              ,'ustar_o_U10N','L_p','Psi_int' ...
              ,'Um_SZ','Wcap_SZ','USd_SZ','Uwd_SZ','US_SZ')
           
    % end loop over mono vs. spec. parameterizations
end



%%
%  nested function - COARE 3.5 neutral drag coefficient
    function dustar = get_ustar(ustara)
    %     gamma
    %     nu
    %     alphax
    %     g
    %     ustara
        z0x=gamma*nu/ustara+alphax*ustara^2/g;
        CDN10sqrt=kappa/log(10/z0x);
        ustarx=CDN10sqrt*U10Nx;
    
        dustar=ustarx-ustara;
    
    end

end
