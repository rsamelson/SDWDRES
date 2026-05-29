
clear;close all

fname='SDWDRES_z_U10_lat_3.mat';
load(fname);
% 
% fname=strcat('SDWDRES_z_U10_lat_',num2str(j_mono_spec),'.mat');
% save(fname,'zz_out','U10N_WD','phi_lat_array' ...
%     ,'ustar_o_U10N','L_p','Psi_int' ...
%     ,'Um_SZ','Wcap_SZ','USd_SZ','Uwd_SZ','US_SZ')

jphi_n=[5 13 21];
njphi=length(jphi_n);
phi_lat_array(jphi_n)

jU10N_n=[16 28 48 60];
nU10N=length(jU10N_n)
U10N_WD(jU10N_n)

fsa=14;

fig_last=figure;

for jjphi=1:njphi

    jphi=jphi_n(jjphi);
    phi_str=num2str(phi_lat_array(jphi));

    subplot(3,njphi,jjphi)
    %
    plot(U10N_WD,abs(Um_SZ(1,:,jphi)),'-b','LineWidth',2)
    hold on
    plot(U10N_WD,abs(Wcap_SZ(1,:,jphi)),'-r','LineWidth',2)
    plot(U10N_WD,abs(USd_SZ(1,:,jphi)),'--r','LineWidth',2)
    plot(U10N_WD,abs(Uwd_SZ(1,:,jphi)),'-g','LineWidth',1)
    plot([0 20],0.03*[0 20],':k')
    xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
    ylabel('Surface drift speed (m s^-^1)','FontSize',fsa)
    pltname=strcat(phi_str,'^o N latitude')
    title(pltname,'FontSize',fsa)
    set(gca,'FontSize',fsa)

    subplot(3,njphi,njphi+jjphi)
    %
    plot(U10N_WD,(180/pi)*angle(Um_SZ(1,:,jphi)),'-b','LineWidth',2)
    hold on
    plot(U10N_WD,(180/pi)*angle(Wcap_SZ(1,:,jphi)),'-r','LineWidth',2)
    plot(U10N_WD,(180/pi)*angle(USd_SZ(1,:,jphi)),'--r','LineWidth',2)
    plot(U10N_WD,(180/pi)*angle(Uwd_SZ(1,:,jphi)),'-g','LineWidth',1)
    plot([0 20],-15*[1 1],':k')
    axis([0 20 -45 0])
    xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
    ylabel('Surface drift angle (^o)','FontSize',fsa)
    pltname=strcat(phi_str,'^o N latitude')
    title(pltname,'FontSize',fsa)
    set(gca,'FontSize',fsa)

end
%
subplot(3,njphi,2*njphi+1)
%
% Resio et al 1999
g=9.81;
Ur=0.516*U10N_WD.^(1.244);
Hs_fd=0.21*Ur.^2/g;
%
plot(U10N_WD,sqrt(Psi_int),'-b','LineWidth',2)
hold on
plot(U10N_WD,2*sqrt(2*Psi_int),'-g','LineWidth',2)
plot(U10N_WD,Hs_fd,'-r','LineWidth',1)
axis([0 20 0 8])
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('\zeta_r_m_s, H_\zeta_r_m_s, H_s (m)','FontSize',fsa)
pltname=strcat('Phillips (1985) + Resio et al. (1999)')
title(pltname,'FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(3,njphi,2*njphi+2)
%
plot(U10N_WD,L_p,'-b','LineWidth',2)
axis([0 20 0 300])
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('L_p (m)','FontSize',fsa)
pltname=strcat('Resio et al. (1999)')
title(pltname,'FontSize',fsa)
set(gca,'FontSize',fsa)
%
subplot(3,njphi,2*njphi+3)
%
plot(U10N_WD,sqrt(ustar_o_U10N./squeeze(US_SZ(1,:,jphi))),'-b','LineWidth',2)
axis([0 20 0 1])
xlabel('U_1_0_N (m s^-^1)','FontSize',fsa)
ylabel('La_t=(u_*/U_S_0)^1^/^2','FontSize',fsa)
pltname=strcat('COARE 3.5')
title(pltname,'FontSize',fsa)
set(gca,'FontSize',fsa)


for jjU10N=1:nU10N

    fig_last=figure;

    jU10N=jU10N_n(jjU10N);
    U10N_str=num2str(U10N_WD(jU10N));

    for jjphi=1:njphi

        jphi=jphi_n(jjphi);
        phi_str=num2str(phi_lat_array(jphi));

    
        subplot(2,njphi,jjphi)
        %
        plot(abs(Um_SZ(:,jU10N,jphi)),zz_out,'-b','LineWidth',2)
        hold on
        plot(abs(Wcap_SZ(:,jU10N,jphi)),zz_out,'-r','LineWidth',2)
        plot(abs(USd_SZ(:,jU10N,jphi)),zz_out,'--r','LineWidth',2)
        plot(abs(Uwd_SZ(:,jU10N,jphi)),zz_out,'-g','LineWidth',1)
        plot(abs(Um_SZ(1,jU10N,jphi)),zz_out(1),'.b','MarkerSize',16)
        plot(abs(Wcap_SZ(1,jU10N,jphi)),zz_out(1),'.r','MarkerSize',16)
        plot(abs(USd_SZ(1,jU10N,jphi)),zz_out(1),'.r','MarkerSize',16)
        plot(abs(Uwd_SZ(1,jU10N,jphi)),zz_out(1),'.g','MarkerSize',16)
        % axis([0 0.035*U10N_WD(jU10N) -15 0])
        xlabel('Surface drift speed (m s^-^1)','FontSize',fsa)
        ylabel('Depth (m)','FontSize',fsa)
        pltname=strcat(U10N_str,' m s^-^1 (',phi_str,'^o N)')
        title(pltname,'FontSize',fsa)
        set(gca,'FontSize',fsa)
    
        subplot(2,njphi,njphi+jjphi)
        %
        plot((180/pi)*angle(Um_SZ(:,jU10N,jphi)),zz_out,'-b','LineWidth',2)
        hold on
        plot((180/pi)*angle(Wcap_SZ(:,jU10N,jphi)),zz_out,'-r','LineWidth',2)
        plot((180/pi)*angle(USd_SZ(:,jU10N,jphi)),zz_out,'--r','LineWidth',2)
        plot((180/pi)*angle(Uwd_SZ(:,jU10N,jphi)),zz_out,'-g','LineWidth',1)
        plot((180/pi)*angle(Um_SZ(1,jU10N,jphi)),zz_out(1),'.b','MarkerSize',16)
        plot((180/pi)*angle(Wcap_SZ(1,jU10N,jphi)),zz_out(1),'.r','MarkerSize',16)
        plot((180/pi)*angle(USd_SZ(1,jU10N,jphi)),zz_out(1),'.r','MarkerSize',16)
        plot((180/pi)*angle(Uwd_SZ(1,jU10N,jphi)),zz_out(1),'.g','MarkerSize',16)
        % axis([-120 0 -15 0])
        xlabel('Surface drift angle (^o)','FontSize',fsa)
        ylabel('Depth (m)','FontSize',fsa)
        pltname=strcat(U10N_str,' m s^-^1 (',phi_str,'^o N)')
        title(pltname,'FontSize',fsa)
        set(gca,'FontSize',fsa)

    end

    fig_last=figure;

    jU10N=jU10N_n(jjU10N);
    U10N_str=num2str(U10N_WD(jU10N));

    for jjphi=1:njphi

        jphi=jphi_n(jjphi);
        phi_str=num2str(phi_lat_array(jphi));

        subplot(2,njphi,jjphi)
        %
        semilogy(abs(Um_SZ(2:end,jU10N,jphi)),-zz_out(2:end),'-b','LineWidth',2)
        hold on
        semilogy(abs(Wcap_SZ(2:end,jU10N,jphi)),-zz_out(2:end),'-r','LineWidth',2)
        semilogy(abs(USd_SZ(2:end,jU10N,jphi)),-zz_out(2:end),'--r','LineWidth',2)
        semilogy(abs(Uwd_SZ(2:end,jU10N,jphi)),-zz_out(2:end),'-g','LineWidth',1)
        %
        semilogy(abs(Um_SZ(1,jU10N,jphi)),-zz_out(2),'.b','MarkerSize',16)
        semilogy(abs(Wcap_SZ(1,jU10N,jphi)),-zz_out(2),'.r','MarkerSize',16)
        semilogy(abs(USd_SZ(1,jU10N,jphi)),-zz_out(2),'.r','MarkerSize',16)
        semilogy(abs(Uwd_SZ(1,jU10N,jphi)),-zz_out(2),'.g','MarkerSize',16)
        % axis([0 0.035*U10N_WD(jU10N) -15 0])
        set(gca,'YDir','reverse')
        xlabel('Surface drift speed (m s^-^1)','FontSize',fsa)
        ylabel('Depth (m)','FontSize',fsa)
        pltname=strcat(U10N_str,' m s^-^1 (',phi_str,'^o N)')
        title(pltname,'FontSize',fsa)
        set(gca,'FontSize',fsa)

        subplot(2,njphi,njphi+jjphi)
        %
        semilogy((180/pi)*angle(Um_SZ(2:end,jU10N,jphi)),-zz_out(2:end),'-b','LineWidth',2)
        hold on
        semilogy((180/pi)*angle(Wcap_SZ(2:end,jU10N,jphi)),-zz_out(2:end),'-r','LineWidth',2)
        semilogy((180/pi)*angle(USd_SZ(2:end,jU10N,jphi)),-zz_out(2:end),'--r','LineWidth',2)
        semilogy((180/pi)*angle(Uwd_SZ(2:end,jU10N,jphi)),-zz_out(2:end),'-g','LineWidth',1)
        %
        semilogy((180/pi)*angle(Um_SZ(1,jU10N,jphi)),-zz_out(2),'.b','MarkerSize',16)
        semilogy((180/pi)*angle(Wcap_SZ(1,jU10N,jphi)),-zz_out(2),'.r','MarkerSize',16)
        semilogy((180/pi)*angle(USd_SZ(1,jU10N,jphi)),-zz_out(2),'.r','MarkerSize',16)
        semilogy((180/pi)*angle(Uwd_SZ(1,jU10N,jphi)),-zz_out(2),'.g','MarkerSize',16)
        % axis([-120 0 -15 0])
        set(gca,'YDir','reverse')
        xlabel('Surface drift angle (^o)','FontSize',fsa)
        ylabel('Depth (m)','FontSize',fsa)
        pltname=strcat(U10N_str,' m s^-^1 (',phi_str,'^o N)')
        title(pltname,'FontSize',fsa)
        set(gca,'FontSize',fsa)

    end

end

for jf=1:fig_last.Number
    figure(jf);
    pause;
end
