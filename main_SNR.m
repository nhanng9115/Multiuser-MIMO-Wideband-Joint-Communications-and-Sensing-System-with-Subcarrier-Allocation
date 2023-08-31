clc;
clear all;
close all

%% system parameters
Nt = 8; % number of transmit antennas
M = 2; % number of users
K = 2; % number of subcarriers
plot_beampattern = 1; plot_MSE = 1; plot_rate = 1;

%% simulation parameters
n_chans = 1;
SNR_dB_vec = -2:2:12;
n_snr = length(SNR_dB_vec);
J = K/2;
CSI_error = 0;
HBF = 1;

T = 181;
Nrf = 6; % number of RF chains
F0 = exp(1i*randn(Nt,Nrf));
rho_vec = [0.2, 0.4];
for rr = 1:length(rho_vec)
    rho = rho_vec(rr);
    
    %% Initiate results holder for parallel running
    % rate
    rate_dig_propose = zeros(n_chans,n_snr);
    rate_hyb_propose = zeros(n_chans,n_snr);
    rate_dig_overlap = zeros(n_chans,n_snr);
    rate_hyb_overlap = zeros(n_chans,n_snr);
    rate_dig_non_overlap = zeros(n_chans,n_snr);
    rate_hyb_non_overlap = zeros(n_chans,n_snr);
    rate_dig_comm = zeros(n_chans,n_snr);
    rate_hyb_comm = zeros(n_chans,n_snr);
    rate_dig_propose_rand_mean = zeros(n_chans,n_snr);
    rate_hyb_propose_rand = zeros(n_chans,n_snr);
    % MSE
    MSE_dig_propose = zeros(n_chans,n_snr);
    MSE_hyb_propose = zeros(n_chans,n_snr);
    MSE_dig_overlap = zeros(n_chans,n_snr);
    MSE_hyb_overlap = zeros(n_chans,n_snr);
    MSE_dig_non_overlap = zeros(n_chans,n_snr);
    MSE_hyb_non_overlap = zeros(n_chans,n_snr);
    MSE_dig_comm = zeros(n_chans,n_snr);
    MSE_hyb_comm = zeros(n_chans,n_snr);
    % beampattern
    beam_benchmark = zeros(T,K,n_chans,n_snr);
    beam_dig_overlap = zeros(T,K,n_chans,n_snr);
    beam_hyb_overlap = zeros(T,K,n_chans,n_snr);
    beam_dig_non_overlap = zeros(T,J,n_chans,n_snr);
    beam_hyb_non_overlap = zeros(T,J,n_chans,n_snr);
    beam_dig_propose = zeros(T,J,n_chans,n_snr);
    beam_hyb_propose = zeros(T,J,n_chans,n_snr);
    
    %% start simulations for SNRs
    for ss = 1:n_snr
        ss
        SNR_dB = SNR_dB_vec(ss);
        Pt = db2pow(SNR_dB);
        % load data
        [Q0, Rate, theta, H, C0, Pd_theta, a] = load_data(Nt,M,K,SNR_dB);
        
        %% start simulations for channels
        for nn = 1:n_chans
            h = H(:,:,:,nn); % current channel
            
            % get Q and rate
            Q0_nn = Q0(:,:,:,nn);
            Rate_nn = Rate(:,nn);
            
            % choose best subcarriers
            [Omg, Omg_rand] = subcarrier_select(K,J,Q0_nn,C0);
            
            % obtain beampatterns, MSE, and rate
            [beam_dig_overlap(:,:,nn,ss), MSE_dig_overlap(nn,ss), rate_dig_overlap(nn,ss),...
                beam_dig_non_overlap(:,:,nn,ss), MSE_dig_non_overlap(nn,ss), rate_dig_non_overlap(nn,ss),...
                beam_dig_propose(:,:,nn,ss), MSE_dig_propose(nn,ss), rate_dig_propose(nn,ss), ...
                beam_hyb_overlap(:,:,nn,ss), MSE_hyb_overlap(nn,ss), rate_hyb_overlap(nn,ss),...
                beam_hyb_non_overlap(:,:,nn,ss), MSE_hyb_non_overlap(nn,ss), rate_hyb_non_overlap(nn,ss),...
                beam_hyb_propose(:,:,nn,ss), MSE_hyb_propose(nn,ss), rate_hyb_propose(nn,ss),...
                rate_dig_propose_rand_mean(nn,ss), rate_hyb_propose_rand(nn,ss),...
                rate_dig_comm(nn,ss),rate_hyb_comm(nn,ss), MSE_dig_comm(nn,ss), MSE_hyb_comm(nn,ss)]...
                = JCAS_design(Nt,M,K,C0,Q0_nn,Pt,Pd_theta,Omg,Omg_rand,rho,a,h,T,F0,HBF,CSI_error);
        end
    end
    
    %% setup for figures
    schemes = {'Communication only', 'Prop. JCAS', 'Prop. JCAS, random $\mathcal{J}$', 'Conv. JCAS, overlap', 'Conv. JCAS, nonoverlap'};
    
    %% plot beampattern ============================================================================================
    if plot_beampattern == 1
        
        % plot beampatterns at SNR = 12 dB
        SNR_to_show_beampatter = 12;
        SNR_indx = find(SNR_dB_vec == SNR_to_show_beampatter);
        
        
        %% average beampattern for all channels, digital -----------------------------------------------
        x_angles = theta*180/pi; % show on x-axis
        % compute beampattern gain
        %beam_benchmark_mean = mean(mean(beam_benchmark(:,:,:,SNR_indx),2),3);
        beam_dig_propose_mean = mean(mean(beam_dig_propose(:,:,:,SNR_indx),2),3);
        beam_dig_overlap_mean = mean(mean(beam_dig_overlap(:,:,:,SNR_indx),2),3);
        beam_dig_non_overlap_mean = mean(mean(beam_dig_non_overlap(:,:,:,SNR_indx),2),3);
        % plot figure
        figure(1)
        plot(x_angles,Pd_theta,'--b','LineWidth',0.5);hold on;
        plot(x_angles,beam_dig_propose_mean,'-','LineWidth',1.5);hold on;
        plot(x_angles,beam_dig_overlap_mean,':','LineWidth',1.5);hold on;
        plot(x_angles,beam_dig_non_overlap_mean,'-.','LineWidth',1.5);hold on;
        legend('Desired beampattern',schemes{2},schemes{4},schemes{5},'Location','Best','fontsize',11,'interpreter','latex')
        xlim([-90, 90])
        xticks([-90:30:90])
        xlabel('Angles $(\theta_t^{\circ})$','fontsize',12,'interpreter','latex');
        ylabel('Normalized beampatter','fontsize',12,'interpreter','latex');
        ylim([0 0.3])
        title('Average beampattern, digital');
        grid on
        
        %% average beampattern for all channels, HBF -----------------------------------------------
        % compute beampattern gain
        beam_hyb_overlap_mean = mean(mean(beam_hyb_overlap(:,:,:,SNR_indx),2),3);
        beam_hyb_non_overlap_mean = mean(mean(beam_hyb_non_overlap(:,:,:,SNR_indx),2),3);
        beam_hyb_propose_mean = mean(mean(beam_hyb_propose(:,:,:,SNR_indx),2),3);
        % plot figure
        figure(2)
        plot(x_angles,Pd_theta,'--b','LineWidth',0.5);hold on;
        plot(x_angles,beam_hyb_propose_mean,'-','LineWidth',1.5);hold on;
        plot(x_angles,beam_hyb_overlap_mean,':','LineWidth',1.5);hold on;
        plot(x_angles,beam_hyb_non_overlap_mean,'-.','LineWidth',1.5);hold on;
        legend('Desired beampattern',schemes{2},schemes{4},schemes{5},'Location','Best','fontsize',11,'interpreter','latex')
        xlim([-90, 90])
        xticks([-90:30:90])
        xlabel('Angles $(\theta_t^{\circ})$','fontsize',12,'interpreter','latex');
        ylabel('Normalized beampatter','fontsize',12,'interpreter','latex');
        ylim([0 0.3])
        grid on
    end
    
    
    
    %% plot rate ============================================================================================
    if plot_rate == 1
        %% fully digital --------------------------------------------------------------------------------------
        % compute rate
        rate_dig_comm_mean = mean(rate_dig_comm,1);
        rate_dig_propose_mean = mean(rate_dig_propose,1);
        rate_dig_propose_rand_mean = mean(rate_dig_propose_rand_mean,1);
        rate_dig_overlap_mean = mean(rate_dig_overlap,1);
        rate_dig_non_overlap_mean = mean(rate_dig_non_overlap,1);
        % plot figure
        figure(3)
        plot(SNR_dB_vec, rate_dig_comm_mean, '-k+','LineWidth',2,'MarkerSize',8); hold on;
        plot(SNR_dB_vec, rate_dig_propose_mean, '-r*','LineWidth',2,'MarkerSize',8); hold on;
        plot(SNR_dB_vec, rate_dig_propose_rand_mean, '--+','Color',[0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8); hold on;
        plot(SNR_dB_vec, rate_dig_overlap_mean, ':bo','LineWidth',2,'MarkerSize',8); hold on;
        plot(SNR_dB_vec, rate_dig_non_overlap_mean, '-.^','Color',[0.4660 0.6740 0.1880],'LineWidth',2,'MarkerSize',8); hold on;
        xlabel('SNR [dB]','fontsize',12,'interpreter','latex');
        ylabel('Achievable rate [bits/s/Hz]','fontsize',12,'interpreter','latex');
        xticks(SNR_dB_vec)
        legend(schemes{1},schemes{2},schemes{3},schemes{4},schemes{5},'Location','Best','fontsize',11,'interpreter','latex')
        grid on
        
        %% hybrid beamforming --------------------------------------------------------------------------------------
        % compute rate
        rate_hyb_comm_mean = mean(rate_hyb_comm,1);
        rate_hyb_propose_mean = mean(rate_hyb_propose,1);
        rate_hyb_propose_rand_mean = mean(rate_hyb_propose_rand,1);
        rate_hyb_overlap_mean = mean(rate_hyb_overlap,1);
        rate_hyb_non_overlap_mean = mean(rate_hyb_non_overlap,1);
        
        % plot figure
        figure(4)
        plot(SNR_dB_vec, rate_hyb_comm_mean, '-k+','LineWidth',2,'MarkerSize',8); hold on;
        plot(SNR_dB_vec, rate_hyb_propose_mean, '-r*','LineWidth',2,'MarkerSize',8); hold on;
        plot(SNR_dB_vec, rate_hyb_propose_rand_mean, '--+','Color',[0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8); hold on;
        plot(SNR_dB_vec, rate_hyb_overlap_mean, ':bo','LineWidth',2,'MarkerSize',8); hold on;
        plot(SNR_dB_vec, rate_hyb_non_overlap_mean,  '-.^','Color',[0.4660 0.6740 0.1880],'LineWidth',2,'MarkerSize',8); hold on;
        xlabel('SNR [dB]','fontsize',12,'interpreter','latex');
        ylabel('Achievable rate [bits/s/Hz]','fontsize',12,'interpreter','latex');
        xticks(SNR_dB_vec)
        legend(schemes{1},schemes{2},schemes{3},schemes{4},schemes{5},'Location','Best','fontsize',11,'interpreter','latex')
        grid on
    end
end