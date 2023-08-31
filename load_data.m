function [Q0, Rate, theta, H, C0, Pd_theta, a] = load_data(Nt,M,K,SNR_dB)

SNR_dB_vec = -2:2:12;

%% load comms
system = strcat(num2str(Nt),'x',num2str(M),'x',num2str(K));
data = load(strcat('./sim_data/data_comm_',system,'.mat'));
H = data.H;
SNR_indx = find(SNR_dB_vec == SNR_dB);
Q0 = data.Q0_all(:,:,:,:,SNR_indx);
Rate = data.R_all(:,:,SNR_indx);

%% load sensing
data_Cbar = load(strcat('./sim_data/data_sens_',system,'.mat'));
C0 = data_Cbar.Cbar(:,:,:,SNR_indx);
Pd_theta = data_Cbar.Pd_theta;
theta = data_Cbar.theta;
a = data_Cbar.a;

end % EOF