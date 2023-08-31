% need: Q0, Rate, theta, H, C0, Pd_theta, a
clear all
SNR_dB_vec = -2:2:12;
%% system parameters
Nt = 16; % number of transmit antennas
M = 4; % number of users
K = 64; % number of subcarriers

if Nt == 16
    %% load comms 1
    system = strcat(num2str(Nt),'x',num2str(M),'x',num2str(K));
    data1 = load(strcat('./sim_data/data_comm_',system,'.mat'));
    H = data1.H;
    Q0_all = data1.Q0_all;
    R_all = data1.R_all;
    
    %% load sensing
    data_Cbar = load(strcat('./sim_data/data_sens_',system,'.mat'));
    Cbar = data_Cbar.Cbar;
    Pd_theta = data_Cbar.Pd_theta;
    theta = data_Cbar.theta;
    a = data_Cbar.a;
else
    %% load comms 1
    system = strcat(num2str(Nt),'x',num2str(M),'x',num2str(K));
    data1 = load(strcat('./sim_data/data_comm_',system,'.mat'));
    H1 = data1.H;
    Q01 = data1.Q0_all;
    Rate1 = data1.R_all;
    
    %% load comms 2
    system = strcat(num2str(Nt),'x',num2str(M),'x',num2str(K));
    data2 = load(strcat('./sim_data/data_comm_more_',system,'.mat'));
    H2 = data2.H;
    Q02 = data2.Q0_all;
    Rate2 = data2.R_all;
    
    %% merge data
    H = cat(4,H1,H2);
    Q0_all = cat(4,Q01,Q02);
    R_all = cat(2,Rate1,Rate2);
    
    %% load sensing
    data_Cbar = load(strcat('./sim_data/data_sens_more_',system,'.mat'));
    Cbar = data_Cbar.Cbar;
    Pd_theta = data_Cbar.Pd_theta;
    theta = data_Cbar.theta;
    a = data_Cbar.a;
end
data_name = strcat('./sim_data/data_',system,'.mat');
save(data_name, 'H','Q0_all', 'R_all', 'theta', 'Cbar', 'Pd_theta', 'a');