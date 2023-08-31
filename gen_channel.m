clear all;
clc;

gen_comm = 1;
gen_sensing = 1;

%% System parameters
Nt = 8; % number of transmit antennas
M = 2; % number of users
Nrf = M; % number of RF chains
K = 2; % number of subcarriers
system_config = strcat(num2str(Nt),'x',num2str(M),'x',num2str(K));
SNR_vec = -2:2:12;
init_scheme = 1; % 0: random, 1: ZF BF

Nr = 1; % number of receive antennas
UPA = 0;
Nray = 4; % number of rays in each cluster
gamma = sqrt((Nt*Nr)/(Nray)); %normalization factor
fc = 28e9; B = 1e9; delay_max = 20e-9;

if gen_comm == 1
    %% Generate channels
    n_trial = 5;
    H = zeros(Nt,M,K,n_trial);
    for ii = 1:n_trial
        ii
        for m = 1:M
            % randomly generate azimuth and elevation angles
            phi_AoD = unifrnd(0,2*pi,1,Nray);
            theta_AoD = unifrnd(-pi/2,pi/2,1,Nray);
            Delay = rand(Nray,1)*delay_max;
            alpha = normrnd(0,sqrt(1/2),[Nray,1]) + 1i*normrnd(0,sqrt(1/2),[Nray,1]);
            for k = 1:K
                fk = fc + B/K * (k-1 - (K-1)/2);
                hk = zeros(Nr,Nt);
                for j = 1:Nray
                    if UPA == 1
                        %Ntv = sqrt(Nt);
                        Ntv = 4;
                        Nth = Nt/Ntv;
                        a = array_response_UPA(Ntv,Nth,Nt,fk,fc,phi_AoD(j),theta_AoD(j));
                    else
                        a = array_response_ULA(Nt,fk,fc,phi_AoD(j));
                    end
                    hk = hk + alpha(j) * exp(-1i*2*pi*fk*Delay(j)) * a';
                end
                H(:,m,k,ii) = hk * gamma;
                %H(:,m,k,ii) = 1/sqrt(2) * (randn(Nt,1) + 1i*randn(Nt,1));
            end
        end
    end
    
    %% Solve effective precoders and rates
    [Q0_all, R_all] = gen_q(Nt,M,K,H,n_trial,SNR_vec,init_scheme);
    save(strcat('./sim_data/data_comm_',system_config,'.mat'));
end

if gen_sensing == 1
    %% Radar detection angles
    delta = pi/180;
    theta = -pi/2:delta:pi/2;
    target_DoA = [-pi/3,-pi/6,pi/6,pi/3];
    beam_width = 5;

    %% Ideal beampattern design
    l = ceil((target_DoA+pi/2*ones(1,length(target_DoA)))/(delta)+ones(1,length(target_DoA)));
    Pd_theta = zeros(length(theta),1);
    for ii = 1:length(target_DoA)
        Pd_theta(l(ii)-(beam_width-1)/2:l(ii)+(beam_width-1)/2,1) = ones(beam_width,1);
    end
    
    %% Solve Cbar
    n_SNR = length(SNR_vec);
    Cbar = zeros(Nt,Nt,K,n_SNR);
    a = zeros(Nt,length(theta),K);
    for k = 1:K
        fk = fc + B/K * (k-1 - (K-1)/2);
        ak = exp(1i* pi* fk/fc * ( [1:Nt].'.*sin(theta)));
        a(:,:,k) = ak;
        for ss = 1:n_SNR
            Cbar(:,:,k,ss) = solve_Cbar(Pd_theta,Nt,ak,theta,db2pow(SNR_vec(ss)));
        end
    end
    save(strcat('./sim_data/data_sens_',system_config,'.mat'));
end
