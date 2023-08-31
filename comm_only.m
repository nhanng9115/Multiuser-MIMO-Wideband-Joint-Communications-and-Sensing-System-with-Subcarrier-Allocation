function [rate_comm_dig_sub, rate_comm_dig_all, rate_comm_hyb_sub, rate_comm_hyb_all] = comm_only(Nt,M,K,Q0,Omg,h,HBF,F0,Pt)

%% digital
Q_dig_sub = zeros(Nt,M,K);
Q_dig_all = zeros(Nt,M,K);
for k = 1:K
    Q_dig_all(:,:,k) = Q0(:,:,k);
    % subcarrier selection: only use opt_radcom for k in Omg
    if ~ismember(k,Omg)
        Q_dig_sub(:,:,k) = Q0(:,:,k);
    end
end
% Compute rate
rate_comm_dig_sub = compute_rate_comm(h,Nt,M,K,Q_dig_sub);
rate_comm_dig_all = compute_rate_comm(h,Nt,M,K,Q_dig_all);

%% HBF
rate_comm_hyb_sub = 0;
rate_comm_hyb_all = 0;

if HBF == 1
    [F_sub, W_sub] = MO_AltMin(Q_dig_sub, F0, Pt);
    [F_all, W_all] = MO_AltMin(Q_dig_all, F0, Pt);
    Q_hyb_sub = zeros(Nt,M,K);
    Q_hyb_all = zeros(Nt,M,K);
    for k = 1:K
        Q_hyb_sub(:,:,k) = F_sub*W_sub(:,:,k);
        Q_hyb_all(:,:,k) = F_all*W_all(:,:,k);
    end
    
    % Compute rate
    rate_comm_hyb_sub = compute_rate_comm(h,Nt,M,K,Q_hyb_sub);
    rate_comm_hyb_all = compute_rate_comm(h,Nt,M,K,Q_hyb_all);
end