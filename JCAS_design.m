function [beam_dig_overlap, MSE_dig_overlap, rate_dig_overlap,...
    beam_dig_non_overlap, MSE_dig_non_overlap, rate_dig_non_overlap,...
    beam_dig_propose, MSE_dig_propose, rate_dig_propose, ...
    beam_hyb_overlap, MSE_hyb_overlap, rate_hyb_overlap,...
    beam_hyb_non_overlap, MSE_hyb_non_overlap, rate_hyb_non_overlap,...
    beam_hyb_propose, MSE_hyb_propose, rate_hyb_propose,...
    rate_dig_propose_rand, rate_hyb_propose_rand,...
    rate_dig_comm, rate_hyb_comm, MSE_dig_comm, MSE_hyb_comm, beam_benchmark] = ...
    JCAS_design(Nt,M,K,Cbar,Q0,Pt,Pd_theta,Omg,Omg_rand,rho,at,h,T,F0,HBF,error)

J = length(Omg);

%% digital -----------------------------------------------------------------------------------------------
beam_benchmark = zeros(T,K);
beam_dig_overlap = zeros(T,K);
beam_dig_non_overlap = zeros(T,K);
beam_dig_propose = zeros(T,K);
beam_dig_comm = zeros(T,K);

Q_dig_overlap = zeros(Nt,M,K);
Q_dig_non_overlap = zeros(Nt,M,K);
Q_dig_propose = Q0;
Q_dig_propose_rand = Q0;
Q_dig_comm = Q0;
% j = 1;
for k = 1:K
    beam_benchmark(:,k) = real(diag(at(:,:,k)'*Cbar(:,:,k)*at(:,:,k)))/(Nt*Pt);
    
    Qbar = Q0(:,:,k);
    
    %% overlap
    Q_dig_overlap(:,:,k) = opt_radcom(Nt,M,rho,Cbar(:,:,k),Qbar,Pt);
    beam_dig_overlap(:,k) = real(diag(at(:,:,k)'*Q_dig_overlap(:,:,k)*Q_dig_overlap(:,:,k)'*at(:,:,k)))/(Nt*Pt);
    
    %% proposed subcarrier selection: only use opt_radcom for k in Omg
    if ismember(k,Omg)
        Q_dig_propose(:,:,k) = Q_dig_overlap(:,:,k);
        beam_dig_propose(:,k) = real(diag(at(:,:,k)'*Q_dig_propose(:,:,k)*Q_dig_propose(:,:,k)'*at(:,:,k)))/(Nt*Pt);
    end
    
    %% random subcarrier: only use opt_radcom for k in Omg
    if ismember(k,Omg_rand)
        Q_dig_propose_rand(:,:,k) = Q_dig_overlap(:,:,k);
    end
    
    %% comm only
    beam_dig_comm(:,k) = real(diag(at(:,:,k)'*Q_dig_comm(:,:,k)*Q_dig_comm(:,:,k)'*at(:,:,k)))/(Nt*Pt);
    
    %% non-overlap: 1->J for sensing, J+1->K: for comm
    if k <= J
        %% sensing only: rho = 1; rho_bar = 0;
        Q_dig_non_overlap(:,:,k) = opt_radcom(Nt,M,1,Cbar(:,:,k),Qbar,Pt);
        beam_dig_non_overlap(:,k) = real(diag(at(:,:,k)'*Q_dig_non_overlap(:,:,k)*Q_dig_non_overlap(:,:,k)'*at(:,:,k)))/(Nt*Pt);
    else
        %% comm only: rho = 0; rho_bar = 1;
        Q_dig_non_overlap(:,:,k) = Qbar;
    end
end


% Compute rate
rate_dig_overlap = compute_rate(h,Nt,M,K,Q_dig_overlap,error);
rate_dig_non_overlap = compute_rate(h,Nt,M,K,Q_dig_non_overlap,error);
rate_dig_propose = compute_rate(h,Nt,M,K,Q_dig_propose,error);
rate_dig_propose_rand = compute_rate(h,Nt,M,K,Q_dig_propose_rand,error);
rate_dig_comm = compute_rate(h,Nt,M,K,Q_dig_comm,error);


% Compute MSE
MSE_dig_overlap = norm(Pd_theta - mean(beam_dig_overlap,2))^2/T;
beam_dig_propose( :, ~any(beam_dig_propose,1) ) = []; % remove non-used subcarriers
MSE_dig_propose = norm(Pd_theta - mean(beam_dig_propose,2))^2/T;
beam_dig_non_overlap( :, ~any(beam_dig_non_overlap,1) ) = []; % remove non-used subcarriers
MSE_dig_non_overlap = norm(Pd_theta - mean(beam_dig_non_overlap,2))^2/T;
MSE_dig_comm = norm(Pd_theta - mean(beam_dig_comm,2))^2/T;

%% HBF -------------------------------------------------------------------------------------------
beam_hyb_overlap = zeros(T,K); beam_hyb_non_overlap = zeros(T,K); beam_hyb_propose = zeros(T,K); beam_hyb_comm = zeros(T,K);
% MSE_hyb_overlap = 0; MSE_hyb_non_overlap = 0; MSE_hyb_propose = 0;
% rate_hyb_overlap = 0; rate_hyb_non_overlap = 0; rate_hyb_propose = 0; rate_hyb_propose_rand = 0; rate_comm_hyb = 0;
if HBF == 1
    [F_overlap, W_overlap] = MO_AltMin(Q_dig_overlap, F0, Pt);
    [F_non_overlap, W_non_overlap] = MO_AltMin(Q_dig_non_overlap, F0, Pt);
    [F_propose, W_propose] = MO_AltMin(Q_dig_propose, F0, Pt);
    [F_propose_rand, W_propose_rand] = MO_AltMin(Q_dig_propose_rand, F0, Pt);
    [F_comm, W_comm] = MO_AltMin(Q_dig_comm, F0, Pt);
    
    Q_hyb_overlap = zeros(Nt,M,K);
    Q_hyb_non_overlap = zeros(Nt,M,K);
    Q_hyb_propose = zeros(Nt,M,K);
    Q_hyb_propose_rand = zeros(Nt,M,K);
    Q_hyb_comm = zeros(Nt,M,K);
    %j = 1;
    for k = 1:K
        Q_hyb_overlap(:,:,k) = F_overlap*W_overlap(:,:,k);
        Q_hyb_propose_rand(:,:,k) = F_propose_rand*W_propose_rand(:,:,k);
        Q_hyb_comm(:,:,k) = F_comm * W_comm(:,:,k);
        beam_hyb_overlap(:,k) = real(diag(at(:,:,k)'*Q_hyb_overlap(:,:,k)*Q_hyb_overlap(:,:,k)'*at(:,:,k)))/(Nt*Pt);
        Q_hyb_propose(:,:,k) = F_propose*W_propose(:,:,k);
        if ismember(k,Omg)
            beam_hyb_propose(:,k) = real(diag(at(:,:,k)'*Q_hyb_propose(:,:,k)*Q_hyb_propose(:,:,k)'*at(:,:,k)))/(Nt*Pt);
        end
        beam_hyb_comm(:,k) = real(diag(at(:,:,k)'*Q_hyb_comm(:,:,k)*Q_hyb_comm(:,:,k)'*at(:,:,k)))/(Nt*Pt);
    end
    
    %% non-overlap
    for k = 1:K
        Q_hyb_non_overlap(:,:,k) = F_non_overlap*W_non_overlap(:,:,k);
        if k <= J
            beam_hyb_non_overlap(:,k) = real(diag(at(:,:,k)'*Q_hyb_non_overlap(:,:,k)*Q_hyb_non_overlap(:,:,k)'*at(:,:,k)))/(Nt*Pt);
        end
    end
    
    % Compute rate
    rate_hyb_overlap = compute_rate(h,Nt,M,K,Q_hyb_overlap,error);
    rate_hyb_non_overlap = compute_rate(h,Nt,M,K,Q_hyb_non_overlap,error);
    rate_hyb_propose = compute_rate(h,Nt,M,K,Q_hyb_propose,error);
    rate_hyb_propose_rand = compute_rate(h,Nt,M,K,Q_hyb_propose_rand,error);
    rate_hyb_comm = compute_rate(h,Nt,M,K,Q_hyb_comm,error);
    
    % Compute MSE
    MSE_hyb_overlap = norm(Pd_theta - mean(beam_hyb_overlap,2))^2/T;
    beam_hyb_propose( :, ~any(beam_hyb_propose,1) ) = [];
    MSE_hyb_propose = norm(Pd_theta - mean(beam_hyb_propose,2))^2/T;
    beam_hyb_non_overlap( :, ~any(beam_hyb_non_overlap,1) ) = [];
    MSE_hyb_non_overlap = norm(Pd_theta - mean(beam_hyb_non_overlap,2))^2/T;
    MSE_hyb_comm = norm(Pd_theta - mean(beam_hyb_comm,2))^2/T;
end
end % EOF