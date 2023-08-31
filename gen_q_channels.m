function [Q0_all, R_all] = gen_q_channels(Nt,M,K,H,Pt,init_scheme,n_chans)

Q0_all = zeros(Nt,M,K,n_chans);
R_all = zeros(K,n_chans);

for nn = 1:n_chans
    %nn
    Hnn = H(:,:,:,nn);
    % Obtain q with SCA
    [q0, x0] = initialize(M,Nt,K,Hnn,Pt,init_scheme);
    
    %% run parallel
    [q, rate] = gen_q_subcarrier(Nt,M,K,Hnn,x0,q0,Pt);
    
    %q: Nt*M x K
    Q0_all(:,:,:,nn) = reshape(q,[Nt,M,K]);
    R_all(:,nn) = rate;
end