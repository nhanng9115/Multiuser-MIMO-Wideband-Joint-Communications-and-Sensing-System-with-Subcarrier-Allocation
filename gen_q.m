function [Q0_all, R_all] = gen_q(Nt,M,K,H,n_chans,SNR_dB_vec,init_scheme)

Q0_all = zeros(Nt,M,K,n_chans,length(SNR_dB_vec));
R_all = zeros(K,n_chans,length(SNR_dB_vec));

for ss = 1:length(SNR_dB_vec)
    ss
    SNR_dB = SNR_dB_vec(ss);
    Pt = db2pow(SNR_dB);
    [Q0_all(:,:,:,:,ss), R_all(:,:,ss)] = gen_q_channels(Nt,M,K,H,Pt,init_scheme,n_chans);
end
