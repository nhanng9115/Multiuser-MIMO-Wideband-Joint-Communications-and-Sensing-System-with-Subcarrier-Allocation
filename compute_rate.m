function rate_mean = compute_rate(h,Nt,M,K,Q,e)

% q: Nt*M x K
% Q: Nt x M
CSI_error = 1;
rate = zeros(M,K);

for k = 1:K
    Qk = Q(:,:,k);
    for m = 1:M
        q_mk = Qk(:,m);
        h_mk = h(:,m,k);
        
        if CSI_error == 1
            eps_mk = e;
            error_mk = eps_mk/sqrt(2) * (randn(Nt,1) + 1i*randn(Nt,1));
            h_mk = h_mk - error_mk;
        end
        
        P_S = abs(h_mk'*q_mk)^2;
        P_I = 0;
        for j = 1:M
            if j ~= m
                q_jk = Qk(:,j);
                P_I = P_I + abs(h_mk'*q_jk)^2;
            end
        end
        SINR_mk = P_S / (P_I + 1);
        rate(m,k) = log2(1 + SINR_mk);
    end
    
end
% rate_mean = mean(mean(rate(:,[Omg(end)+1:end]),2));
rate_mean = sum(mean(rate,2));

end