function rate_mean = compute_rate_ICSI(h,Nt,M,K,Q)

% q: Nt*M x K
% Q: Nt x M

rate = zeros(M,K);

for k = 1:K
    Qk = Q(:,:,k);
    for m = 1:M
        q_mk = Qk(:,m);
        error_mk = eps_mk/sqrt(2) * (randn(1,Nt) + 1i*randn(1,Nt));
        h_mk = h(:,m,k) + error_mk;
        P_S = abs(h_mk'*q_mk)^2;
        P_I = 0;
        for j = 1:M
            if j ~= m
                q_jk = Qk(:,j);
                P_I = P_I + abs(h(:,m,k)'*q_jk)^2;
            end
        end
        SINR_mk = P_S / (P_I + 1);
        rate(m,k) = log2(1 + SINR_mk);
    end
    
end
% rate_mean = mean(mean(rate(:,[Omg(end)+1:end]),2));
rate_mean = sum(mean(rate,2));

end