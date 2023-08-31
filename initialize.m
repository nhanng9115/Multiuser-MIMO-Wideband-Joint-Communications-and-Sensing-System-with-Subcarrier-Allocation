function [q0, x0] = initialize(M,Nt,K,h,Pt,init_scheme)

if init_scheme == 0
    %% Random initial point ----------------------------------------
    n_rand = 1;
    x0_tmp = zeros(M,K,n_rand);
    q0_tmp = sqrt(Pt/2) * (randn(M*Nt,K,n_rand)+ 1j*randn(M*Nt,K,n_rand));
    sum_rate = zeros(M,K,n_rand);
    for ii = 1:n_rand
        for k = 1:K
            for m = 1:M
                hmk = h(:,m,k);
                Hmk = hmk*hmk';
                
                % construct Hi_hat and Hi_bar
                Hmk_hat = zeros(M*Nt,M*Nt);
                Hmk_hat(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = Hmk;
                
                Hmk_bar = kron(eye(M),Hmk);
                Hmk_bar(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = zeros(Nt,Nt);
                
                x0_tmp(m,k,ii) = real(q0_tmp(:,k,ii)'*Hmk_hat*q0_tmp(:,k,ii) / ...
                    (q0_tmp(:,k,ii)'*Hmk_bar*q0_tmp(:,k,ii) + 1));
                sum_rate(m,k,ii) = log2(1 + x0_tmp(m,k,ii));
            end
        end
    end
    sum_rate_all = sum(sum(sum_rate,1),2);
    [rate_rand,ii_max] = max(sum_rate_all(:));
    q0 = q0_tmp(:,:,ii_max);
    x0 = x0_tmp(:,:,ii_max);
else
    %% ZF initialization
    Q_ZF = zeros(Nt,M,K);
    for k = 1:K
        Hk = h(:,:,k)';
        Q_ZF_tmp = Hk'*(Hk*Hk')^(-1);
        channel_power = vecnorm(Q_ZF_tmp).^2;
        pp = waterfilling(Pt,M,channel_power);
        pp = pp ./ vecnorm(Q_ZF_tmp).^2;
        %sum(pp)
        P = diag(sqrt(pp));
        Q_ZF(:,:,k) = Q_ZF_tmp * P;
        % check power
        if abs(norm(Q_ZF(:,:,k),'fro')^2 - Pt) >= 1e-4
            error('wrong power')
        end
    end
    q_ZF = reshape(Q_ZF,[M*Nt,K]);
    
    % test rate
    rate_ZF = zeros(M,K);
    x0_ZF = zeros(M,K);
    for k = 1:K
        for m = 1:M
            hmk = h(:,m,k);
            Hmk = hmk*hmk';
            
            % construct Hi_hat and Hi_bar
            Hmk_hat = zeros(M*Nt,M*Nt);
            Hmk_hat(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = Hmk;
            
            Hmk_bar = kron(eye(M),Hmk);
            Hmk_bar(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = zeros(Nt,Nt);
            
            x0_ZF(m,k) = real(q_ZF(:,k)'*Hmk_hat*q_ZF(:,k) / ...
                (q_ZF(:,k)'*Hmk_bar*q_ZF(:,k) + 1));
            rate_ZF(m,k) = log2(1 + x0_ZF(m,k));
        end
    end
    sum_rate_ZF = sum(sum(rate_ZF,1),2);
    q0 = q_ZF;
    x0 = x0_ZF;
end
end

