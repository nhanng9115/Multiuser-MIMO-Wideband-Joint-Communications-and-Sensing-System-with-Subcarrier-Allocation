function [q0, rate] = gen_q_SCA(Nt,M,h,x0,q0,Pt)

opt_solver = sdpsettings('solver','mosek','verbose',0);

scale_factor = 0.1;
h = scale_factor*h;
sigma2 = scale_factor^2;
% sigma2 = 1;

convergence = 0;
obj_old = 0;
count = 0;
count_max = 100;
while convergence == 0
    
    yalmip('clear')
    
    %% variables
    q = sdpvar(Nt*M,1,'full','complex');
    x = sdpvar(M,1,'full','real');
    
    %% objective function
    obj = sum(log2(1 + x),'all');
    
    %% constraints
    F = [];
    F = [F, x >= 0];
    F = [F, sum(q.^2,1) <= Pt]; % normalize to get SNR as the max power
    % second constrant
    %for k = 1:K
    %F = [F, norm(q(:,k))^2 <= SNR]; % normalize to get SNR as the max power
    for m = 1:M
        hmk = h(:,m);
        Hmk = hmk*hmk';
        
        % construct Hi_hat and Hi_bar
        Hmk_hat = zeros(M*Nt,M*Nt);
        Hmk_hat(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = Hmk;
        
        Hmk_bar = kron(eye(M),Hmk);
        Hmk_bar(Nt*(m-1)+1:Nt*m,Nt*(m-1)+1:Nt*m) = zeros(Nt,Nt);
        
        F1(m) = q'*Hmk_bar*q + sigma2;
        F_qol(m) = q0'*Hmk_hat*q0/x0(m)^2 * x(m) - 2*real(q0'*Hmk_hat*q)/x0(m); % (17)
    end
    %end
    F = [F, F1 + F_qol <= 0]; % (18)
    
    %% Solving
    optimize(F, -obj, opt_solver);
    
    %% update variables
    q0 = double(q);
    x0 = double(x);
    
    %% compute objective value
    obj_new = double(obj);
    if abs(obj_new - obj_old) <= 1e-3 || count >= count_max
        convergence = 1;
    else
        obj_old = obj_new;
        count = count + 1;
    end
%     q0
end
count
rate_all = log2(1 + x0);
rate = sum(rate_all,1); % rate for subcarriers
power = vecnorm(q0).^2;
end