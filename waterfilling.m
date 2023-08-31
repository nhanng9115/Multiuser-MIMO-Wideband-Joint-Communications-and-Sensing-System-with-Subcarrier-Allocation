function p = waterfilling(P0,M,channel_power)

% Bisection search for mu
mu_low = 0; % Initial low
mu_high = M/(P0 + sum(channel_power)); % Initial high

stop_threshold = 1e-10; % Stop threshold

% Iterate while low/high bounds are further than stop_threshold
while(abs(mu_low - mu_high) > stop_threshold)
    mu = (mu_low + mu_high) / 2; % Test value in the middle of low/high

    % Solve the power allocation
    p = 1/mu - channel_power; 
    p = max(p,0); % Consider only positive power allocation
    
    % Test sum-power constraints
    if (sum(p) > P0) % Exceeds power limit => lower the upper bound
        mu_low = mu;
    else % Less than power limit => increase the lower bound
        mu_high = mu;
    end
    %sum(p)
end
%P0