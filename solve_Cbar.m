function R = solve_Cbar(Pd_theta,N,a,theta,Pt)
    
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
cvx_solver sedumi
cvx_begin quiet
variable bt
variable R(N,N) hermitian semidefinite
expression u(length(theta))
for i=1:length(theta)
    u(i)=(bt*Pd_theta(i) - a(:,i)'*R*a(:,i));
end
minimize norm(u,2)
subject to
diag(R)==ones(N,1)*Pt/N;
% trace(R)==power;
cvx_end

% bt
end
