function Q = opt_radcom(Nt,n_MS,weight_sens,Cbar,Q0,Pt)

rho = weight_sens;
rho1 = (1-rho);
X0 = Q0;
C0 = Cbar;
% Create the problem structure.
MO = obliquecomplexfactory(Nt,n_MS);
problem.M = MO;

% Define the problem cost function and its gradient.
problem.cost = @cost;
    function f = cost(X)
        f = rho*norm(X*X' - C0,'fro')^2 + rho1*norm(X - X0,'fro')^2;
    end
problem.grad = @(X) problem.M.egrad2rgrad(X,egrad(X));
    function g = egrad(X)
        g = 4*rho*(X*X' - C0)*X + 2*rho1*(X - X0);
    end

% checkgradient(problem);
% warning('off', 'manopt:getHessian:approx');

% Execute the optimization
Qtmp = conjugategradient(problem);
Q = sqrt(Pt/n_MS)*Qtmp;

% check power constraint
if abs(norm(Q,'fro')^2 - Pt) > 1e-4
    error('check power constraint!!!!!!!!!');
end

%test
% term1 = rho*norm(Q*Q' - Cbar,'fro')^2;
% term2 = rho1*norm(Q - X1,'fro')^2;
% ratio = term1 / term2;

end