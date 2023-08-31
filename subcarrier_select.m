function [Omg, Omg_rand] = subcarrier_select(K,J,Q0,C0)

norm_list = zeros(K,1);
for k = 1:K
    norm_list(k) = norm(Q0(:,:,k)*Q0(:,:,k)' - C0(:,:,k),'fro')^2;
end
[~, Omg] = mink(norm_list,J); 
Omg_rand = randperm(K,J); %Omg_not = setdiff([1:K],Omg);

end % EOF