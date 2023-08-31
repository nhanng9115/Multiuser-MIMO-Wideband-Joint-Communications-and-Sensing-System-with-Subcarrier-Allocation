function [q, rate] = gen_q_subcarrier(Nt,M,K,h,x0,q0,Pt)
q = zeros(Nt*M,K);
rate = zeros(K,1);
for k = 1:K
    [q(:,k), rate(k)] = gen_q_SCA(Nt,M,h(:,:,k),x0(:,k),q0(:,k),Pt);
end

end