function [FRF,FBB] = MO_AltMin(Fopt, FRF, Pt)
[Nt, Ns, K] = size(Fopt);
Nrf = size(FRF,2);

%% test new initialization
FF = reshape(Fopt,Nt,Ns*K);
Q = FF*FF';
% Q = Fopt(:,:,K/2)*Fopt(:,:,K/2)';
[U,S,V] = eig(Q);
FRF = exp(1i*angle(U(:,[1:Nrf])));

Nrf = size(FRF,2);
FBB = zeros(Nrf, Ns, K);
count = 0;
if Nt > Nrf % HBF
    y = [];
    dis = [];
    %FRF = exp(1i*unifrnd(0,2*pi,Nt,NRF));
    while (isempty(y) || abs(y(1)-y(2)) > 1e-3) && (count <= 100)
    %for ii = 1:20
        y = [0,0];
        for k = 1:K
            FBB(:,:,k) = pinv(FRF) * Fopt(:,:,k);
            y(1) = y(1) + norm(Fopt(:,:,k) - FRF * FBB(:,:,k),'fro')^2;
        end
        dis = [dis, y(1)];
        [FRF, y(2)] = sig_manif(Fopt, FRF, FBB);
        %abs(y(1)-y(2))
        count = count + 1;
    end
else
    FRF = eye(Nrf);
    FBB = Fopt;
end
% count
% figure(10)
% plot([1:length(dis)], dis/(K/2)); hold on;

for k = 1:K
    FBB(:,:,k) = sqrt(Pt) * FBB(:,:,k) / norm(FRF * FBB(:,:,k),'fro');
%     if abs(norm(FRF * FBB(:,:,k),'fro')^2 - Pt) > 1e-2
%         error('check power constraint !!!!!!!!!!!!')
%     end
end

end