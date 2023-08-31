function a = array_response_UPA(Nv,Nh,N,fk,fc,phi,theta)

% Nh = 2; Nv = N/Nh;
% for m= 0:Nh-1
%     for n= 0:Nv-1
%         y(m*Nv+n+1) = exp( 1i* pi* fk/fc * ( m*sin(a1)*sin(a2) + n*cos(a2) ) );
%     end
% end
% y = y.'/sqrt(N);

% Generate array response vectors
% a = zeros(N,1);
% indice_h = zeros(Nh,1);
% indice_v = zeros(1,Nv);

indice_h = 0:Nh-1;
indice_v = 0:Nv-1;

angle_sinsin = indice_h.' .* sin(phi)*sin(theta);
angle_cos = indice_v .* cos(theta);

angle_sinsin_grid = repmat(angle_sinsin, 1, Nv);
angle_cos_grid = repmat(angle_cos, Nh, 1);
angle_grid = angle_sinsin_grid + angle_cos_grid;
angle_vec = angle_grid(:);
% print(shape((1/sqrt(N))*exp(1j*pi*fk/fc*angle_vec)))
% print(shape(a))
a = (1/sqrt(N))*exp(1j*pi*fk/fc*angle_vec);
end