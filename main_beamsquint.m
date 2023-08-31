% close all
clear all;
clc

% N = 128; K = 3; fc = 28e9; B = 1e9;
N = 128; K = 3; fc = 300e9; B = 60e9;
% N = 32; K = 3; fc = 300e9; B = 30e9;

theta_l = 0.5;
theta = 0.4:0.00001:0.6;

colors = {'b','r','g'};
%% ULA
figure
for k = 1:K
    fk = fc + B*(2*k-1-K)/(2*K);
    kappa = fk/fc;
    
    x = kappa*theta - theta_l;
    
    D = @(x) 1/N * abs((sin(N/2*pi.*x)) ./ (sin(pi/2.*x)));
    
    plot(theta,D(x),'LineWidth',2); hold on
end
% fig = gcf;
% fig.Position(3:4)=[400,300];
xlabel('$\vartheta_{pk} = \frac{f_k}{f_{\mathrm{c}}} \sin(\phi_p)$','Interpreter','latex', 'FontSize', 12)
ylabel('Normalized array gain','FontSize', 12)
% legend('$f_1$','$f_{\mathrm{c}}$','$f_K$','Interpreter','latex')
% title_str = strcat('ULA, $N_r=',num2str(N),'$, $f_c=',num2str(fc/1e9),'$ GHz, $B=',num2str(B/1e9),'$ GHz');
% title(title_str,'Interpreter','latex', 'FontSize', 12);
grid on

% %% UPA
% theta_h = 0.4:0.002:0.6;
% theta_v = 0.4:0.002:0.6;
% 
% figure
% 
% Nh = sqrt(N); Nv = N/Nh;
% 
% for k = 1:K
%     fk = fc + B*(2*k-1-K)/(2*K);
%     kappa = fk/fc;
%     [x_h,x_v] = meshgrid(theta_h , theta_v);
%     
%     D = 1/N * abs((sin(Nh/2*pi*(kappa*x_h - theta_l))) ./ (sin(pi/2.*(kappa*x_h - theta_l)))) .* ...
%         abs((sin(Nv/2*pi.*(kappa*x_v - theta_l))) ./ (sin(pi/2.*(kappa*x_v - theta_l))));
%     
%     surf(x_h,x_v,D, 'FaceColor',colors{k}); hold on
% end
% fig = gcf;
% fig.Position(3:4)=[400,300];
% 
% xlabel('$\theta^{\mathrm{r}}_{qk}$','Interpreter','latex', 'FontSize', 12)
% ylabel('$\phi^{\mathrm{r}}_{qk}$','Interpreter','latex', 'FontSize', 12)
% zlabel('Normalized array gain','FontSize', 12)
% % legend('$f_1$','$f_{\mathrm{c}}$','$f_K$','Interpreter','latex')
% 
% % title_str = strcat('UPA, $N_r=',num2str(N),'$, $f_c=',num2str(fc/1e9),'$ GHz, $B=',num2str(B/1e9),'$ GHz');
% % title(title_str,'Interpreter','latex', 'FontSize', 12);
