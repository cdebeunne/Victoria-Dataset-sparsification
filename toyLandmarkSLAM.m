
%% Create a toy example

% We create a toy pose graph to understand basic concept
%
%  l0 <        > l1      
%       \    /           
%        \  /       
%          x0     
%        /  \     
%       /    \     
%   l3 <      > l2  
% The observation model is z = R' * (l - x)

%% Problem statement 
sigma = 0.05;
C = [sigma^2 0;
    0 sigma^2];

z_0_0 = [0, 1];
z_0_1 = [0.5, 0.5];
z_0_2 = [0.5, -0.5];
z_0_3 = [-0.5, -0.5];

x0 = [0, 0, pi/4];
l0 = [0, 1];
l1 = [0.5, 0.5];
l2 = [0.5, -0.5];
l3 = [-0.5, -0.5];

% Formula of the 2D Jacobian
R = @(x) ([cos(x) -sin(x); sin(x) cos(x)]);
R_d = @(x) ([-sin(x) cos(x); -cos(x) -sin(x)]);
J = @(x, l) ([-R(x(3))', R_d(x(3)) * (l - x(1:2))', R(x(3))]);

% Information Matrix for 2 pose and 3 landmarks
J_0_0 = J(x0, l0);
J_0_0_aug = [J_0_0(:,1:3), J_0_0(:,4:5), zeros(2,6)];
J_0_1 = J(x0, l1);
J_0_1_aug = [J_0_1(:,1:3), zeros(2,2), J_0_1(:,4:5), zeros(2,4)];
J_0_2 = J(x0, l2);
J_0_2_aug = [J_0_2(:,1:3), zeros(2,4), J_0_1(:,4:5), zeros(2,2)];
J_0_3 = J(x0, l3);
J_0_3_aug = [J_0_3(:,1:3), zeros(2,6), J_0_1(:,4:5)];

I = J_0_0_aug' * C * J_0_0_aug + J_0_1_aug' * C * J_0_1_aug +...
    J_0_2_aug' * C * J_0_2_aug + J_0_3_aug' * C * J_0_3_aug;

%% Marginalization

% We remove x0 and l0
idx_to_remove = 1:3;

% We keep l0, l1, l2, l3
idx_to_keep = 4:11;

I_rr = I(idx_to_remove, idx_to_remove);
I_kk = I(idx_to_keep, idx_to_keep);
I_kr = I(idx_to_keep, idx_to_remove);

I_marg = I_kk + I_kr * inv(I_rr) * I_kr';

%% Sparsification unary factor

% Compute covariance matrix
Sigma_marg = inv(I_marg);

% Compute information matrix of each unary factor (Mazuran et al.)
I_l0 = inv(Sigma_marg(1:2, 1:2));
I_l1 = inv(Sigma_marg(3:4, 3:4));
I_l2 = inv(Sigma_marg(5:6, 5:6));
I_l3 = inv(Sigma_marg(7:8, 7:8));


% Build the global info mat
I_spars = zeros(8);
I_spars(1:2, 1:2) = I_l0;
I_spars(3:4, 3:4) = I_l1;
I_spars(5:6, 5:6) = I_l2;
I_spars(7:8, 7:8) = I_l3;

d_kl = 0.5 * (trace(I_spars * Sigma_marg) - log(det(I_spars * Sigma_marg)) - 8)

%% Sparsification lmk chain

% Compute covariance matrix
Sigma_marg = inv(I_marg);

% Compute information matrix of each unary factor (Mazuran et al.)
J_l0 = [eye(2,2), zeros(2,6)];
I_l0 = inv(Sigma_marg(1:2, 1:2));

J_l0_l1 = [eye(2,2), -eye(2,2), zeros(2,4)];
I_l0_l1 = inv(J_l0_l1 * Sigma_marg * J_l0_l1');

J_l1_l2 = [zeros(2,2), eye(2,2), -eye(2,2), zeros(2,2)];
I_l1_l2 = inv(J_l1_l2 * Sigma_marg * J_l1_l2');

J_l2_l3 = [zeros(2,4), eye(2,2), -eye(2,2)];
I_l2_l3 = inv(J_l2_l3 * Sigma_marg * J_l2_l3');

% Build the global info mat
J = [J_l0; J_l0_l1; J_l1_l2; J_l2_l3];
Omega = [I_l0, zeros(2,6);
    zeros(2,2), I_l0_l1, zeros(2,4);
    zeros(2,4), I_l1_l2, zeros(2,2);
    zeros(2,6), I_l2_l3];
I_spars2 = J' * Omega * J

d_kl = 0.5 * (trace(I_spars2 * Sigma_marg) - log(det(I_spars2 * Sigma_marg)) - 8);


figure;
cmap = [0.6 0.6 0.6;turbo];
imagesc(log10(abs(I_marg)), [-15, -2]);
colormap(cmap)
yticks({})
xticks({})
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf, 'dense.svg');
figure;
imagesc(log10(abs(I_spars)), [-15, -2]);
colormap(cmap)
yticks({})
xticks({})
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf, 'unary.svg');
figure;
imagesc(log10(abs(I_spars2)), [-15, -2]);
colormap(cmap)
yticks({})
xticks({})
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf, 'chain.svg');
figure;
imagesc(log10(abs(I)), [-15, -2]);
yticks({})
xticks({})
colormap(cmap)
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf, 'full.svg');





