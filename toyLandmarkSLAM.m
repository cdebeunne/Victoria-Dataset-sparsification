
%% Create a toy example

% We create a toy pose graph to understand basic concept
%
%  l0 <        > l1 <     
%       \    /      |      
%        \  /       | 
%          x0 --> x1 --> x2     
%        /  \      |
%       /    \     | 
%   l3 <      > l2 < 
% The observation model is z = R' * (l - x)

%% Formula of the 2D Jacobian

R = @(x) ([cos(x) -sin(x); sin(x) cos(x)]);
R_d = @(x) ([-sin(x) cos(x); -cos(x) -sin(x)]);
J = @(x, l) ([-R(x(3)), R_d(x(3)) * (l - x(1:2))', R(x(3))]);

dt = 0.1;
J_mot = @(x1,x2,u) ([-1 0 sin(x1(3))*u(1)*dt 1 0 0; ...
    0 -1 cos(x1(3))*u(1)*dt 0 1 0; ...
    0 0 -1 0 0 1]);
f_mot = @(x, u) [R(x(3)) * [u(1) * dt ; 0] + x(1:2)'; ...
    x(3) + u(2) * dt]';

%% Problem statement 
sigma = 0.05;
C = ones(2,2) * sigma^2;

x0 = [rand, rand, rand * pi];
R0 = R(x0(3));

u0 = [0.1, 0.05];
x1 = f_mot(x0, u0);
R1 = R(x1(3));

u1 = [0.1, 0.05];
x2 = f_mot(x1, u1);
R2 = R(x2(3));

l0 = [0, 1];
l1 = [0.5, 0.5];
l2 = [0.5, -0.5];
l3 = [-0.5, -0.5];

z_0_0 = (R0 * (l0 - x0(1:2))')';
z_0_1 = (R0 * (l1 - x0(1:2))')';
z_0_2 = (R0 * (l2 - x0(1:2))')';
z_0_3 = (R0 * (l3 - x0(1:2))')';

z_1_1 = (R0 * (l1 - x1(1:2))')';
z_1_2 = (R0 * (l2 - x1(1:2))')';
z_1_3 = (R0 * (l3 - x1(1:2))')';

z_2_2 = (R0 * (l2 - x2(1:2))')';
z_2_3 = (R0 * (l3 - x2(1:2))')';


% Information Matrix for 2 pose and 3 landmarks
J_0_0 = J(x0, l0);
J_0_0_aug = [J_0_0(:,1:3), zeros(2, 6) ,J_0_0(:,4:5), zeros(2,6)];
J_0_1 = J(x0, l1);
J_0_1_aug = [J_0_1(:,1:3), zeros(2, 6) ,zeros(2,2), J_0_1(:,4:5), zeros(2,4)];
J_0_2 = J(x0, l2);
J_0_2_aug = [J_0_2(:,1:3), zeros(2, 6) ,zeros(2,4), J_0_1(:,4:5), zeros(2,2)];

J_1_1 = J(x1, l1);
J_1_1_aug = [zeros(2, 3), J_1_1(:,1:3), zeros(2, 3), J_1_1(:,4:5), zeros(2,6)];
J_1_2 = J(x1, l2);
J_1_2_aug = [zeros(2, 3), J_1_2(:,1:3), zeros(2, 3), zeros(2,2) ,J_1_2(:,4:5), zeros(2,4)];
J_1_3 = J(x1, l3);
J_1_3_aug = [zeros(2, 3), J_1_3(:,1:3), zeros(2, 3), zeros(2,6), J_1_3(:,4:5)];

J_2_2 = J(x2, l2);
J_2_2_aug = [zeros(2, 6), J_2_2(:,1:3), zeros(2,2) ,J_2_2(:,4:5), zeros(2,4)];
J_2_3 = J(x2, l3);
J_2_3_aug = [zeros(2, 6), J_2_3(:,1:3), zeros(2,6), J_2_3(:,4:5)];

C_mot = 5 * ones(3);
J_mot_0_1 = [J_mot(x0, x1, u0), zeros(3,3), zeros(3, 8)];
J_mot_1_2 = [zeros(3,3), J_mot(x1, x2, u1), zeros(3, 8)];

I = J_0_0_aug' * C * J_0_0_aug + J_0_1_aug' * C * J_0_1_aug +...
    J_0_2_aug' * C * J_0_2_aug +...
    J_1_1_aug' * C * J_1_1_aug + J_1_2_aug' * C * J_1_2_aug +...
    J_1_3_aug' * C * J_1_3_aug +...
    J_2_2_aug' * C * J_2_2_aug + J_2_3_aug' * C * J_2_3_aug +...
    J_mot_0_1' * C_mot * J_mot_0_1 + J_mot_1_2' * C_mot * J_mot_1_2;

Jac = [
    J_mot_0_1;
    J_mot_1_2;
    J_0_0_aug;
    J_0_1_aug;
    J_0_2_aug;
    J_1_1_aug;
    J_1_2_aug;
    J_1_3_aug;
    J_2_2_aug;
    J_2_3_aug];


x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 600; % Width of figure
height = 400; % Height of figure (by default in pixels)
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
figure('position', [x y width height]);
H = abs(I)>0;
% H(17, 4) = 1;
% H(17, 7) = 1;
imagesc(H, [0, 1]);
set(gca,'xaxisLocation','top')
hold on
g_y=0.5:17.5; % user defined grid Y [start:spaces:end]
g_x=0.5:17.5; % user defined grid X [start:spaces:end]
for i=1:length(g_x)
   plot([g_x(i) g_x(i)],[g_y(1) g_y(end)], 'Color', [0,0,0, 0.5]) %y grid lines
   hold on    
end
for i=1:length(g_y)
   plot([g_x(1) g_x(end)],[g_y(i) g_y(i)], 'Color', [0,0,0, 0.5]) %x grid lines
   hold on    
end
% colorMap = [linspace(1,0,256)', ones(256,2)];
colorMap = [linspace(1,0,256)', linspace(1,0,256)', ones(256,1)];
colormap(colorMap);
xticks([2 5 8, 10.5, 12.5, 14.5, 16.5])
xticklabels({'$\mathbf{x}_0$','$\mathbf{x}_1$','$\mathbf{x}_2$', '$\mathbf{l}_0$', '$\mathbf{l}_1$', '$\mathbf{l}_2$', '$\mathbf{l}_3$'})
yticks([2 5 8, 10.5, 12.5, 14.5, 16.5])
yticklabels({'$\mathbf{x}_0$','$\mathbf{x}_1$','$\mathbf{x}_2$', '$\mathbf{l}_0$', '$\mathbf{l}_1$', '$\mathbf{l}_2$', '$\mathbf{l}_3$'})
axis equal tight

 
figure('position', [x y width height * 1.222222]);
Jacs = abs(Jac)>0;
imagesc(Jacs, [0, 1]);
set(gca,'xaxisLocation','top')
hold on
g_y=0.5:22.5; % user defined grid Y [start:spaces:end]
g_x=0.5:17.5; % user defined grid X [start:spaces:end]
for i=1:length(g_x)
   plot([g_x(i) g_x(i)],[g_y(1) g_y(end)], 'Color', [0,0,0, 0.5]) %y grid lines
   hold on    
end
for i=1:length(g_y)
   plot([g_x(1) g_x(end)],[g_y(i) g_y(i)], 'Color', [0,0,0, 0.5]) %x grid lines
   hold on    
end
% colorMap = [linspace(1,0,256)', ones(256,2)];
colorMap = [linspace(1,0,256)', linspace(1,0,256)', ones(256,1)];
colormap(colorMap);
xticks([2 5 8, 10.5, 12.5, 14.5, 16.5])
xticklabels({'$\mathbf{x}_0$','$\mathbf{x}_1$','$\mathbf{x}_2$', '$\mathbf{l}_0$', '$\mathbf{l}_1$', '$\mathbf{l}_2$', '$\mathbf{l}_3$'})
yticks([]);
axis equal tight

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
cmap = jet;
imagesc(log10(abs(I_marg)), [-21, -2]);
colormap(cmap)
yticks({})
xticks({})
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf, 'dense.svg');
figure;
imagesc(log10(abs(I_spars)), [-21, -2]);
colormap(cmap)
yticks({})
xticks({})
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf, 'unary.svg');
figure;
imagesc(log10(abs(I_spars2)), [-21, -2]);
colormap(cmap)
yticks({})
xticks({})
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf, 'chain.svg');
figure;
imagesc(log10(abs(I)), [-21, -2]);
yticks({})
xticks({})
colormap(cmap)
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf, 'full.svg');





