%% Create a toy example

% We create a toy pose graph to understand basic concept
%
%  l0 <        > l1      
%       \    /    ^       
%        \  /     |  
%          x0     x1
%           \     |
%            \    ˇ 
%              > l2  
% The observation model is z = R' * (l - x)

sigma = 0.05;
C = [sigma^2 0;
    0 sigma^2];

z_0_0 = [0, 1];
z_0_1 = [0.5, 0.5];
z_0_2 = [0.5, -0.5];

z_1_1 = [-0.5, 0.5];
z_1_2 = [-0.5, -0.5];

l0 = [0, 1];
l1 = [0.5, 0.5];
l2 = [0.5, -0.5];

x0 = [0, 0, 0];
x1 = [1, 0, 0];

% Formula of the 2D Jacobian
R = @(x) ([cos(x) -sin(x); sin(x) cos(x)]);
R_d = @(x) ([-sin(x) cos(x); -cos(x) -sin(x)]);
J = @(x, l) ([-R(x(3))', R_d(x(3)) * (l - x(1:2))', R(x(3))]);


% Information Matrix for 2 pose and 3 landmarks
J_0_0 = J(x0, l0);
J_0_0_aug = [J_0_0(:,1:3), zeros(2,3), J_0_0(:,4:5), zeros(2,4)];
J_0_1 = J(x0, l1);
J_0_1_aug = [J_0_1(:,1:3), zeros(2,5), J_0_1(:,4:5), zeros(2,2)];
J_0_2 = J(x0, l2);
J_0_2_aug = [J_0_2(:,1:3), zeros(2,7), J_0_1(:,4:5)];
J_1_1 = J(x1, l1);
J_1_1_aug = [zeros(2,3), J_1_1(:,1:3), zeros(2,2), J_1_1(:,4:5), zeros(2,2)];
J_1_2 = J(x1, l2);
J_1_2_aug = [zeros(2,3), J_1_2(:,1:3), zeros(2,4), J_1_2(:,4:5)];

I = J_0_0_aug' * C * J_0_0_aug + J_0_1_aug' * C * J_0_1_aug +...
    J_0_2_aug' * C * J_0_2_aug + J_1_1_aug' * C * J_1_1_aug +...
    J_1_2_aug' * C * J_1_2_aug;

% We remove x0 and l0
idx_to_remove = [1:3, 7:8];

% We keep x1, l1 and l2
idx_to_keep = [4:6, 9:12];

I_rr = I(idx_to_remove, idx_to_remove);
I_kk = I(idx_to_keep, idx_to_keep);
I_kr = I(idx_to_keep, idx_to_remove);

I_marg = I_kk + I_kr * inv(I_rr) * I_kr'

