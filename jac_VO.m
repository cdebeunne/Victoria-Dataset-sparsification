J_0_1 = [-50 0 -100 0 -50 0 0 -50 0;
         0 125 0 -25 0 -50 -25 0 -50];

J_0_1_aug = [J_0_1(:,1:6), zeros(2,3), J_0_1(:,7:9), zeros(2,3)];

J_0_2 = [50 0 -100 0 -50 0 0 -50 0;
         0 125 0 25 0 -50 25 0 -50];

J_0_2_aug = [J_0_2(:,1:6), zeros(2,6), J_0_2(:,7:9)];

J_0_0 = [25  0   -100   0   -50   0     0   -50  0;
         0 106.25  0   12.5  0   -50  12.5   0  -50];

J_0_0_aug = [J_0_0, zeros(2,6)];

I = J_0_0_aug' * J_0_0_aug + J_0_1_aug'*J_0_1_aug + J_0_2_aug'*J_0_2_aug;
I = eye(15) + I;

I_kr = I(10:end, 1:9);
I_rk = I(1:9, 10:end);
I_rr = I(1:9, 1:9);
I_kk = I(10:end, 10:end);

[V, D] = svd(I_rr);
S_inv = zeros(9,9);
S_inv(1:6,1:6) = inv(D(1:6,1:6)); 
Sigma_rr = V * S_inv * V';
Sigma_rr = inv(I_rr);

I_k = I_kk - I_kr * Sigma_rr * I_kr'

J = [eye(3,3) zeros(3,6);
     eye(3,3) -eye(3,3) zeros(3,3);
     zeros(3,3) eye(3,3) -eye(3,3)]