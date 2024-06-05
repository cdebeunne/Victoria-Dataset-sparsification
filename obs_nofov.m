%% Define toy example 

% Transform between cameras
T_c2_c1 = [-1 0 0 -1;
    0 -1 0 0;
    0 0 1 0;
    0 0 0 1];
R_c2_c1 = T_c2_c1(1:3, 1:3);
t_c2_c1 = T_c2_c1(1:3, 4);

% Define two landmarks
l1 = rand(3,1); %[0.5,0.5,1]';
l2 = rand(3,1); %[-0.5,0.5,0]';

% Define a rover motion
yaw = pi / 4;
R_c1p_w = [cos(yaw) sin(yaw) 0;
    -sin(yaw) cos(yaw) 0;
    0 0 1];
t_c1p_w = [1,0,0]';
t_c1p_w = t_c1p_w / norm(t_c1p_w);
T_c1p_w = [R_c1p_w, t_c1p_w;
    0 0 0 1];
T_c1_w = eye(4);

%% Define observation model and jacobians 

obs = @(T, x) ((T(1:3,1:3) * x + T(1:3,4)) / norm(T(1:3,1:3) * x + T(1:3,4)));
J_norm = @(x) (eye(3,3) / norm(x) - x * x' / norm(x)^3);
J_e1p_lambda = @(l, lambda) (J_norm(R_c1p_w * l + lambda * t_c1p_w) * t_c1p_w);
J_e1p_l1 = @(l, lambda) (J_norm(R_c1p_w * l + lambda * t_c1p_w) * R_c1p_w);
J_e2p_lambda = @(l, lambda) (J_norm(R_c2_c1 * R_c1p_w * l + R_c1p_w * t_c2_c1 + lambda * t_c1p_w) * t_c1p_w);
J_e2p_l2 = @(l, lambda) (J_norm(R_c2_c1 * R_c1p_w * l + R_c1p_w * t_c2_c1 + lambda * t_c1p_w) * R_c2_c1 * R_c1p_w);
J_e1_l1 = @(l) (J_norm(l));
J_e2_l2 = @(l) (J_norm(R_c2_c1 * l + t_c2_c1));

ray1 = obs(T_c1_w, l1);
ray1p = obs(T_c1p_w,l1);
ray2 = obs(T_c2_c1 * T_c1_w, l2);
ray2p = obs(T_c2_c1 * T_c1p_w, l2);

%% Compute errors and jacss

lambda_est = 0.8;
T_c1p_w_est = [R_c1p_w, lambda_est * t_c1p_w;
    0 0 0 1];
T_c2p_w_est = T_c2_c1 * T_c1p_w_est;
l1_est = triangulate(ray1, [0,0,0], R_c1p_w'* ray1p, - lambda_est * R_c1p_w' * t_c1p_w);
l2_est = triangulate(R_c2_c1' * ray2, - R_c2_c1' * t_c2_c1, T_c2p_w_est(1:3,1:3)' * ray2p, -T_c2p_w_est(1:3,1:3)' * T_c2p_w_est(1:3, 4));

for i=1:10

    T_c1p_w_est = [R_c1p_w, lambda_est * t_c1p_w;
    0 0 0 1];
    T_c2p_w_est = T_c2_c1 * T_c1p_w_est;

    % errors
    ray2p_est = obs(T_c2_c1 * T_c1p_w_est, l2_est);
    ray1p_est = obs(T_c1p_w_est, l1_est);
    e1 = ray1 - obs(T_c1_w, l1_est);
    e1p = ray1p - ray1p_est;
    e2 = ray2 - obs(T_c2_c1 * T_c1_w, l2_est);
    e2p = ray2p - ray2p_est;
    e = [e1; e1p; e2; e2p];
    
    % order of the FIM lambda, l1, l2
    J_e1 = [zeros(3,1), J_e1_l1(l1), zeros(3,3)];
    J_e1p = [J_e1p_lambda(l1_est, lambda_est), J_e1p_l1(l1_est, lambda_est), zeros(3,3)];
    J_e2 = [zeros(3,1), zeros(3,3), J_e2_l2(l2)];
    J_e2p = [J_e2p_lambda(l2_est, lambda_est), zeros(3,3), J_e2p_l2(l2_est, lambda_est)];
    J = [J_e1; J_e1p; J_e2; J_e2p];
    I = J' * J; 
    
    % GN step
    grad = 2 * e' * J;
    dx = -0.5 * inv(I) * grad';
    
    % Update 
    lambda_est = lambda_est + dx(1);
    l1_est = l1_est + dx(2:4);
    l2_est = l2_est + dx(5:7);
    disp(norm(e))
end



