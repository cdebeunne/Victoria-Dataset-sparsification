function I_joint = computeJointMarginalMatInf(I, i, j)
% Compute the joint Marginal information matrix from an information matrix
% and two nodes using regularization from Bianco et. al (2013)

n = length(I)/3;

% reorder the information matrix

% put i in first position
k=1;
% exchange rows
rows_k = I((k-1)*3+1:k*3, :);
rows_i = I((i-1)*3+1:i*3, :);
I((k-1)*3+1:k*3, :) = rows_i;
I((i-1)*3+1:i*3, :) = rows_k;
% exchange columns
cols_k = I(:,(k-1)*3+1:k*3);
cols_i = I(:, (i-1)*3+1:i*3);
I(:,(k-1)*3+1:k*3) = cols_i;
I(:, (i-1)*3+1:i*3) = cols_k;

% put j in second position
k=2;
% exchange rows
rows_k = I((k-1)*3+1:k*3, :);
rows_j = I((j-1)*3+1:j*3, :);
I((k-1)*3+1:k*3, :) = rows_j;
I((j-1)*3+1:j*3, :) = rows_k;
% exchange columns
cols_k = I(:,(k-1)*3+1:k*3);
cols_j = I(:, (j-1)*3+1:j*3);
I(:,(k-1)*3+1:k*3) = cols_j;
I(:, (j-1)*3+1:j*3) = cols_k;

% compute schur complement
I_aa = I(1:3*2, 1:3*2);
I_bb = I(3*2+1:3*n, 3*2+1:3*n);
I_ab = I(1:3*2,3*2+1:3*n);

% Here we use pseudo inverse as I_bb can be low rank
I_joint = I_aa - I_ab * pinv(I_bb) * I_ab';

end

