function MI = computeMutualInfo(I, i, j)
% Compute the mutual information between variables i and j wrt to the
% distribution of information matrix I
% We use the approximation of Bianco et. al (2013)

I_joint = computeJointMarginalMatInf(I, i, j);
I_ii = I_joint(1:3,1:3);
I_ij = I_joint(1:3,4:6);
I_jj = I_joint(4:6,4:6);
MI = 0.5 * log ( det(I_ii + eye(3,3)) / ...
                det(I_ii - I_ij * pinv(I_jj) * I_ij' + eye(3,3)));
end

