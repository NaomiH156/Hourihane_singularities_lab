function q = adjoint_inv(g)
R = g(1:3,1:3);
p = g(1:3,4);
p_hat = skew_matrix(p);

q(1:3,1:3) = R.';
q(1:3,4:6) = -1*skew_matrix((p_hat*R))*R.';
q(4:6,4:6) = R.';
end

