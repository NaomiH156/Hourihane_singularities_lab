function w_hat = skew_matrix(w)
%accepts a 3x1 vector and returns the skew symmetrix 3x3 matrix
w_hat = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

