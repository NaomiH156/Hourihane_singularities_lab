function ew = Rodriguez(w,theta)
%uses Rodriguez' formula to calculate exp(omega-hat * theta)
    w_hat = skew_matrix(w);
    ew = eye(3) + w_hat*sin(theta) + (w_hat^2)*(1-cos(theta));
end

