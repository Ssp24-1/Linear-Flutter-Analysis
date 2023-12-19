function [Ms, Ks] = GetStructuralMatrices(ms, wh, w_theta, w_beta, a, b, c, x_theta, x_beta, r_theta2, r_beta2)

Ms = ms*b*b*[1  x_theta x_beta;
    x_theta r_theta2    (r_beta2 + x_beta*(c-a));
    x_beta  (r_beta2 + x_beta(c-a))   r_beta2];

Ks = ms*b*b*[wh^2   0   0;
    0   (r_theta2*w_theta*w_theta)    0;
    0   0   (r_beta2*w_beta*w_beta)];

