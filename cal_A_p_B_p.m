function [ A, B ] = cal_A_p_B_p( M, nu, h, zeta, Q_m, xi_b, c_b)

Gamma = zeros(1, M); 
A     = zeros(1, M); 
Gamma(M) = 1;
A(M)     = - Q_m(M) / (2 * xi_b * c_b);

for idx = 1: 1: M - 1
    p = M - idx;
    GammaNumerator   = A(p + 1) * exp( nu(p +1) * h(p + 1) ) * ( zeta(p) * nu(p) - zeta(p + 1) * nu(p + 1) ... 
        + ( zeta(p) * nu(p) + zeta(p + 1) * nu(p + 1) ) * Gamma(p + 1) * exp( - 2 * nu(p + 1) * h(p + 1) ) ) ...
        + zeta(p) * nu(p) * ( Q_m(p + 1) - Q_m(p) ) / ( xi_b * c_b );
    GammaDenominator = A(p + 1) * exp( nu(p +1) * h(p + 1) ) * ( zeta(p) * nu(p) + zeta(p + 1) * nu(p + 1) ... 
        + ( zeta(p) * nu(p) - zeta(p + 1) * nu(p + 1) ) * Gamma(p + 1) * exp( - 2 * nu(p + 1) * h(p + 1) ) ) ...
        + zeta(p) * nu(p) * ( Q_m(p + 1) - Q_m(p) ) / ( xi_b * c_b );
    Gamma(p)         =  GammaNumerator ./ GammaDenominator;

    A(p) =  ( A(p + 1) * exp( nu(p +1) * h(p + 1) ) / ( 2 * zeta(p) * nu(p) ) ) ...
        * ( zeta(p) * nu(p) + zeta(p + 1) * nu(p + 1) + ( zeta(p) * nu(p) - zeta(p + 1) * nu(p + 1) ) ...
        * Gamma(p + 1) * exp( - 2 * nu(p + 1) * h(p + 1) ) ) ...
        + ( Q_m(p + 1) - Q_m(p) ) / ( 2 * xi_b * c_b );
end

end

