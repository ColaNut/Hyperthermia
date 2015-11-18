function [ A, B ] = cal_A_p_B_p( M, nu, h, zeta, Q_m, xi_b, c_b)

A       = zeros(1, M); 
B       = zeros(1, M); 
A(M)    = - Q_m(M) / (2 * xi_b * c_b);
B(M)    = A(M);

for idx = 1: 1: M - 1
    p = M - idx;
    A(p) = A(p + 1) * exp( nu(p +1) * h(p + 1) ) * ( 0.5 + ( zeta(p + 1) * nu(p + 1) ) / ( 2 * zeta(p) * nu(p) ) ) ... 
         + B(p + 1) * exp( (- 1) *  nu(p + 1) * h(p + 1) ) * ( 0.5 -  ( zeta(p + 1) * nu(p + 1) ) / ( 2 * zeta(p) * nu(p) ) ) ...
         + ( Q_m(p + 1) - Q_m(p) ) / ( 2 * xi_b * c_b );
    B(p) = A(p + 1) * exp( nu(p +1) * h(p + 1) ) * ( 0.5 - ( zeta(p + 1) * nu(p + 1) ) / ( 2 * zeta(p) * nu(p) ) ) ... 
         + B(p + 1) * exp( (- 1) *  nu(p + 1) * h(p + 1) ) * ( 0.5 +  ( zeta(p + 1) * nu(p + 1) ) / ( 2 * zeta(p) * nu(p) ) ) ...
         + ( Q_m(p + 1) - Q_m(p) ) / ( 2 * xi_b * c_b );
end

end

