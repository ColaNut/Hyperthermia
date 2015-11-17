function [ k ] = cal_k( N, epsilon, sigma, mu, omega_0 )

k_1 = zeros(N + 2, 1);
k_2 = zeros(N + 2, 1);
k   = zeros(N + 2, 1);

%  some thing is wrong in the derivation
k_1 = omega_0 .* sqrt( (mu .* epsilon) / 2 ) .* sqrt( sqrt( 1 + (sigma ./ (omega_0 * epsilon)).^2 ) + 1);

k_2 = omega_0 .* sqrt( (mu .* epsilon) / 2 ) .* sqrt( sqrt( 1 + (sigma ./ (omega_0 * epsilon)).^2 ) - 1);

k = k_1 - j * k_2;

end