clc; clear;
digits;

% parameters
% 20.58 / 14 = 1.47 (divide lung into 14 layer)
M               = 18; % [ skin, muscle, lung * 14, muscle, skin ] % we assume there is air outside the skin (may be water).
rho             = [1.05, 1.06, 0.27, 0.25, 0.25, 0.25, 0.25, 0.26, 0.28, 0.30, 0.32, 0.34, 0.37, 0.39, 0.40, 0.41, 1.06, 1.05]' .* 1000; % g/cm^3 -> kg/m^3
LungLayerThick  = repmat( [1.47], 1, 14 );
LungLayerZeta   = repmat( [0.0143], 1, 14 );
LungLayerCapa   = repmat( [3886], 1, 14 );
LungLayerQ_m    = repmat( [4200], 1, 14 );
thickness       = [0.02, 2.09, LungLayerThick, 2.67, 0.19]' ./ 100; % cm -> m,
zeta            = [0.31, 0.1314, LungLayerZeta, 0.1314, 0.31]';
capa            = [3150, 3396, LungLayerCapa, 3396, 3150]';
Q_m             = [1125, 4200, LungLayerQ_m, 4200, 1125]';
xi_b            = 1.59;
c_b             = 3600;
nu              = sqrt(xi_b * c_b ./ zeta);
T_0             = 38;
% calculate the accumulated depth
AccuDepth           = zeros(M, 1);
% the accumulated depth is used for idx = 2: N + 1;
for idx = 1: 1: M
    AccuDepth(idx:M) = AccuDepth(idx:M) - thickness(idx);
end

Gamma = zeros(1, M); 
A     = zeros(1, M); 

[Gamma, A] = cal_Gamma_p_AND_A_p( M, nu, thickness, zeta, Q_m, xi_b, c_b);

Depth       = (-25.55) / 100;   % 30   cm
Interval    = 0.01  / 100;   % 0.01 cm
Height      = 0     / 100;   % 1    cm

% calculate the rho and sigma array
% the value of rho, and sigma at -30 cm is set to be that of 'M_SAR'-th layer
A_array         = [A(M)];
nu_array        = [nu(M)];
AccuDepth_array = [AccuDepth(M)];
Gamma_array     = [Gamma(M)];
Q_m_array       = [Q_m(M)];
for idx = 1: 1: M
    p   = M + 1 - idx;
    A_array = [A_array, repmat( A(p), 1, int32( thickness(p) / Interval) )];
    nu_array = [nu_array, repmat( nu(p), 1, int32( thickness(p) / Interval) )];
    AccuDepth_array = [AccuDepth_array, repmat( AccuDepth(p), 1, int32( thickness(p) / Interval) )];
    Gamma_array = [Gamma_array, repmat( Gamma(p), 1, int32( thickness(p) / Interval) )];
    Q_m_array = [Q_m_array, repmat( Q_m(p), 1, int32( thickness(p) / Interval) )];
    % if idx == M
    %     A           = A_tmp';
    %     nu          = nu_tmp';
    %     AccuDepth   = AccuDepth_tmp';
    %     Gamma       = Gamma_tmp';
    %     Q_m         = Q_m_tmp';
    % end
end

T_p = zeros( (Height - Depth) / Interval + 1, 1);
z = [Depth : Interval : Height]';
T_p = A_array' .* exp( nu_array' .* (z - AccuDepth_array') ) .* ( 1 + Gamma_array' .* exp( (- 2) * nu_array' .* (z - AccuDepth_array') ) ) + T_0 + Q_m_array' / (xi_b * c_b);

% plot T_p
clf;
set(figure(1),'defaulttextinterpreter','latex');
plot( z, T_p, 'LineWidth',2 );
axis([-0.2555, 0, 0, 40]);
set(gca,'fontsize',14);
xlabel('z/m','FontSize',18);
ylabel('T/($^\circ$C)','FontSize',18);

% add vertical line
hold on
Max = 40;
for idx = 1: 1: 18
    x = AccuDepth(idx) * ones(1, length(z));
    y = 0: Max / (length(z) - 1) : Max;
    plot(x, y, 'r--','LineWidth',1 );
end
hold off