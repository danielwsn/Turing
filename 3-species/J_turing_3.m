function Ju = J_turing_3(u, y, opt)

%% prepare data and parameters
N = opt.N;
h = opt.L/N;
u1 = u(1:N);
u2 = u(N+1:2*N);
% u3 = u(2*N+1:end);
y1 = y(1:N);
y2 = y(N+1:2*N);
y3 = y(2*N+1:end);

%% calculate
Jy1 = (-opt.tk12 * u2.^2 -opt.k12 -opt.k13).* y1 + ...
      (-2*opt.tk12*u1.*u2+3*opt.tk21*u2.^2+opt.k21).*y2 + ...
      opt.k31*y3;
Jy3 = opt.k13*y1 + opt.k23*y2-(opt.k31+opt.k32)*y3;
Jy2 = -Jy1-Jy3+opt.D2 * My_Laplace(y2, h);
Jy1 = Jy1 + opt.D1 * My_Laplace(y1, h);
Jy3 = Jy3 + opt.D3 * My_Laplace(y3, h);

Ju = [Jy1; Jy2; Jy3];

end

function F = My_Laplace(u, h)
du = diff(u);
F = ([du; 0] - [0; du]) / h^2;
end