function Jtu = Jt_turing_3(u, y, opt)

%% prepare data and parameters
N = opt.N;
h = opt.L/N;
u1 = u(1:N);
u2 = u(N+1:2*N);
% u3 = u(2*N+1:end);
y1 = y(1:N);
y2 = y(N+1:2*N);
y3 = y(2*N+1:end);

a11 = -opt.tk12 * u2.^2 -opt.k12 -opt.k13;
a12 = -2*opt.tk12*u1.*u2+3*opt.tk21*u2.^2+opt.k21;

%% calculate
Jty1 = opt.D1 * My_Laplace(y1, h) + ...
       a11.* y1 + (-a11-opt.k13).*y2 + opt.k13*y3;
Jty2 = opt.D2 * My_Laplace(y2, h) + ...
       a12.* y1 + (-a12-opt.k23).*y2 + opt.k23*y3;
Jty3 = opt.D3 * My_Laplace(y3, h) + ...
       opt.k31*y1 + opt.k32*y2-(opt.k31+opt.k32)*y3;
   
Jtu = [Jty1; Jty2; Jty3];
end

function F = My_Laplace(u, h)
du = diff(u);
F = ([du; 0] - [0; du]) / h^2;
end