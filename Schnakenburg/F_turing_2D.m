function F = F_turing_2D(x,opt)
% x is a 2*N^2 vector
N = opt.N;
h = opt.L/N;
u = reshape(x(1:N^2),N,N);
v = reshape(x(N^2+1:2*N^2),N,N);
eta = opt.eta;
a = opt.a;
b = opt.b;
d = opt.d;
%% the force
F1 = eta*(a-u+(u.^2).*v);
F2 = eta*(b-(u.^2).*v);

F1 = F1 + 1 * My_Laplace2D(u, h, N);
F2 = F2 + d * My_Laplace2D(v, h, N);

%% output
F1= reshape(F1,N^2,1);
F2= reshape(F2,N^2,1);
F = [F1; F2];% out put 2*N^2 vector
end

function F = My_Laplace2D(u, h, N)
du = diff(u);
ddu = diff(u');
z = zeros(1,N);
F = ([du; z] - [z; du]+[ddu; z]'-[z; ddu]') / h^2;
end