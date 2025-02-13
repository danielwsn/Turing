function F = F_turing_3_reaction(u, opt)
% three species
%% seperate the field
N = opt.N;
h = opt.L/N;
u1 = u(1:N);
u2 = u(N+1:2*N);
u3 = u(2*N+1:end);

%% the force
F1 = -opt.tk12 * u1 .* u2.^2 + opt.tk21 * u2.^3 - opt.k12 * u1 + ...
    opt.k21 * u2;
F2 = -F1 - opt.k23 * u2 + opt.k32 * u3;
F1 = F1 - opt.k13 * u1 + opt.k31 * u3;
F3 = -F1 - F2;
%% output
F = [F1; F2; F3];
end
