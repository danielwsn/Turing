function opt = option_set(N, gamma, tau)
% construct the parameters opt

%% flexiable parameters
opt.N = N;
opt.gamma1 = gamma;
opt.tau = tau;
opt.tk21 = 1.2024e-4 * opt.gamma1;
opt.k13 = 0.0139 * opt.tau;
opt.k31 = 0.0416 * opt.tau;
opt.k32 = 0.00139 * opt.tau;
opt.k23 = 0.0139 * opt.tau;



%% some fixed paramters
opt.D1 = 1.8;
opt.D2 = 0.012;
opt.D3 = 1.8;
opt.utotal = 1200;
opt.L = 6;

opt.k12 = 0.5;
opt.k21 = 3.6;
opt.tk12 = 1.67e-5;






end