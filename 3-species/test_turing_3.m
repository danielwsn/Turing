% 3-species model
% clear
F_func = @F_turing_3;

%% parameters
N = 120; % grid



gamma = 0.02; % 0.011 and 0.02
tau = 1; 
opt = option_set(N, gamma, tau);
opt.dt1 =2.5e-5;   
opt.dt2 = 1e-4;     
opt.l = 1e-4;
opt.max = 8e8;
opt.eps1 = 50e-8;
opt.eps2 = 2e-8;
opt.s = 0; % save or not
opt.m=1e4;

k = 0;

%% inital state
% u1 = 440 * ones(N, 1);
% u2 = 466 * ones(N, 1);
% u3 = 400 * ones(N, 1);
% u0 = [u1; u2; u3];


% y = linspace(0, 6*pi, N)';
% u1 = 300 - 200 * cos(y);
% u2 = 400 - 300 * cos(y);
% u3 = opt.utotal - u1 - u2;
% u0 = [u1; u2; u3];


% %  mesh refinement q times
% N = 120;
% q=5;
% u00=S0_1h.x;
% u001=u00(1:N);
% u002=u00(N+1:2*N);
% u003=u00(2*N+1:3*N);
% u001n=zeros(q*N,1);u002n=zeros(q*N,1);u003n=zeros(q*N,1);
% for i=1:N
%     u001n(q*i-q+1:q*i)=u001(i);
%     u002n(q*i-q+1:q*i)=u002(i);
%     u003n(q*i-q+1:q*i)=u003(i);
% end
% u0=[u001n;u002n;u003n];
% N=q*N;
% opt.N=N;


u0=perf.x;


% shift
% u00=S0_3h.x;
% xshift=6;
% % u0=u00([121-xshift:120,1:120-xshift,241-xshift:240,121:240-xshift,361-xshift:360,241:360-xshift]');
% % u0=[u00(1+xshift:120);u00(120)*ones(xshift,1); u00(121+xshift:240);u00(240)*ones(xshift,1); u00(241+xshift:360);u00(360)*ones(xshift,1)];


%% initiate vector V0
V0 = maxmode2(F_func, u0, k, opt);

%% iterate
[perf, info] = GHiSD2(F_func, u0, V0, k, opt);

%% plot
figure()
semilogy(info.nF,'LineWidth',2)


ux=perf.x;
plot(linspace(0,6,N),ux(1+N:2*N));
xticks([0 2 4 6])
xticklabels({'0','2','4','6'})
ylim([0,1600]);
yticks([500 1000 1500])

% xtickangle(0)
% ytickangle(0)



% cal_index2(ux,F_func,opt)




