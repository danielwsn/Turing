% solve the homogeneous state

%% set parameter
N = 120;
gamma = 0.02; % 0.011 and 0.02
tau = 0.01; % from 1 to 100
% opt = option_set(N, gamma, tau);

%% 
alpha = opt.k13/(opt.k31+opt.k32);
beta = opt.k23/(opt.k31+opt.k32);

p = [opt.tk21+(1+beta)*opt.tk12/(1+alpha),...
     -opt.tk12*opt.utotal/(1+alpha),...
     (1+beta)*(opt.k12+opt.k13-opt.k31*alpha)/(1+alpha)+opt.k21+opt.k31*beta,...
     opt.utotal*(opt.k31*alpha-opt.k12-opt.k13)/(1+alpha)];
 
%% plot
% x=300:1:500;
% y=polyval(p,x);
% figure()
% plot(x,y)

%% solve
r = roots(p);
disp(r)

%% 
u2 = r(1);
u1 = (opt.utotal-(1+beta)*u2)/(1+alpha);
u3 = opt.utotal-u1-u2;

uH = [u1*ones(N,1); u2*ones(N,1);u3*ones(N,1)];
g = F_turing_3(uH, opt);
disp(['norm(F) = ',num2str(norm(g))]);

%% cal index


[ind,V]=cal_index2(uH,@F_turing_3,opt);

%% prepare initial data
% perf.x = u;
% V2 = maxmode2(@F_turing_3, perf.x, 0, opt);
% perf.V = V2;

%% calculate the matrix

% syms ka la
% a11 = -opt.tk12 * u2^2 -opt.k12-opt.k13-opt.D1*ka-la;
% a22 = 2*opt.tk12*u1*u2-3*opt.tk21*u2^2-opt.k21-opt.k23-opt.D2*ka-la;
% a33 = -opt.k31-opt.k32-opt.D3*ka-la;
% a12 = -2*opt.tk12*u1*u2+3*opt.tk21*u2^2+opt.k21;
% a13 = opt.k31;
% a21 = opt.tk12*u2^2+opt.k12;
% a23 = opt.k32;
% a31 = opt.k13;
% a32 = opt.k23;
% A = [a11,a12,a13;a21,a22,a23;a31,a32,a33;];
% z = -det(A);
% % vpa(z)
% vpa(subs(z,{la},{0}))
% 
% % 0.03888 1.0201048001457585432135033443046 2.6629304610632683389383843383342 -0.00000000000000007100546996886284555244505037559
% 0.03888 5.915499165835976866426904052787 19.882953178333193280417194174827  0.00000000000000058549114004139799970648796549032



ux=uH;
plot(linspace(0,6,N),ux(1+N:2*N));
xticks([0 2 4 6])
xticklabels({'0','2','4','6'})
yticks([500 1000 1500])
xlim([0,6]);
ylim([0,1600]);
