%calculate action of a given path: a shift of pattern in Hilbert space
F_func = @F_turing_3;
% shift,time T = 1, velocity v = shift/T 


N = 5*120;
opt.N = N;
% N = 120;
h=6/N;
shiftleft=-50;
shiftright=50;
num_shift=shiftright-shiftleft+1;
gamma=0.02;


%% solution to shift
% u0=S0_h0h.x;
u0=perf.x;

taulist=[1, 0.1, 0];  % set of preset coupling parameter tau
Smat=zeros(length(taulist),num_shift);
shiftvec= shiftleft:1:shiftright;
for i=1:length(taulist)

    tau=taulist(i);
    opt = option_set(N, gamma, tau);
    opt.dt1 = 2.5e-5;
    opt.dt2 = 6e-4;  % N = 400, opt.dt2 = 5e-5;
    opt.l = 1e-4;
    opt.max = 20e7;
    opt.eps1 = 4000e-8;
    opt.eps2 = 2e-8;
    opt.s = 0; % save or not
    opt.m=1e4;
    [perf, info] = GHiSD2(F_func, u0, zeros(3*N,1), 0, opt);
    x0=perf.x;
    S=zeros(1,num_shift);
    for j=shiftvec
        ntime=j;  
        s = ntime * h ;   % shift displacement
        S(j+1-shiftleft) = ComMid(x0,ntime,opt)*abs(s);
    end
    Smat(i,:)=S;
end




% plot
ql=shiftleft*h;
qr=shiftright*h;
x0=ql:0.0001:qr;
for i=1:length(taulist)
    x = shiftvec*h;
    y = Smat(i,1:num_shift);
    y0 = interp1(x,y,x0,'spline'); 
    plot(x0,y0);
    xlim([ql,qr])
    xticks([-0.3 0 0.3]);
    xlabel('\Delta');
    ylabel('action');
    ylim([0,3.5e5])    % different for solutions
    hold on
end
legend('\tau=1','\tau=0.1','\tau=0');










% Integration, compound midpoint
function res = ComMid(x0,n,opt)
if n==0
    res=0;
elseif n>0
    res=calcFint(x0,0,opt);
    for i=1:n-1
        res=res+2*calcFint(x0,i,opt);
    end
    res=res+calcFint(x0,n,opt);
    res=res/(2*n);
else
    n=-n;
    res=calcFint(x0,0,opt);
    for i=1:n-1
        res=res+2*calcFint(x0,-i,opt);
    end
    res=res+calcFint(x0,-n,opt);
    res=res/(2*n);
end
end



% calculate the integrand at shift d
% d : # of shift grid
function f=calcFint(x0,d,opt)
N=opt.N;
% format long;
F_func = @F_turing_3;
u1=x0(1:N);
u2=x0(N+1:2*N);
u3=x0(2*N+1:3*N);
if d>0
    ualpha1=[ones(d,1)*u1(1);u1(1:N-d)]; 
    ualpha2=[ones(d,1)*u2(1);u2(1:N-d)];
    ualpha3=[ones(d,1)*u3(1);u3(1:N-d)];
else
    d=-d;
    ualpha1=[u1(d+1:N);ones(d,1)*u1(N)];  
    ualpha2=[u2(d+1:N);ones(d,1)*u2(N)];
    ualpha3=[u3(d+1:N);ones(d,1)*u3(N)];
end
ualpha=[ualpha1;ualpha2;ualpha3];


% Derivative of ualpha
dualpha1=[ualpha1(2:N);ualpha1(N)]-[ualpha1(1);ualpha1(1:N-1)];
dualpha2=[ualpha2(2:N);ualpha2(N)]-[ualpha2(1);ualpha2(1:N-1)];
dualpha3=[ualpha3(2:N);ualpha3(N)]-[ualpha3(1);ualpha3(1:N-1)];
dualpha=[dualpha1;dualpha2;dualpha3];

f=norm(F_func(ualpha,opt))*norm(dualpha)+dualpha'*F_func(ualpha,opt);
end





























