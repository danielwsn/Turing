function [perf,info]=GHiSD2(F_func,x0,V0,k,opt)
% add a v0=ones(n,1)

%% prepare
N=opt.N;
x = x0;
n = length(x0);
dt = opt.dt1;
ds = opt.dt2;
l = opt.l;
pp=0;

% gamma=0.9;  %heavy ball , gamma=0 back to zero momentum

gamma=0;



%% initialize of V

[~,kk]=size(V0);
V = V0;
V = My_orth(V);
if opt.s
    perf.x=zeros(n,opt.max/opt.m);
end
info.nF=zeros(1,opt.max);

ii=opt.max; % save the # of steps
x_old=x0;


%% iterate

for i=1:opt.max
    if (i>opt.max*3.9/4 && pp==0 && norm(F_func(x,opt))>10*opt.eps1) %|| (i>7*opt.max/8 & pp==1 & norm(F_func(x,opt))>2*opt.eps1)
        dt=dt/1000;
        ds=ds/1;
        pp=pp+1;
    end
    F = F_func(x,opt);
    nF=norm(F);
    info.nF(i)=nF;
    if nF<opt.eps1  %&& i>opt.max*1/2
        ii=i-1;
        disp(['find! i=', num2str(i)])
        break;
    elseif nF>5e2
        F=F/nF * 5e2;
    end
    x_oldtemp=x;
    x = x+dt*(F-2*V(:,2:k+1)*(V(:,2:k+1)'*F))+gamma*(x-x_old);
    x_old=x_oldtemp;
    for j=2:kk
        V(:,j)=V(:,j)+ds*(F_func(x+l*V(:,j),opt)-...
               F_func(x-l*V(:,j),opt))/(2*l);
    end
    V = My_orth(V);
    if opt.s
        perf.x(:,i)=x;
    end



    %% plot k-saddle
    if mod(i,opt.m)==0  
        if opt.s  % save or not
            perf.x(:,i/opt.m)=x;    
        end
        plot(x(1+N:2*N))
        xlim([1,N]);
        title(['step=',num2str(i),' g=',num2str(nF)])
        drawnow
    end



end






%% output
if ~opt.s
    perf.x=x;
else
    perf.x=perf.x(:,1:floor(ii/opt.m));
end
perf.V=V;
info.step=ii;
info.nF=info.nF(1:ii);

end

function nV=My_orth(V)
[~,m]=size(V);
if m==1
    nV=V;
    return 
end
nV=V;
nV(:,1)=V(:,1)/norm(V(:,1));
for j=2:m
    q=V(:,j);
    for t=1:j-1
        q=q-q'*nV(:,t)*nV(:,t);
    end
    q=q/norm(q);
    for t=1:j-1
        q=q-q'*nV(:,t)*nV(:,t);
    end
    q=q/norm(q);
    nV(:,j)=q;
end

end
