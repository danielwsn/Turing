function [perf,info]=HiOSD_NG(F_func,x0,V0,k,opt)
%% prepare
x=x0;
n = length(x0);
dt = opt.dt1;    
ds = opt.dt2;   
l = opt.l;


%% initialize of V
% V = randn(n,k);
[~,kk]=size(V0);
V = V0;
V = My_orth(V);
if opt.s
    perf.x=zeros(n,opt.max/opt.m);   
end
info.nF=zeros(1,opt.max);
ii=opt.max; % save the # of steps
dtt=0;
pp=0;
%% iterate
for i=1:opt.max
    if (i>opt.max/3 && pp==0 && norm(F_func(x,opt))>10*opt.eps1) %|| (i>9*opt.max/10 & pp==1 & norm(F_func(x,opt))>2*opt.eps1)
        dt=dt/100;
        ds=ds/1;
        pp=pp+1;
    end
    F= F_func(x,opt);
    nF=norm(F);
    info.nF(i)=nF;
    if nF<opt.eps1  % stop or not
        ii=i-1;
        fprintf('find!,F=%f\n',nF);
        
        break;
    end
    if nF>1e1 
        F=F/nF;     
    end
    x = x+dt*(F-2*V(:,1:k)*(V(:,1:k)'*F));
    for j=1:kk
        V(:,j)=V(:,j)+ds*(F_func(x+l*V(:,j),opt)-...
               F_func(x-l*V(:,j),opt))/(2*l);
    end
    V = My_orth(V);
%    dtt=dtt+dt;


%     if opt.s
%         perf.x(:,i)=x;
%     end
    %% plot       
    if mod(i,opt.m)==0
        if opt.s
            perf.x(:,i/opt.m)=x;   
        end
        pcolor(reshape(x(1:opt.N^2),opt.N,opt.N))
            axis equal
            axis off
            colormap(jet)
            colorbar
            shading interp
            title([num2str(opt.dt1*i),' g=',num2str(nF)])          
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
if m==0
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
