function V=maxmode2(F_func,x,k,opt)
% initiate vector V0
% for conservation, add a v0=ones(n,1)
N=opt.N;
n=length(x);
V=[ones(n,1),randn(n,k)]; 
V=My_orth(V);
if k==0
    return;
end
tau=opt.dt2;
l=1e-4;
maxstep=8e7;
for i=1:maxstep
    VV=V;
    for j=2:k+1
        V(:,j)=V(:,j)+tau*(F_func(x+l*V(:,j),opt)-F_func(x-l*V(:,j),opt))/(2*l);
    end
    V=My_orth(V);

    if mod(i,5e4)==0
        plot(V(1+N:2*N,2))
        xlim([1,N]);
        title(['eigstep=',num2str(i),' ; norm(VV-V)=',num2str(norm(VV-V))])
        drawnow
    end

    if norm(VV-V)<opt.eps2
        disp(['eigstop=',num2str(i)]);
        break;
    end
end

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