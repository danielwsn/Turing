function V=maxmode(F_func,x0,k,opt)    
N=opt.N;
n=length(x0);
x=x0;
% k=k+1;     
V=randn(n,k); V=orth(V);     

if k==0
    return;
end
tau=opt.dt2;
l=1e-4;
maxstep=1e6;

for i=1:maxstep
    VV=V(:,1);
    for j=1:k  
        V(:,j)=V(:,j)+tau*(F_func(x+l*V(:,j),opt)-F_func(x-l*V(:,j),opt))/(2*l);
    end
    
    V=orth(V);      
    
    if norm(VV-V(:,1))<opt.eps2
        fprintf('i_stop = %d \n',i);
        break;
    end
    
    if mod(i,1e3)==0
        
        pcolor(reshape(V(1:N^2,1),N,N))
        colormap(jet)
        colorbar
        shading interp
        title([num2str(i),'; ',num2str(norm(VV-V(:,1)))])
        drawnow
        
    end
    
end
end



