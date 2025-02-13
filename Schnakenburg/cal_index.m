function [ind,V]=cal_index(x,F_func,opt)
% calculate index at x
g=F_func(x,opt);
ind=-1;
N=opt.N;
tau=opt.dt2/5;
l=1e-4;
n=length(x);
maxstep=5e6;
found=0;
k=0;
V=zeros(n,0);
if norm(g)>opt.eps1*2
    disp('It is not a critical point')
    return
end






while ~found
    
    v=randn(n,1);
    v=orthv2V(v,V);
%     figure()


    for j=1:maxstep
        vv=v;
        v=v+tau*(F_func(x+l*v,opt)-F_func(x-l*v,opt))/(2*l);
        v=orthv2V(v,V);
        
        if norm(v-vv)<opt.eps2
           break; 
        end
        
        if mod(j,1e4)==0
            pcolor(reshape(v(1:N^2),N,N))
            colormap(jet)
            colorbar
            shading interp
            title(['vector',num2str(j),' ; ',num2str(norm(v-vv))])
            drawnow
        end
    end
    
    
    k=k+1;
    fprintf('k = %d ;',k);fprintf('norm(v-vv) = ');
    disp(norm(v-vv))
    
    
    
    V=[V,v];
    q=v'*(F_func(x+l*v,opt)-F_func(x-l*v,opt))/(2*l);
    fprintf('q is %f\n',q);
    fprintf('\n');
    fprintf('\n');
    if q<-1e-5
        found=1;
    end
    
    
end


B=zeros(k);

for i=1:k
    for j=1:k
        q=V(:,i)'*...
       (F_func(x+l*V(:,j),opt)-F_func(x-l*V(:,j),opt))/(2*l);
        B(i,j)=q;
    end
end


lambda=eig(B);
disp(lambda')
ind=sum(real(lambda)>0);       % index of saddles
fprintf('index2=%g\n',ind);
fprintf('！！！！！！！！！！！！！！！！！！！！\n');
end



function v=orthv2V(v,V)
[~,k]=size(V);
for i=1:k
     v=v-v'*V(:,i)*V(:,i); % orthogonal
end
v=v/norm(v);
end
