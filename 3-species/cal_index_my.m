function [ind,V]=cal_index_my(x,F_func,opt)
% calculate index at x

N=opt.N;
g=F_func(x,opt);
ind=-1;
tau=opt.dt2;
l=1e-4;
n=length(x);
maxstep=5e7;
found=0;
k=0;
V=zeros(n,0);

if norm(g)>opt.eps1*100
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
        plot(v(1+N:2*N,1));
                xlim([1,N]);
        title(['step=',num2str(j),'; ',num2str(norm(vv-v))]);
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
ind=sum(real(lambda)>0);       %index of saddle
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
