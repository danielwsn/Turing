function [ind,V]=cal_index2(x,F_func,opt)
% calculate index at x
N=opt.N;
g=F_func(x,opt);
ind=0;
if norm(g)>opt.eps1*100
    disp('It is not a critical point')
    V=[];
    return
end


tau=opt.dt2;
l=1e-4;
n=length(x);
maxstep=8e7;
found=0;
k=0;
V=ones(n,1);
V = V/norm(V);
while  ~found %|| k<1 
    v=randn(n,1);
    v=orthv2V(v,V);
    for j=1:maxstep
        vv=v;
        v=v+tau*(F_func(x+l*v,opt)-F_func(x-l*v,opt))/(2*l);
        v=orthv2V(v,V);
        if mod(j,5e4)==0
            plot(v(N+1:2*N));
            title(['eigstep= ',num2str(j),' ; norm(v-vv)= ', num2str(norm(v-vv))]);
            drawnow
            xlim([1,N]);
            if norm(v-vv)<opt.eps2
                break
            end
        end
    end
    k=k+1;
    V=[V,v];
    q=v'*(F_func(x+l*v,opt)-F_func(x-l*v,opt))/(2*l);
    fprintf('%.9f\n',q);

    if q<-1e-6
        found=1;
    end
end

B=zeros(k);
for i=2:k+1
    for j=2:k+1
        q=V(:,i)'*...
       (F_func(x+l*V(:,j),opt)-F_func(x-l*V(:,j),opt))/(2*l);
        B(i,j)=q;
    end
end
lambda=eig(B);
disp('eigenvs are:')
fprintf('%.9f\n',lambda');
% disp(lambda')
ind=sum(real(lambda)>0);
fprintf('index2=%g\n',ind)

end

function v=orthv2V(v,V)
[~,k]=size(V);
for i=1:k
     v=v-v'*V(:,i)*V(:,i); % orthogonal
end
v=v/norm(v);
end
