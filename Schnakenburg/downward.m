% downward search 
%k to m, m<k

k=1;

m=0;


delta=2e-2;



x=perf.x; % start point
V=perf.V; % dim=k


% r=5;
% x=saddleC1p(r).p;
% V=saddleC1p(r).V;


S=zeros(2*k,4); % record
for i=1:2*k
    sink(i).p=zeros(2*N^2,1);
    sink(i).V=zeros(2*N^2,m+1);
end

for i=1:k
    if i>m 
        v=V(:,1:m);
    else
        v=V(:,1:m+1);
        v(:,[i])=[];
    end
    
    
    
    %  + direction
    [perfm1,infom1]=HiOSD_NG(F_func,x+delta*V(:,i),v,m,opt);
    %[ind,~]=cal_index(perfm1.x,F_func,opt);
    [ind,VV]=cal_index(perfm1.x,F_func,opt);
    figure()
    semilogy(infom1.nF,'LineWidth',2)
    title(['norm of force, k=',num2str(k), ' m=',num2str(m),' i=',num2str(i),'']);
    xlabel('step')
    ylabel('norm')
    if ind>-1
        S(2*i-1,1)=k;S(2*i-1,2)=m;S(2*i-1,3)=ind;S(2*i-1,4)=i;
        sink(2*i-1).s=S(2*i-1,:);  
        sink(2*i-1).p=perfm1.x;             
        sink(2*i-1).V=VV;
    end
    
    
    
    
    %  - direction
    [perfm2,infom2]=HiOSD_NG(F_func,x-delta*V(:,i),v,m,opt);
    %[ind,~]=cal_index(perfm2.x,F_func,opt);
    [ind,VV]=cal_index(perfm2.x,F_func,opt);
    figure()
    semilogy(infom2.nF,'LineWidth',2)
    title(['norm of force, k=',num2str(k), ' m=',num2str(m),' i=',num2str(i),'-']);
    xlabel('step')
    ylabel('norm')
    if ind>-1
        S(2*i,1)=k;S(2*i,2)=m;S(2*i,3)=ind;S(2*i,4)=i;
        sink(2*i).s=S(2*i-1,:);  
        sink(2*i).p=perfm2.x;            
        sink(2*i).V=VV;
    end

end

j=1;
u=sink(j).p(1:N^2);



pcolor(reshape(u,N,N))
            axis equal
            axis off
            colormap(jet)
%             colorbar
            shading interp
%             title('u')
       caxis([0.7,1.8])
            drawnow
