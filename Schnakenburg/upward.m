k=1;
m=k+1;   %  m>k, set m=k+1
n=4;   %n>m

delta=5e-2;



x=saddle.p;  % k-saddle to start
Vgo=maxmode(F_func,x(:),n,opt);   
v=Vgo;  

for j=1:n-k
    plus(j).p=zeros(2*N^2,1);
    minus(j).p=zeros(2*N^2,1);
    plus(j).V=zeros(2*N^2,n+1);
    minus(j).V=zeros(2*N^2,n+1);
end
for j=2:n-k
    
    v=[Vgo(:,1:k),Vgo(:,j+k)];
    
    
    
    %  + direction
    [perfm1,infom1]=HiOSD_NG(F_func,x+delta*Vgo(:,j-1+m),v,m,opt);
    %[ind,~]=cal_index(perfm1.x,F_func,opt);
    [ind,VV]=cal_index(perfm1.x,F_func,opt);
    
    if ind>-1
        plus(j).p=perfm1.x;           
        plus(j).V=VV;
    end
    
    
    
    
    % - direction
    [perfm2,infom2]=HiOSD_NG(F_func,x-delta*Vgo(:,j-1+m),v,m,opt);
    %[ind,~]=cal_index(perfm2.x,F_func,opt);
    [ind,VV]=cal_index(perfm2.x,F_func,opt);
 
    if ind>-1
        minus(j).p=perfm2.x;           
        minus(j).V=VV;
    end
    
    
    
end

j=1;
u=plus(j).p(1:N^2);



pcolor(reshape(u,N,N))
            axis equal
            axis off
            colormap(jet)
            colorbar
            shading interp
            title('u')
            %caxis([0.2,1.8])
            drawnow
