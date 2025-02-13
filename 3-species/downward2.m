% from high2low layer by layer for Turing model

F_func = @F_turing_3;

delta = 1;

%% initial point
u0=S1_th.x;


[~,V]=cal_index2(u0,F_func,opt);
dimV=size(V,2)-1;    % dimV = index+1
V0 = V(:,1:dimV);

x_downward=struct([]);


count=0;
k = 0;  % target index
opt.s=0;

%% iterate

for i = 2:dimV

    % plus direction
    u_plus = u0 + delta * V0(:,i);
    for j=[2:i-1,i+1:dimV]
        Vdown = V0(:,[1,j]);
        [perf, ~] = GHiSD2(F_func, u_plus, Vdown, k, opt);
        [ind,V]=cal_index2(perf.x,F_func,opt);
        
        if ind==1
            count=count+1;
            x_downward(count).x=perf.x;
            x_downward(count).V=V;
        end
    end
    

    % minus direction
    u_minus = u0 - delta * V0(:,i);
    for j=[2:i-1,i+1:dimV]
        Vdown = V0(:,[1,j]);
        [perf, ~] = GHiSD2(F_func, u_minus, Vdown, k, opt);
        [ind,V]=cal_index2(perf.x,F_func,opt);
        
        if ind==1
            count=count+1;
            x_downward(count).x=perf.x;
            x_downward(count).V=V;
        end
    end

end 

ux=x_downward(1).x;
plot(ux(1+N:2*N))
xlim([1,N]);