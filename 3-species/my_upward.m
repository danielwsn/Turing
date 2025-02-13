% index-0 to index-1

delta = 10;

%% initial point
u0 = S0_1.x;
m=1;  % m>=1
V0=maxmode2(F_func,u0,m,opt);  %  dim=m+1 



x_upward=struct([]);
count=0;

opt.s=0;



%% iterate
k = 1;  % target index, pre-fixed
for i = 2


    % plus direction
    u_plus = u0 + delta * V0(:,i);

    Vup = V0(:,[1,i]);
    [perf, ~] = GHiSD2(F_func, u_plus, Vup, k, opt);
    [ind,V]=cal_index2(perf.x,F_func,opt);
    
    if ind==1
        count=count+1;
        x_upward(count).x=perf.x;
        x_upward(count).V=V;
    end



    % minus direction
    u_minus = u0 - delta * V0(:,i);

    Vup = V0(:,[1,i]);
    [perf, ~] = GHiSD2(F_func, u_minus, Vup, k, opt);
    [ind,V]=cal_index2(perf.x,F_func,opt);
    
    if ind==1
        count=count+1;
        x_upward(count).x=perf.x;
        x_upward(count).V=V;
    end

end 

ux=x_upward(2).x;
plot(linspace(0,6,N),ux(1+N:2*N));
xticks([0 2 4 6])
xticklabels({'0','2','4','6'})
yticks([500 1000 1500])
% xtickangle(0)
% ytickangle(0)




