% 2D schnakenburg
F_func=@F_turing_2D;


%% parameters
opt.N=32; 
% N = opt.N;
opt.L=1;
% opt.h=1/opt.N;
% opt.y=(1:opt.N)/opt.N;


opt.dt1=5e-4;      

opt.dt2=5e-6;     


opt.l=1e-4;
opt.eps1=1e-4;    % force threshold
opt.eps2=3e-8;    % eigen threshold
opt.max=4.5e6;
opt.m=2e3;    %check point
opt.s=0; % save the dynamic or not 

%schnakenberg model parameter
opt.d=37;
opt.eta=200;
opt.a=1/3;   
opt.b=2/3;   


k=1;  % index to search

%% initialize
x=(1:opt.N)/opt.N;
N = opt.N;
[xx,yy]=meshgrid(x,x);


% u0 = (opt.a + opt.b) * ones(N,N);
% v0 = opt.b * ones(N,N);
% x0 = [u0(:);v0(:)];

  
% u0 = 1.2*ones(N,N);
% v0 = 0.65*ones(N,N);
% u0(opt.N/4:3*opt.N/4,opt.N/4:3*opt.N/4)=0.7;
% x0 = [u0(:);v0(:)];        %%

   
% u0 = 1.5*ones(N,N);
% v0 = 0.5*ones(N,N);
% u0(opt.N/4:3*opt.N/4,opt.N/4:3*opt.N/4)=0.3;
% x0 = [u0(:);v0(:)];          %%

   
% u0 = 0.8*ones(N,N);
% v0 = 0.6*ones(N,N);
% u0(opt.N/4:3*opt.N/4,opt.N/4:3*opt.N/4)=1.1;
% x0 = [u0(:);v0(:)];          %%

 
% u0 = 0.8*ones(N,N);
% v0 = 0.7*ones(N,N);
% x0 = [u0(:);v0(:)];       


% u0 = rand(N,N);
% v0 = rand(N,N);
% x0 = [u0(:);v0(:)];


%% build unstable eigen space
V0=maxmode(F_func,x0(:),k,opt);     


%% HiSD
[perf, info]=HiOSD_NG(F_func,x0(:),V0,k,opt);


%% plot
if ~isempty(info.nF)
    figure()
    semilogy(info.nF,'LineWidth',2)
    title('norm of force')
    xlabel('step')
    ylabel('norm')
    drawnow
end

u=perf.x(1:N^2);
v=perf.x(N^2+1:2*N^2);

pcolor(reshape(u,N,N))
            axis equal
            axis off
            colormap(jet)
            %colorbar
            shading interp
%             title('u')
            caxis([0.7,1.8])
            drawnow




%% calculate index
[ind,V]=cal_index(perf.x(:,end),F_func,opt);     