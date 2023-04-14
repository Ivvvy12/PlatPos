clear all ,clf,clc

global ns nt T L LL m0 alfa beta rho N p0 z0 eta k 

%parameters
ns = 10; % # of space nodes 
nt = 50; % # of time nodes
T = 2;   %terminal
L = 5;
num=5;
LL = num*L;
m0 = 1/(ns); %uniform distribution
alfa = 0.9;
beta = 0.1;
rho = 0.25;
N = 50;
p0 = 1;
z0 = 0.1;
eta = 0.01;
k = 1;


%mesh
dx = L/ns;
dt = T/nt;
x = 0:dx:LL;
x1 = 0:dx:(num/2)*L;
x2 = ((num/2)*L + dx):dx:(num/2+1)*L;
x3 = ((num/2+1)*L + dx):dx:LL;
t = 0:dt:T;
LT = (N^(1/alfa) + t).^alfa; %free boundary

% bc & ic
init1 = ones(1,length(x1)).* 0;
init2 = ones(1,length(x2)).* m0; %initial condition
init3 = ones(1,length(x3)).* 0;
init = [init1 init2 init3];
bd = ones(1,length(t)).*0; %boundary condition

Nt = (N^(1/alfa)+t).^alfa;
N1 = alfa*(N^(1/alfa)+t).^(alfa-1); %derivative of Nt
v_or = t.* 0.01;% first guess for v_or (?m0)


%step 1 z(t)
z = Z_eq(N,alfa,t,k,x,dx,init,dt,Nt,N1,v_or)

%step 2 HJB equation
nu = HJB(N,alfa,beta,eta,rho,t,k,x,dx,init,dt,Nt,N1,v_or,z,T,z0,p0)

% step3 v_star
% v = ones(1,length(x)).* 0.2;
%P_bar = p0-eta*(z-z0)
%v = N1./(2*rho).*(-P_bar);
v_star = v_star(N1,rho,p0,eta,z,z0,beta,t,nu,dx)

%step4 transport equation m(t,x)
%M = M_ftbs(init,x,t,N1,Nt,v_star,dx,dt)
%M = M_fd_central(init,x,t,N1,Nt,v_star,dx,dt)
M = M_upwind(init,x,t,N1,Nt,v_star,dx,dt)


%step5
v_new = Int(t,x,dx,M,v_star)

%step6 iteration
error_inf =[ max(abs(v_new - v_or))]
error2 = [norm(v_new - v_or,'fro')]

M_storage = [M]
v_storage = []
iteration = 2;
for it = 1 : iteration;

    v_storage = [v_storage v_or];
    v_or = v_new;
    
    z = Z_eq(N,alfa,t,k,x,dx,init,dt,Nt,N1,v_or);
    
    %HJB
    nu = HJB(N,alfa,beta,eta,rho,t,k,x,dx,init,dt,Nt,N1,v_or,z,T,z0,p0)


    %v_star = v_star(N1,rho,p0,eta,z,z0,beta,t,nu,dx)
    v_eq1=[]
for i=1:length(t)
dnul=(nu(2,i)-nu(1,i))/dx;
dnur=(nu(end,i)-nu(end-1,i))/dx;
dnu = [dnul central_diff(nu,dx) dnur];
v1 = N1./(2*rho).*(exp(-beta.*t(i)).*dnu-p0-eta.*(z-z0))

v_eq1 =[v_eq1;v1];
end

v_star=v_eq1

    %M = M_ftbs(init,x,t,N1,Nt,v_star,dx,dt)
    %M = M_fd_central(init,x,t,N1,Nt,v_star,dx,dt)
    M = M_upwind(init,x,t,N1,Nt,v_star,dx,dt)
    M_storage = [M_storage M]
    
    v_new = Int(t,x,dx,M,v_star);
    error_inf = [error_inf max(abs(v_new - v_or))]; % error in infinity norm
    error2 = [error2 norm(v_new - v_or,'fro')]; % error in L2 norm
    
end







