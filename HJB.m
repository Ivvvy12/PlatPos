function nu = HJB(N,alfa,beta,eta,rho,t,k,x,dx,init,dt,Nt,N1,v_or,z,T,z0,p0)

%LT = (N^(1/alfa) + t).^alfa; %free boundary
%dx = LT/ns;
%dt = T/nt;
%x = 0:dx:LT;
%t = 0:dt:T;

% bc & ic
h = (1+exp(x)).^(-1);
h0 = 0
Terminal = exp(-beta * T) * h;
bdl = exp(beta .* -t) .* h0;
bdh = exp(beta .* -t) .* ((N.^(1/alfa)+t).^alfa).^(1/2);

wl = Terminal;
W = [Terminal];
ww = 1 : length(x) ;
nu1 = 1 : length(x) ;
Xp=[x];
Tp=[t];
nu=[]
l = x
z1 = central_diff(z,dt);
zl = (z(2)-z(1))/dt;
zr = (z(end)-z(end-1))/dt
z1 = [zl z1 zr];

for i = length(t) - 1:-1:1; %time
    %L = (N^(1/alfa) + t(i))^alfa;
    %dx = L/ns;
    %x = 0:dx:L;
    Xp=[x;Xp];
    Tp=[t;Tp];
    wl1 = exp(-beta * t(i)) * 0;
    wl2 = exp(-beta * t(i)) * t(i)^alfa;
    for j = 2 : length(x) - 1; %space
        %Nt = (N^(1/alfa)+t(i))^alfa;
        %N1 = alfa*(N^(1/alfa)+t(i))^(alfa-1); %derivative
        %z = sin(t(i)^2);
        %z1= cos(t(i)^2)*2*t(i);  %derivative
        
        P_bar = exp(-beta * t(i))* (p0 - eta * (z(i)-z0));
        P_bar1 = (-beta) * P_bar - exp(-beta * t(i)) * eta * z1(i)
        
        
        G1 = l(i) * exp(-beta * t(i));
        G2 = x(j) * P_bar1;
        G3 = (x(j)*N1(i))/Nt(i) * exp(-beta * t(i)).*(p0 - eta * (z(i)-z0));
        G4 = rho*beta * x(j)^2 *exp(-beta * t(i))/Nt(i);
        G = G1+G2+G3+G4;
        
        ww(j) = -(exp(-beta * t(i)) * N1(i)/(4*rho) *  (wl(j) - wl(j-1))/dx * (wl(j+1) - wl(j))/dx + G);
        ww(j) = wl(j) - ww(j) * dt;
        nu1(j) = ww(j) + x(j) * P_bar - rho* x(j)^2 *exp(-beta * t(i))/Nt(i)
        nu1(1)= wl1 + x(j) * P_bar - rho* x(j)^2 *exp(-beta * t(i))/Nt(i)
        nu1(end)= wl2 + x(j) * P_bar - rho* x(j)^2 *exp(-beta * t(i))/Nt(i)
    end
    %bc
    
    
    %ww(1) = wl1
    %ww(length(ww)) = wl2;
%     ww = [wl1 ww wl2]
%     wl=ww;
    %W =[W;wl];?
    nu=[nu1;nu] 
end

nu



