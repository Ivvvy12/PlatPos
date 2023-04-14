function v_eq = v_star(N1,rho,p0,eta,z,z0,beta,t,nu,dx)

v_eq1=[]
for i=1:length(t)
dnul=(nu(2,i)-nu(1,i))/dx;
dnur=(nu(end,i)-nu(end-1,i))/dx;
dnu = [dnul central_diff(nu,dx) dnur];
v1 = N1./(2*rho).*(exp(-beta.*t(i)).*dnu-p0-eta.*(z-z0))

v_eq1 =[v_eq1;v1];
end

v_eq=v_eq1

