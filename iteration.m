function [v_eq, M, error_inf, error_l2]=itertation(itertion,v_new,v_or,M,N,alfa,t,k,x,dx,init,dt,Nt,N1,rho,)

error_inf =[ max(abs(v_new - v_or))]
error2 = [norm(v_new - v_or,'fro')]

M_storage = [M]
v_storage = []
iteration = 20;
for it = 1 : iteration;
    v_storage = [v_storage v_or];
    v_or = v_new;
    
    z = Z_eq(N,alfa,t,k,x,dx,init,dt,Nt,N1,v_or);
    
    %HJB
    
    v_star = v_star(N1,rho,p0,eta,z,z0,beta,t,nu);
    
    M = M_ftbs(init,x,t,N1,Nt,v_star);
    
    v_new = Int(t,x,dx,M);

