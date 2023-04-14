function z = Z_eq1(N,alfa,t,k,X,dx,init,dt,Nt,N1,v_or)


z = [];
for i = 1 : length(t) ; %time
        p11 = 1 / N;%(N?)
        p12 = 0;
        for j = 2 : length(X(i,:))-1; % Simpson's Rule
            p12 = p12 + dx/3 * X(i,j-1)*init(j-1) + 4/3 * dx * X(i,j)*init(j) + dx/3 * X(i,j+1)*init(j+1) ;
        end
        p1 = -p11 * p12;
        
        p21 = k;
        p22 = 0;
        if i == 1;
            p22 = 0;
        else;
            for j = 1 : i-1;
                p22 = p22 + 1/2*(v_or(j)/Nt(j) + v_or(j+1)/Nt(j+1)) *dt;
            end
        end
        p2 = -p21 * p22;
        zz = (1+p1+p2)*Nt(i);
        z = [z zz];
end
