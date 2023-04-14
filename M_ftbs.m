function M = M_ftbs(init,x,t,N1,Nt,v_star,dx,dt)

ml = init;
M = [ml];
mm = 1 : length(x) ;
ml_der=[]
for i = 1 : length(t) - 1; %time
    for j = 2 : length(x); %space
        %Nt = (N^(1/alfa)+t(i))^alfa;
        %N1 = alfa*(N^(1/alfa)+t(i))^(alfa-1); %derivative of Nt
        
        d11 = ml(j);
        d12 = N1(i)/Nt(i);
        ml1 = (v_star(i,j)-v_star(i,j-1))/dx;
        d1= -(ml1+d12) * d11;
         
        d21 = v_star(i);
        d22 = x(j)*N1(i)/Nt(i);
        d23 = (ml(j)-ml(j-1))/dx;
        d2 = -(d21+d22)*d23;
         
        mm(j) = (d1+d2)*dt+ml(j);
    end
    %bc
    mm(1) = 0;
    ml=mm;
    M =[M;ml];
end
