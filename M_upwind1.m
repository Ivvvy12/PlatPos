function M = M_upwind1(init,X,t,N1,Nt,v_star,dx,dt)

ml = init;
M = [ml];
mm = 1 : length(X(1,:)) ;
ml_der=[]
for i = 1 : length(t) - 1; %time
    for j = 2 : length(X(i,:))-1; %space
        %Nt = (N^(1/alfa)+t(i))^alfa;
        %N1 = alfa*(N^(1/alfa)+t(i))^(alfa-1); %derivative of Nt
        
        d11 = ml(j);
        d12 = N1(i)/Nt(i);
        ml1 = (v_star(i,j)-v_star(i,j-1))/dx;
        d1= -(ml1+d12) * d11;
         
        d21 = v_star(i);
        d22 = X(i,j)*N1(i)/Nt(i);
        if d21+d22 > 0 
            d23 = (ml(j)-ml(j-1))/dx;
        else 
            d23 = (ml(j+1)-ml(j))/dx;
        end
            d2 = -(d21+d22)*d23;
         
        mm(j) = (d1+d2)*dt+ml(j);
    end
    %bc
    mm(1) = 0;
    mm(end)=mm(end-1)
    ml=mm;
    M =[M;ml];
end
