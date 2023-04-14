function v_new = Int(t,x,dx,M,v)
v_new = [];
for i = 1 : length(t);
    vv = 0;
    for j = 2 : length(x)-1; % Simpson's Rule
        vv = vv + dx/3 * M(i,j-1) + 4/3 * dx * M(i,j) + dx/3 * M(i,j+1) ;
    end
    vv = v(i) * vv;
    v_new = [v_new  vv];
end
