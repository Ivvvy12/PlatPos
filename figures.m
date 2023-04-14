% plot
figure(1)% animation of density m r.s.t time t in last iertation
for i=1:length(t);
    plot(x,M(i,:),'-o','LineWidth',4)
    axis ([0 25 0 0.101])
    pause(0.1)
end
title('density M r.s.t time t')


figure(2) %strategy v
plot(v_new,'-o','LineWidth',5)


