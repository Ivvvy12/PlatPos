function u1 = central_diff(u,du)
% do not contain the diff of bd

index = length(u);
u1=[]
for i = 2:(index-1)
    uu=(u(i+1)-u(i-1))/(2*du)
    u1=[u1 uu]
end
u1
    