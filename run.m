
xv = -250:1:250;
yv = 0:1:500;
for i = 1:size(xv,2)
    for j = 1:size(yv,2)
        velocity(i,j) = allen(xv(1,i),yv(1,j),4,500,1000000,1);
    end
end
velocity(isnan(velocity)) = 0;

surf(yv,xv,velocity,LineStyle="none");
%hold on
%plot(x,BC(x)),