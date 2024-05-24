function v = velocity(x,y)

velocity = zeros(1,size(x,2));
for i = 1:size(x,2)
    velocity(i) = allen(x(1,i),y(1,i),4,500,1000000,1);
end

velocity(isnan(velocity)) = 0;

v = velocity;

end

