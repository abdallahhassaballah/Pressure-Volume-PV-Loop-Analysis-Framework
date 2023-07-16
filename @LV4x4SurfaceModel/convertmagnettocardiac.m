function points =convertmagnettocardiac(transform, data, n)
a1 = transform(1,:)
a2 = transform(2,:)
a3 = transform(3,:)
t = transform(4,:)
%convert each point
for i = 1:n
    temp = data(i,:) - t_;
    points(i,:) = (temp(1)*a1 +temp(2)*a2 +temp(3)*a3)
end
end