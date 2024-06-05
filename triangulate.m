function p = triangulate(ray1,t1,ray2,t2)
%TRIANGULATE Summary of this function goes here
%   Detailed explanation goes here

% Fill A and B for first ray
A = [ray1(1)^2 - 1, ray1(1) * ray1(2), ray1(1) * ray1(3);
    ray1(1) * ray1(2), ray1(2)^2 - 1, ray1(2) * ray1(3);
    ray1(1) * ray1(3), ray1(2) * ray1(3), ray1(3)^2 - 1];
B = [(ray1(1)^2 - 1)*t1(1) + ray1(1) * ray1(2) * t1(2) + ray1(1) * ray1(3) * t1(3);
    ray1(1) * ray1(2) * t1(1) + (ray1(2)^2 - 1)*t1(2) + ray1(2) * ray1(3) * t1(3);
    ray1(1) * ray1(3) * t1(1) + ray1(3) * ray1(2) * t1(2) + (ray1(3)^2 - 1)*t1(3)];

% Fill A and B for second ray
A = A + [ray2(1)^2 - 1, ray2(1) * ray2(2), ray2(1) * ray2(3);
    ray2(1) * ray2(2), ray2(2)^2 - 1, ray2(2) * ray2(3);
    ray2(1) * ray2(3), ray2(2) * ray2(3), ray2(3)^2 - 1];
B = B + [(ray2(1)^2 - 1)*t2(1) + ray2(1) * ray2(2) * t2(2) + ray2(1) * ray2(3) * t2(3);
    ray2(1) * ray2(2) * t2(1) + (ray2(2)^2 - 1)*t2(2) + ray2(2) * ray2(3) * t2(3);
    ray2(1) * ray2(3) * t2(1) + ray2(3) * ray2(2) * t2(2) + (ray2(3)^2 - 1)*t2(3)];

p = inv(A) * B;
 
end

