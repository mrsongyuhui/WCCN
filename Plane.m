function [i] = Plane(N)
for i = 0:sqrt(N)
    if i^(2)+i+1 == N
        break;
    end   
end