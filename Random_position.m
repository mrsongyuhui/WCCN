function [positions] = Random_position(base_station_location, point_number, r_c, r_h)
y = real(base_station_location{1});
x = imag(base_station_location{1});
% generate MS randomly
i=1;
while i <= point_number
    point_x = (2*rand(1)-1)*r_c;
    point_y = (2*rand(1)-1)*r_c;
    condition1 = abs(point_x) + abs(point_y)/sqrt(3) <= r_c;
    condition2 = point_x^(2)+point_y^(2) >= r_h^(2);
    condition3 = abs(point_y) <= sqrt(3)/2*r_c;
    if  condition1 && condition2 && condition3
        positions(i) = y + point_y + 1i*(x + point_x);
        i=i+1;
    end
end
end