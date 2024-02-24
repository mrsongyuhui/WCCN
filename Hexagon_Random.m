function [coordinate] = Hexagon_Random(r_c, r_h, reuse_factor, point_number, layer_number )
coordinate = cell(1,layer_number);
coordinate{1} = num2cell([0;0+1i*0]);

% on the line
for i = 2:layer_number
    theta = 2*pi/(6*(i-1));
    gap = sqrt(3)/2*2*(i-1)*r_c;
    coordinate_temp = cell(2,6*(i-1));
    for j = 1:(i-1):6*(i-1)+1
        x = sin(theta*(j-1))*gap;
        y = cos(theta*(j-1))*gap;
        coordinate_temp{2,j} = (  y+1i*x );
    end
    coordinate{i} = coordinate_temp;
end

% not on the line
for i = 3:layer_number
    coordinate_temp = coordinate{i};
    for j = 1:(i-1):6*(i-1)
        for k = 1:(i-1)
            coef0 = (i-1-k)/(i-1);
            coef1 = k/(i-1);
            coordinate_temp{2,j+k} = coef0*coordinate_temp{2,j} + coef1*coordinate_temp{2,j+i-1};
        end
    end
    coordinate{i} = coordinate_temp;
end

% delete
for i = 2:layer_number
    coordinate_temp = coordinate{i};
    coordinate_temp_final = coordinate_temp(:,1:6*(i-1));
    coordinate{i} = coordinate_temp_final;
end

% reuse factor
for i = 2:layer_number
    coordinate_temp = coordinate{i};
    for j = 1:6*(i-1)
        temp = coordinate_temp{2,j};
        y = real(temp);
        x = imag(temp);
        u = round(x*2/3/r_c);
        v = round( (y-u*r_c/2*sqrt(3))*cos(pi/6) *2/3/r_c);
        p = Plane(reuse_factor);
        number = mod( (p+1)*u+v,reuse_factor);
        coordinate_temp{1,j} = number;
    end
    coordinate{i} = coordinate_temp;
end

% generate MS
for i = 1:layer_number
    coordinate_temp = coordinate{i};
    if i == 1
        j = 1;
        base_station_location = coordinate_temp(2,j);
        coordinate_temp{3,j} = Random_position(base_station_location, point_number, r_c, r_h);
    end
    for j = 1:6*(i-1)
        base_station_location = coordinate_temp(2,j);
        coordinate_temp{3,j} = Random_position(base_station_location, point_number, r_c, r_h);
    end
    coordinate{i} = coordinate_temp;
end