function points = generateRandomPointsOnSurface(N)
    points = zeros(0, 3);  % May not always get N points due to rejection
    while size(points, 1) < N
        x = 2 * rand() - 1 ;  
        y = 2 * rand() - 1 ;  
        
      
        z4 = (1 - (x)^4 - (y)^4);
        
        
        if z4 >= 0
            % Calculate real z values
            z = nthroot(z4, 4);
            z_values = [z, -z];
            
            % Store valid (x, y, z) combinations
            for z_val = z_values
                if size(points, 1) < N
                    points(end+1, :) = [x, y, z_val];  % Add new point
                else
                    break;
                end
            end
        end
        y = 2 * rand() - 1 ;  
        z = 2 * rand() - 1 ;  
        
      
        x4 = (1 - (y)^4 - (z)^4);
        
        
        if x4 >= 0
            % Calculate real x values
            x = nthroot(x4, 4);
            x_values = [x, -x];
            
            % Store valid (x, y, z) combinations
            for x_val = x_values
                if size(points, 1) < N
                    points(end+1, :) = [x_val, y, z];  % Add new point
                else
                    break;
                end
            end
        end
        z = 2 * rand() - 1 ;  
        x = 2 * rand() - 1 ;  
        
      
        y4 = (1 - (z)^4 - (x)^4);
        
        
        if y4 >= 0
            % Calculate real z values
            y = nthroot(y4, 4);
            y_values = [y, -y];
            
            % Store valid (x, y, z) combinations
            for y_val = y_values
                if size(points, 1) < N
                    points(end+1, :) = [x, y_val, z];  % Add new point
                else
                    break;
                end
            end
        end
    end
end
