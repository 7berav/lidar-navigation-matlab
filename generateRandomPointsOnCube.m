function points = generateRandomPointsOnCube(N)
    points = zeros(0, 3);  % May not always get N points due to rejection
    while size(points, 1) < N
        x = 2 * rand() - 1 ;  
        y = 2 * rand() - 1 ;  
        z = 1;
        z_values = [z, -z];
            
            % Store valid (x, y, z) combinations
        for z_val = z_values
            if size(points, 1) < N
                points(end+1, :) = [x, y, z_val];  % Add new point
            else
                break;
            end
        end
        y = 2 * rand() - 1 ;  
        z = 2 * rand() - 1 ;  
        
      
        x1 = 1;
                
        if x >= 0
            % Calculate real x values
            x = nthroot(x1, 1);
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
        
      
        y1 = 1;
        
        
        if y1 >= 0
            % Calculate real z values
            y = nthroot(y1, 1);
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
