function points = generateRandomPointsOnHexagonPrism(N)
    points = zeros(0, 3);  % May not always get N points due to rejection
    while size(points, 1) < N
        x = 1 * rand() ;  
        y = 1 * rand() ;  
        
      
        z2 = 1;
        
        
        if z2 > 0 
            % Calculate real z values
            z = 1;
            z_values = [z, -z];
            
            % Store valid (x, y, z) combinations
            for z_val = z_values
                if size(points, 1) < N
                    points(end+1, :) = [x-y*1/2, y*sqrt(3)/2, z_val];  % Add new point
                else
                    break;
                end
                if size(points, 1) < N
                    points(end+1, :) = [x-y*1/2, -y*sqrt(3)/2, z_val];  % Add new point
                else
                    break;
                end
                if size(points, 1) < N
                    points(end+1, :) = [-x*1/2-y*1/2, x*sqrt(3)/2-y*sqrt(3)/2, z_val];  % Add new point
                else
                    break;
                end

            end
        end
%앞뒤면
        x = 1 * rand() - 0.5 ;  
        z = 2 * rand() - 1 ;  
        
      
        y2 =3/4;
        
        
        if y2 >= 0
            % Calculate real x values
            y = nthroot(y2, 2);
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

        %옆면
        x = 1 * rand() - 0.5 ;  
        z = 2 * rand() - 1 ;  
        
      
        y2 =3/4;
        
        
        if y2 >= 0
            % Calculate real x values
            y = nthroot(y2, 2);
            y_values = [y, -y];
            
            % Store valid (x, y, z) combinations
            for y_val = y_values
                if size(points, 1) < N
                    points(end+1, :) = [x*1/2-y_val*sqrt(3)/2, x*sqrt(3)/2+y_val*1/2, z];  % Add new point
                else
                    break;
                end

            end
        end
        
        x = 1 * rand() - 0.5 ;  
        z = 2 * rand() - 1 ;  
        
      
        y2 =3/4;
        
        
        if y2 >= 0
            % Calculate real x values
            y = nthroot(y2, 2);
            y_values = [y, -y];
            
            % Store valid (x, y, z) combinations
            for y_val = y_values
                if size(points, 1) < N
                    points(end+1, :) = [x*1/2+y_val*sqrt(3)/2, -x*sqrt(3)/2+y_val*1/2, z];  % Add new point
                else
                    break;
                end
            end
        end
        

    end
end
