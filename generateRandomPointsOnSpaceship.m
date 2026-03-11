function points = generateRandomPointsOnSpaceship(N)
% RANDOMSPACESHIPSURFACE  N points on a “spaceship” surface:
%   - bottom: hemisphere (radius=1, z<=0)
%   - top:    cone (vertex at z=1, half-angle=45°, z>=0)
%
% Input:
%   N : total number of points
% Output:
%   points : N×3 array of [x,y,z] on the composite surface

    R = 1;
    C = [0,0,2];
    R_hemi = sqrt(C(3)^2+R^2);
    cos_cut = (0 - C(3)) / R_hemi;   % = -2/sqrt(5)
    A_hemi = 2*pi*R^2* (1 + cos_cut);           % hemisphere
    s = sqrt(R^2 + R^2);         % slant height of cone
    A_cone = pi * R * s;         % lateral area of cone
    
    % allocate counts
    N_hemi = round(N * A_hemi/(A_hemi + A_cone));
    N_cone = N - N_hemi;
    
    %---- 1) sample on hemisphere z<=0 ----

    % 균등 표면 샘플링: cosθ ~ Uniform([−1, cos_cut])
    cos_theta = -1 + (cos_cut + 1) * rand(N_hemi,1);

    sin_theta = sqrt(1 - cos_theta.^2);

    phi   = 2*pi*rand(N_hemi,1);
    x_hemi = C(1) + R_hemi * (sin_theta .* cos(phi));
    y_hemi = C(2) + R_hemi * (sin_theta .* sin(phi));
    z_hemi = C(3) + R_hemi *   cos_theta;

    
    %---- 2) sample on cone surface z>=0 ----
    % For lateral surface, uniform sampling:
    %  - phi ~ U(0,2π)
    %  - height coordinate: choose t in [0,1] with density ∝ perimeter(t)
    %    perimeter at height z = 2π * r(z), r(z)= (1 - z) * tan(α)
    alpha = pi/4;  % 45°
    % r(z) = (1 - z)*tan(alpha) = (1 - z)*1
    % so perimeter ∝ (1-z), CDF of t=z is F(z)=1 - (1-z)^2 -> invert
    v = rand(N_cone,1);
    z_cone = 1 - sqrt(1 - v);
    % radius at z
    r_cone = (1 - z_cone) * 1;
    phi2 = 2*pi*rand(N_cone,1);
    x_cone = r_cone .* cos(phi2);
    y_cone = r_cone .* sin(phi2);
    
    %---- combine ----
    points = [x_hemi, y_hemi, z_hemi; 
              x_cone, y_cone, z_cone];
end