function [vxl_val, vxl_idx, vxl_dis, V_within_radius] = sphere_vxl(center_xyz, radius, V)

    % [x, y, z] = meshgrid(1:size(V,2), 1:size(V,1), 1:size(V,3));
    % distances = sqrt((x - center_xyz(2)).^2 + (y - center_xyz(1)).^2 + (z - center_xyz(3)).^2);

    [x, y, z] = meshgrid(1:size(V,1), 1:size(V,2), 1:size(V,3));
    distances = sqrt((x - center_xyz(1)).^2 + (y - center_xyz(2)).^2 + (z - center_xyz(3)).^2);
    distances = permute(distances, [2 1 3]);

    V_within_radius = distances < radius;
    D_within_radius = V_within_radius.*distances;
    % D_within_radius(center_xyz(1), center_xyz(2), center_xyz(3)) = 1;

    vxl_idx = find(V_within_radius);
    vxl_val = V(vxl_idx);
    vxl_dis = D_within_radius(vxl_idx);

end