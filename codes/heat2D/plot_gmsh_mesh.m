function plot_gmsh_mesh(filename)
%PLOT_GMSH_MESH  Plot a 2D Gmsh mesh exported in .m format.
%
%   plot_gmsh_mesh('mesh.m')
%
%   The function automatically detects triangular (3-node) or
%   quadrilateral (4-node) elements.

if nargin < 1
  error('Please provide the gmsh .m mesh file name.');
end

% Run the .m file to obtain the 'msh' structure
run(filename);

if ~exist('msh', 'var')
  error('The file did not define a variable named "msh".');
end

% Coordinates
coords = msh.POS(:, 1:2);   % Only x,y used for 2D

% --- Detect element type ---
if isfield(msh, 'TRIANGLES') && ~isempty(msh.TRIANGLES)
  conn = msh.TRIANGLES(:, 1:3);   % only the node indices
  elemType = 'tri';
elseif isfield(msh, 'QUADS') && ~isempty(msh.QUADS)
  conn = msh.QUADS(:, 1:4);
  elemType = 'quad';
else
  error('No TRIANGLES or QUADS found in the mesh file.');
end

% --- Plot mesh ---
figure; hold on;

switch elemType
  case 'tri'
    triplot(conn, coords(:,1), coords(:,2), 'k-');
  case 'quad'
    % Draw each quad as 4 edges
    for e = 1:size(conn,1)
      nodes = conn(e, :);
      x = coords(nodes, 1);
      y = coords(nodes, 2);
      patch(x, y, 'w', 'EdgeColor', 'k', 'FaceColor', 'none');
    end
end

axis equal;
xlabel('x'); 
ylabel('y');
title(['Mesh plot: ', elemType]);

end