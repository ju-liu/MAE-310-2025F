function mesh = read_gmsh_mesh(filename, id_domain, id_dirichlet, id_neumann)
%   READ_GMSH_MESH  Read 2D mesh exported from Gmsh in .m format.
%
%   mesh = read_gmsh_mesh('mesh.m');
%   mesh = read_gmsh_mesh('mesh.m', id_domain, id_dirichlet, id_neumann);
%
%   The function reads:
%       TRIANGLES : 3-node triangles (with tags)
%       QUADS     : 4-node quads     (with tags)
%       LINES     : boundary edges
%
%   Outputs (mesh struct):
%       coords          : [N x 2]
%       tri             : [Nt x 3] domain triangles
%       quad            : [Nq x 4] domain quads
%       edge_dirichlet  : edges with tag = id_dirichlet
%       edge_neumann    : edges with tag = id_neumann

% --------- set default tags ---------
if nargin < 2 || isempty(id_domain)
  id_domain = 1;
end
if nargin < 3 || isempty(id_dirichlet)
  id_dirichlet = 2;
end
if nargin < 4 || isempty(id_neumann)
  id_neumann = 3;
end

% --------- load msh struct from .m file ---------
clear msh;        % avoid stale msh
run(filename);    % must define msh

if ~exist('msh','var')
  error('The file "%s" did not define a variable named "msh".', filename);
end

% --------- node coordinates ---------
mesh.coords = msh.POS(:, 1:2);

% --------------------------------------------------
%  TRIANGLES (3-node)
% --------------------------------------------------
mesh.tri = [];
if isfield(msh, 'TRIANGLES') && ~isempty(msh.TRIANGLES)
  % TRIANGLES: [n1 n2 n3 tag]
  tri_all = msh.TRIANGLES(:, 1:4);
  mask_tri = (tri_all(:,4) == id_domain);
  mesh.tri  = tri_all(mask_tri, 1:3);
end

% --------------------------------------------------
%   QUADS (4-node)
% --------------------------------------------------
mesh.quad = [];
if isfield(msh, 'QUADS') && ~isempty(msh.QUADS)
  % QUADS: [n1 n2 n3 n4 tag]
  quad_all = msh.QUADS(:, 1:5);
  mask_quad = (quad_all(:,5) == id_domain);
  mesh.quad = quad_all(mask_quad, 1:4);
end

% --------------------------------------------------
%   Boundary edges (LINES)
% --------------------------------------------------
mesh.edge_dirichlet = [];
mesh.edge_neumann   = [];

if isfield(msh, 'LINES') && ~isempty(msh.LINES)
  % LINES: [n1 n2 tag]
  edge_all = msh.LINES(:, 1:3);

  % Dirichlet
  mask_d = (edge_all(:,3) == id_dirichlet);
  mesh.edge_dirichlet = edge_all(mask_d, 1:2);

  % Neumann
  mask_n = (edge_all(:,3) == id_neumann);
  mesh.edge_neumann   = edge_all(mask_n, 1:2);
end
end

% EOF