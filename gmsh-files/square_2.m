%% Matlab mesh
%% square, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 9;
msh.POS = [
0 0 0;
1 0 0;
1 1 0;
0 1 0;
0.4999999999986942 0 0;
1 0.4999999999986942 0;
0.500000000002059 1 0;
0 0.500000000002059 0;
0.5000000000016824 0.5000000000016824 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 1 5 2
 5 2 2
 2 6 2
 6 3 2
 3 7 2
 7 4 2
 4 8 3
 8 1 3
];
msh.TRIANGLES =[
 8 7 4 1
 5 8 1 1
 5 2 6 1
 7 8 9 1
 8 5 9 1
 5 6 9 1
 6 3 9 1
 3 7 9 1
];
