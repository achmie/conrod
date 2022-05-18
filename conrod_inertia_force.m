% Assumptions:
% 1) The center of connection rod head is in point (0,0,0).
% 2) Axis are assigned as follows: X (s1), Y (s2), Z (s3).
% 3) Plane (x,y,0) is a symmetry plane of connection rod.

% Create model and read geometry from STL.
model = createpde('structural', 'static-solid');
importGeometry(model, 'Suzuki.stl');
pdegplot(model, 'FaceLabels', 'on', 'FaceAlpha', 0.5);
% Create mesh of tetrahedrons with maximum edge length equal to "Hmax".
mesh = generateMesh(model, 'GeometricOrder', 'linear', 'Hmax', 8);
pdeplot3D(model);
% Export mesh as: points (p), edges (e) and tetrahedrons (t).
[p,e,t] = meshToPet(mesh);

% Connection rod length [mm].
L = 100;
% Density rho = 0.0077224 [g/mm^3].
rho = 0.0077224;
% Starting volume [mm^3].
V = 0;
% Starting mass.
mc = 0;
% Starting values of moments m_{c,1} and m_{c,2};
mc1 = 0;
mc2 = 0;

% Numerical integration - summation by tetrahedrons.
for i=1:length(t(1,:))
    % Take tetrahedron vertices numbers.
    ti = t((1:4), i);
    % Convert vertices numbers to (x,y,z) coordinates:
    %   Ax, Bx, Cx, Dx
    %   Ay, By, Cy, Dy
    %   Az, Bz, Cz, Dz
    pi = [p(:,ti(1)), p(:,ti(2)), p(:,ti(3)), p(:,ti(4))];
    % Convert vertices A,B,C,D to 3 vectors with origin at A.
    wi = pi(:,2:4) - pi(:,1);
    % Volume of tetrahedron ABCD as 1/6 of vectors volume.
    Vi = abs(det(wi')) / 6; 
    V = V + Vi;
    % Compute tetrahedron mass.
    mci = rho*Vi; 
    mc = mc + mci;
    % Tetrahedron ABCD centroid.
    d = [
        sum(pi(1,:))/4;
        sum(pi(2,:))/4;
        sum(pi(3,:))/4
    ];
    % Integration by s1*dm and s2*dm
    mc1 = mc1 + d(1)*mci; 
    mc2 = mc2 + d(2)*mci;
end

% Print volume in [mm^3]
V
% Print mass in [g]
mc
% Print moment m_{c,1} in [mm*g]
mc1
% Print moment m_{c,2} in [mm*g]
mc2
