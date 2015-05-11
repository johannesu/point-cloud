addpath([fileparts(mfilename('fullpath')),filesep,'..']);

%% Generate data
% Chooice n random points on the sphere
rng(2)
n = 1000;

th = 2*pi*rand(1,n);
ph = asin(-1+2*rand(1,n));
r = 1/4;

[x,y,z] = sph2cart(th,ph,r);
points = [x;y;z];

% noise
points = points + 0.02.*([1;1;1]*randn(1,size(points,2)));
neighborhood = PointCloud3D.generate_neighborhood(points, 1.5e-1);

%% Setup problem
data_weight = 10;
tol = 10;
P = PointCloud3D(points, neighborhood, data_weight, tol);
P.tangent_width = 0.025;
start_solution = P.assignments;

%% Solve
figure(1); clf; hold on;
P.fancyplot('Start solution');
view([120 30]); 

P.optimize(5);

figure(2); clf; hold on;
P.fancyplot('After optimization');
view([120 30]);