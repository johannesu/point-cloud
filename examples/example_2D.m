close all;
addpath([fileparts(mfilename('fullpath')),filesep,'..']);

%% Genearte data
rng(0);

phi = linspace(0, 2*pi, 39);
phi = phi(1:end-1);
points = [cos(phi); sin(phi)];

n = diag([1 1])*points;
n = n./([1;1]*sqrt(sum(n.^2)));

% Adding noise
points = points+n*0.025.*([1;1]*randn(1,size(n,2)));

% Neighborhood structure
neighborhood.ind1 = [[1:size(points,2)]'; [2:size(points,2)]'; 1];
neighborhood.ind2 = [[2:size(points,2)]'; 1;[1:size(points,2)]'];
neighborhood.value = ones(size(neighborhood.ind1));

% Settings
tol = 100;
data_weight = 0.5;

P = PointCloud2D(points, neighborhood, data_weight, tol);
P.input_point_size = 5;
intial_assignments = P.assignments;
iterations = 10;


figure(1); clf;
P.assignments = intial_assignments;
P.cost_function = 'linear';
P.optimize(iterations);
P.plot('Linear regularization: sharp corners');

figure(2); clf,
P.assignments = intial_assignments;
P.cost_function = 'length';
P.optimize(iterations);
P.plot('Length regularization: shrinking bias');

figure(3); clf;
P.assignments = intial_assignments;
P.cost_function = 'quadratic';
P.optimize(iterations);
P.plot('Quadratic regularization: growing bias');

figure(4); clf;
P.assignments = intial_assignments;
P.cost_function = 'default';
P.optimize(iterations);
P.plot('Weighted quadratic regularization');