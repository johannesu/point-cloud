classdef (Abstract) PointCloud < handle
    
    
    methods (Abstract)
        generate_global_proposals(self, num_random_proposals);
    end
    
    properties (Hidden)
        pairwise_cost;
        max_neighbors;
    end
    
    properties (SetAccess = protected)
        points;
        dimensions;
        neighborhood;
        num_points;
    end
    
    properties (Access = public)
        assignments;
        data_weight;
        tol;
        
        % For plot function
        tangent_width;
        tangent_thickness;
        tangent_color;
        point_size;
        point_color;
        input_point_size;
        input_point_color;
        
        verbose;
        
        cost_function;
        eps;
    end
    
    methods
        function self = PointCloud(points, neighborhood, data_weight, tol)
            
            self.points = points;
            self.num_points = size(self.points,2);
            
            self.assignments = [points; -sum(points.*points)];
            self.neighborhood = neighborhood;
            self.data_weight = data_weight;
            self.tol = tol;
            
            % Defaults
            self.tangent_width = 0.075;
            self.tangent_thickness = 1;
            self.tangent_color = [0 0 0];
            self.point_size = 1;
            self.point_color = [0 0 1];
            self.input_point_size = 1;
            self.input_point_color = [1 0 0];
        
            self.cost_function = 'default';
            self.eps = 1e-7;
            
            self.verbose = false;
            
            % Calculate most neighbors any point
            self.max_neighbors = 0;
            for i = 1:self.num_points
                self.max_neighbors = max(sum(neighborhood.ind1 == i), self.max_neighbors);
                self.max_neighbors = max(sum(neighborhood.ind2 == i), self.max_neighbors);
            end
            
            % Add path to include
            addpath(fileparts(mfilename('fullpath')));
            addpath([fileparts(mfilename('fullpath')),filesep,'include']);
            
            %% Compile if need be
            compile_script('energy');
            compile_script('local_optimization');
            
            % Compile if need be
            sources = {['QPBO-v1.3.src' filesep 'QPBO.cpp'], ...
                ['QPBO-v1.3.src' filesep 'QPBO_extra.cpp'], ...
                ['QPBO-v1.3.src' filesep 'QPBO_maxflow.cpp'], ...
                ['QPBO-v1.3.src' filesep 'QPBO_postprocessing.cpp']};
            
            compile_script('fusion', sources);
        end
        
        % Plot functions
        function plot_energy(self)
            [E,U,B] = self.energy;
            title(sprintf('Total cost: %2.3f -- data: %2.3f,  regularization: %2.3f', E,U,B));
        end
        
        
        % Fuse one proposal
        function fused_variables = binary_fusion(self, proposal)
            
            % If the proposals is only one plane
            if (size(proposal,2) == 1)
                proposal = repmat(proposal, [1 self.num_points]);
            end
            
            if (all(size(proposal) ~= size(self.assignments)))
                error('Proposal should be same size as assignment, or only one tangent plane');
            end
            
            improve = false;
            [S, E, LB, num_unlabelled] = fusion_mex(...
                self.cost_function, ...
                self.dimensions, ...
                self.assignments,  ...
                proposal, ...
                self.points, ...
                self.connectivity - 1, ...
                self.neighborhood.value, ...
                self.data_weight, ...
                self.tol, ...
                self.eps, ...
                self.verbose, ...
                improve);
            
            fused_variables = sum(S);
            
            if (self.verbose)
                fprintf('Fusion info -- fused %d variables, energy: %g lb: %g unlabeled: %d.\n', ...
                    fused_variables, E, LB, num_unlabelled);
                
                self.plot();
                drawnow();
            end
            
            % Update
            self.assignments(:,S==1) = proposal(:,S==1);
        end
        
        % fuse in random order
        function batch_binary_fusion(self, proposals)
            
            if (~isa(proposals,'cell'))
                error('Proposals must be a cell array');
            end
            
            for i = 1:randperm(numel(proposals))
                self.binary_fusion(proposals{i});
            end
            
        end
        
        function fuse_until_convergence(self, proposals)
            
            if (~isa(proposals,'cell'))
                error('Proposals must be a cell array');
            end
            
            fused = 1;
            while (fused > 0)
                fused = 0;
                for i = 1:randperm(numel(proposals))
                    fused = fused + self.binary_fusion(proposals{i});
                end
            end
        end
        
        function set.cost_function(self, cost_function)
            switch(cost_function)
                case {'linear','quadratic', 'default', 'length'}
                    self.cost_function = cost_function;
                otherwise
                    error('Regularization can be either linear,quadratic or default');
            end
        end
        
        % Normalize input
        function set.assignments(self, assignments)
            self.assignments = self.normalize_assignments(assignments);
        end
        
        function [energy,data_cost,pairwise_cost, data_term, pairwise_terms] = energy(self, assignments)
            
            if (nargin < 2)
                assignments = self.assignments;
            else
                assignments = self.normalize_assignments(assignments);
            end
            
            
            if (nargout < 4)
                [energy, data_cost, pairwise_cost] = ...
                    energy_mex(...
                    self.cost_function, ...
                    self.dimensions, ...
                    assignments,  ...
                    self.points, ...
                    self.connectivity - 1, ...
                    self.neighborhood.value, ...
                    self.data_weight, ...
                    self.tol, ...
                    self.eps, ...
                    self.verbose);
            else
                [energy, data_cost, pairwise_cost, data_term, pairwise_terms] = ...
                    energy_mex(...
                    self.cost_function, ...
                    self.dimensions, ...
                    assignments,  ...
                    self.points, ...
                    self.connectivity - 1, ...
                    self.neighborhood.value, ...
                    self.data_weight, ...
                    self.tol, ...
                    self.eps, ...
                    self.verbose);
            end
            
        end
        
        function num = num_pairwise(self)
            num = numel(self.neighborhood.ind1);
        end
        
        function local_optimization(self, max_iterations)
            if nargin < 2
                max_iterations = 5000;
            end
            
            
            new_assignments = local_optimization_mex( ...
                self.cost_function, ...
                self.dimensions, ...
                self.assignments,  ...
                self.points, ...
                self.connectivity - 1, ...
                self.neighborhood.value, ...
                self.data_weight, ...
                self.tol, ...
                self.eps, ...
                self.verbose, ...
                int32(max_iterations));
            
            if (self.energy(new_assignments) < self.energy)
                self.assignments = new_assignments;
            end
        end
        
        function conn = connectivity(self)
            conn = uint32([self.neighborhood.ind1'; self.neighborhood.ind2']);
        end
        
        % Some local optimization and fusion
        function optimize(self, iterations, num_global_props)
            
            if nargin < 3
                num_global_props = 500;
            end
            
            for iter = 1:iterations
                % Generate some proposals
                mean_local = self.generate_mean_proposals();
                mean_random = self.generate_global_proposals(num_global_props);
                
                %% Optimize
                self.local_optimization;
                
                self.batch_binary_fusion(mean_local);
                self.batch_binary_fusion(mean_random);
                
                
                if (self.verbose)
                    self.plot();
                    drawnow();
                end
            end
        end
        
        % Geneate new proposals by averaging over neighbors
        function proposals = generate_mean_proposals(self)
            
            % Intiallize with current assignments
            % for neighborhood with less than max_neighor neighbors
            % this will be the placeholder assigment.
            
            proposals = cell(self.max_neighbors,1);
            for i = 1:self.max_neighbors
                proposals{i} = self.assignments;
            end
            
            for i = 1:self.num_points
                
                %Mean of neghboring planes
                p = 0;
                for j = self.neighborhood.ind2(self.neighborhood.ind1 == i)';
                    p = p+1;
                    p1 = self.points(:,i);
                    p2 = self.points(:,j);
                    
                    q = (p1+p2)./2;
                    
                    n1 = self.assignments(1:end-1,i);
                    
                    n2 = self.assignments(1:end-1,j);
                    d2 = self.assignments(end,j);
                    
                    if abs(acos(n1'*n2)) > pi/2
                        n2 = -self.assignments(1:end-1,j);
                        d2 = -self.assignments(end,j);
                    end
                    dm2 = -(d2+q'*n2);
                    
                    nnew = mean([n1 n2],2);
                    nnew = nnew./norm(nnew);
                    dnew = -dm2-nnew'*q;
                    
                    proposals{p}(:,i) = [nnew; dnew];
                end
            end
        end
    end
    
    methods (Access = protected)
        
        function assignments = normalize_assignments(~, assignments)
            denom = sqrt(sum(assignments(1:end-1,:).^2,1));
            assignments = assignments./repmat(denom,[size(assignments,1) 1]);
        end
    end
end
