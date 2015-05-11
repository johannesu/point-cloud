classdef PointCloud3D < PointCloud
    
    properties
       point_cut_size;
       displacement_color;
       displacement_thickness;
    end
    
    methods (Static)
        function neighborhood = generate_neighborhood(points, maxdist)
            
            tet = delaunay(points(1,:),points(2,:),points(3,:));
            ind = [tet(:,[1 2]); tet(:,[1 3]); tet(:,[1 4]); tet(:,[2 3]); tet(:,[2 4]); tet(:,[2 4])];
            ind = [ind; ind(:,[2 1])];
            ind = unique(ind,'rows');
            
            x = points(:,ind(:,1));
            y = points(:,ind(:,2));
            eukldist = sqrt(sum((x-y).^2));
            
            neighborhood.ind1 = ind(eukldist < maxdist,1);
            neighborhood.ind2 = ind(eukldist < maxdist,2);
            nrneighb = accumarray(neighborhood.ind1,1);
            neighborhood.value = 1./nrneighb(neighborhood.ind1);
        end
        
        function T = plane_projection(pi)
            %change coordinates
            T1 = eye(4);
            T1(1:3,4) = pi(4)*pi(1:3);
            [U,~,V] = svd(pi(1:3));
            if det(U*V) > 0
                T2 = [U*V zeros(3,1); zeros(1,3) 1];
            else
                T2 = [-U*V zeros(3,1); zeros(1,3) 1];
            end
            
            T = T2'*T1;
        end
        
        % planepoints is logical index
        function points = project_points_to_tagent_plane(points, pi)
            
            pi = pi./norm(pi(1:3));
            points = points - pi(1:3,:)*(pi'*[points; ones(1,size(points,2))]);
            
            T = PointCloud3D.plane_projection(pi);
            
            points = T*[points; ones(1,size(points,2))];
            points= points./(ones(4,1)*points(end,:));
            
            %points in plane
            points = points(2:3,:);
        end
        
        function proposals = project_proposal_to_tangent_plane(proposals, pi)
            pi = pi./norm(pi(1:3));
            T = PointCloud3D.plane_projection(pi);
            
            proposals = inv(T)'*proposals;
            proposals = proposals./([1;1;1;1]*sqrt((proposals(1,:).^2+proposals(2,:).^2+proposals(3,:).^2)));
            
            %line intersection between assigned patches and cut
            proposals = proposals(2:4,:);
            proposals = proposals./([1;1;1]*sqrt((proposals(1,:).^2+proposals(2,:).^2)));
        end
    end
    
    methods
        function self = PointCloud3D(points, neighborhood, data_weight, tol)
            
            if (size(points,1) ~= 3)
                error('Points should be formated as 3 x num_points');
            end
            
            %Remove points that don't have any neighbor
            keepind = unique(neighborhood.ind1);
            invkeepind = 1:size(points,2);
            invkeepind(keepind) = 1:length(keepind);
            points = points(:,keepind);
            neighborhood.ind1 = invkeepind(neighborhood.ind1);
            neighborhood.ind2 = invkeepind(neighborhood.ind2);
            keepind = (neighborhood.ind1 ~= 0) & (neighborhood.ind2 ~= 0);
            neighborhood.ind1 = neighborhood.ind1(keepind)';
            neighborhood.ind2 = neighborhood.ind2(keepind)';
            
            self = self@PointCloud(points, neighborhood, data_weight, tol);
            self.dimensions = uint32(3);
            
            % Settings
            self.point_cut_size = 1;
            self.displacement_color = [1 0.4 0.4];
            self.displacement_thickness = 1;
        end
        
        function plot_input_points(self)
            hold on;
            plot3( self.points(1,:), self.points(2,:),self.points(3,:), ...
                '.' , 'MarkerSize', self.input_point_size, 'Color', self.input_point_color);
        end
        
        function plot_points(self)
            
            q = self.points_from_proposal();
            plot3(q(1,:),q(2,:),q(3,:), '.', 'MarkerSize', self.point_size, 'Color', self.point_color);
        end
        
        function plot_neighborhood_structure(self)
            % Neighboring points
            p1 = self.points(:,self.neighborhood.ind1);
            n1 = self.assignments(1:3,self.neighborhood.ind1);
            d1 = self.assignments(4,self.neighborhood.ind1);
            q1 = p1 - repmat(d1+sum(n1.*p1),[3 1]).*n1;
            
            p2 = self.points(:,self.neighborhood.ind2);
            n2 = self.assignments(1:3,self.neighborhood.ind2);
            d2 = self.assignments(4,self.neighborhood.ind2);
            q2 = p2 - repmat(d2+sum(n2.*p2),[3 1]).*n2;
            
            for i = 1:size(q1,2)
                plot3([q1(1,i) q2(1,i)],[q1(2,i) q2(2,i)],[q1(3,i) q2(3,i)],'-b');
            end
        end
        
        function plot_displacement(self)
            q = points_from_proposal(self);
            plot3([self.points(1,:)' q(1,:)']',[self.points(2,:)' q(2,:)']',[self.points(3,:)' q(3,:)']','-', ...
                'color', self.displacement_color, 'linewidth', self.displacement_thickness);
        end
        
        function plot_tangents(self)
            % Display tangents
            q = points_from_proposal(self);
            n = self.assignments(1:3,:);
            v1 = [n(2,:); -n(1,:); zeros(1,size(n,2))];
            v1 = v1./repmat(sqrt(sum(v1.^2)),[3 1]);
            v2 = cross(v1,n,1);
            v1 = self.tangent_width*v1;
            v2 = self.tangent_width*v2;
            
            X = [q(1,:)+v1(1,:)+v2(1,:); q(1,:)+v1(1,:)-v2(1,:); q(1,:)-v1(1,:)-v2(1,:); q(1,:)-v1(1,:)+v2(1,:)];
            Y = [q(2,:)+v1(2,:)+v2(2,:); q(2,:)+v1(2,:)-v2(2,:); q(2,:)-v1(2,:)-v2(2,:); q(2,:)-v1(2,:)+v2(2,:)];
            Z = [q(3,:)+v1(3,:)+v2(3,:); q(3,:)+v1(3,:)-v2(3,:); q(3,:)-v1(3,:)-v2(3,:); q(3,:)-v1(3,:)+v2(3,:)];
            %fill3(X,Y,Z,'b');
            fill3(X,Y,Z,Z,'edgecolor',[0.25 0.25 0.25]);
                    
            axis equal;
            drawnow;
        end
        
        function q = points_from_proposal(self, proposal)
            
            if nargin < 2
                proposal = self.assignments;
            end
            
            n = proposal(1:3,:);
            d = proposal(4,:);
            q = self.points - ([1;1;1]*(sum(n.*self.points)+d)).*n;
        end
        

        
        %% Proposals
        function proposals = generate_global_proposals(self, number_of_proposals)
            
            qq = self.points_from_proposal();
            
            proposals = {};
            while size(proposals,2) < number_of_proposals/2;
                randpt = ceil(self.num_points*rand(3,1));
                p1 = self.points(:,randpt(1));
                p2 = self.points(:,randpt(2));
                p3 = self.points(:,randpt(3));
                
                randpt = ceil(self.num_points*rand(3,1));
                p1 = self.points(:,randpt(1));
                p2 = self.points(:,randpt(2));
                p3 = self.points(:,randpt(3));
                
                v = null([p1 p2 p3; ones(1,3)]');
                if size(v,2) == 1
                    v = v./norm(v(1:3));
                    proposals{end+1} = v;
                end
            end
            while size(proposals,2) < number_of_proposals/2;
                randpt = ceil(self.num_points*rand(3,1));
                p1 = qq(:,randpt(1));
                p2 = qq(:,randpt(2));
                p3 = qq(:,randpt(3));
                
                randpt = ceil(self.num_points*rand(3,1));
                p1 = self.points(:,randpt(1));
                p2 = self.points(:,randpt(2));
                p3 = self.points(:,randpt(3));
                
                v = null([p1 p2 p3; ones(1,3)]');
                if size(v,2) == 1
                    v = v./norm(v(1:3));
                    proposals{end+1} = v;
                end
            end
        end
        
        % Returns logical index
        function indices = index_of_points_close_to_plane(self, pi, distance)
            
            % Data
            pi = pi./norm(pi(1:3));
            q = self.points_from_proposal();
            
            indices = abs(pi'*[q; ones(1,size(q,2))]) < distance;
            
            if (~any(indices))
                warning('No points at distance to plane');
            end
        end
        
        function plot_points_close_to_plane(self, pi, distance, color)
            
            q = self.points_close_to_plane(pi,distance);
            plot3(q(1,:),q(2,:),q(3,:),'.','color', color,'MarkerSize', self.point_size);
        end
        
        function q = points_close_to_plane(self, pi, distance)
            planepoints = self.index_of_points_close_to_plane(pi, distance);
            q = self.points_from_proposal();
            q = q(:,planepoints);
        end
        
        function plot_cut(self,pi, distance, color)
                        
            q = self.points_from_proposal();
            planepoints = self.index_of_points_close_to_plane(pi, distance);
            
            qpi = PointCloud3D.project_points_to_tagent_plane(q(:,planepoints), pi);
            pointspi = PointCloud3D.project_points_to_tagent_plane(self.points(:,planepoints),pi);
            proposalspi = PointCloud3D.project_proposal_to_tangent_plane(self.assignments(:,planepoints), pi);
            
            
            for i = 1:size(qpi,2);
                plot(pointspi(1,i),pointspi(2,i),'.', 'color', color ,'MarkerSize', self.point_cut_size);
            end
            
            sz = self.tangent_width;
            for i = 1:size(qpi,2);
                p = qpi(:,i);
                n = proposalspi(1:2,i);
                d = proposalspi(3,i);
                qq = p - (n'*p+d)*n;
                v = [n(2); -n(1)];
                plot([qq(1)-sz*v(1) qq(1)+sz*v(1)],[qq(2)-sz*v(2) qq(2)+sz*v(2)], ...
                    'color', self.tangent_color, 'linewidth', self.tangent_thickness);
            end
            
        end
        
        function pi = get_random_hyperplane(self)
            q = self.points_from_proposal;
            ids = randperm(size(q,2),3);
            pi = null([q(:,ids); ones(1,3)]');
        end

        % Faster default plot
        function plot(self, text)
            clf; hold on;

            if (nargin < 2)
                text = '';
            end
            
            self.plot_input_points();
            self.plot_points();
            legend('Initial points', 'Estimated points','Location','northoutside');
            
            self.plot_energy();
            xlabel(text);
            axis equal;
        end

        % Slower plot showing the tangent plane, for large problem this
        % figure will be interoable slow in matlab
        function fancyplot(self, text)
            clf; hold on;

            if (nargin < 2)
                text = '';
            end
            
            self.plot_displacement(); 
            self.plot_tangents();
                          
            legend('Point displacement','Location','northoutside');
            self.plot_energy();
            xlabel(text);
            axis equal;
        end
    end
end