classdef PointCloud2D < PointCloud
    
    properties
        
    end
    
    methods (Static)
        function neighborhood = generate_neighborhood(points, maxdist)
            
            tri = delaunay(points(1,:),points(2,:));
            ind =[tri(:,1:2); tri(:,2:3); tri(:,[3 1])];
            ind = [ind; ind(:,[2 1])];
            ind = unique(ind,'rows');
            
            x = points(:,ind(:,1));
            y = points(:,ind(:,2));
            eukldist = sqrt(sum((x-y).^2));
            
            
            neighborhood.ind1 = ind(eukldist < maxdist,1);
            neighborhood.ind2 = ind(eukldist < maxdist,2);
            
            nrneighb = accumarray(neighborhood.ind1,1,[size(points,2) 1]);
            neighbdist = sqrt(sum((points(:,neighborhood.ind1)-points(:,neighborhood.ind2)).^2));
            meanneighbdist = accumarray(neighborhood.ind1,neighbdist,[size(points,2) 1])./nrneighb;
            neighborhood.value = 1./nrneighb(neighborhood.ind1).*meanneighbdist(neighborhood.ind1);
        end
    end
    
    methods
        function self = PointCloud2D(points, neighborhood, data_weight, tol)
            
            if (size(points,1) ~= 2)
                error('Points should be formated as 2 x num_points');
            end
            
            self = self@PointCloud(points, neighborhood, data_weight, tol);
            self.dimensions = uint32(2);
        end
        
        function q = points_from_proposal(self, proposal)
            
            if nargin < 2
                proposal = self.assignments;
            end
            
            n = proposal(1:2,:);
            d = proposal(3,:);
            q = self.points - ([1;1]*(sum(n.*self.points)+d)).*n;
        end
        
        function plot_neighborhood_structure(self)
            
            % Highlight truncated terms
            [~,~,~,~, pairwise_terms] = self.energy();
            cost = pairwise_terms./self.neighborhood.value;
            thresholeded_terms = (cost >= self.tol);
            
            % Neighboring points
            p1 = self.points(:,self.neighborhood.ind1);
            n1 = self.assignments(1:2,self.neighborhood.ind1);
            d1 = self.assignments(3,self.neighborhood.ind1);
            q1 = p1 - repmat(d1+sum(n1.*p1),[2 1]).*n1;
            
            p2 = self.points(:,self.neighborhood.ind2);
            n2 = self.assignments(1:2,self.neighborhood.ind2);
            d2 = self.assignments(3,self.neighborhood.ind2);
            q2 = p2 - repmat(d2+sum(n2.*p2),[2 1]).*n2;
            
            for i = 1:size(q1,2)
                if (thresholeded_terms(i))
                    plot([q1(1,i) q2(1,i)],[q1(2,i) q2(2,i)],'-r');
                else
                    plot([q1(1,i) q2(1,i)],[q1(2,i) q2(2,i)],'-g');
                end
            end
        end
        
        function plot_input_points(self)
            plot( self.points(1,:), self.points(2,:),'.','Color', self.input_point_color, 'MarkerSize', self.input_point_size);
        end
        
        function plot_points(self)
            % Each point
            p = self.points;
            n = self.assignments(1:2,:);
            d = self.assignments(3,:);
            q = p - repmat(d+sum(n.*p),[2 1]).*n;
            plot(q(1,:),q(2,:),'.', 'Color', self.point_color, 'MarkerSize', self.point_size);
        end
        
        function plot_tangents(self)
            for i = 1:self.num_points;
                p = self.points(:,i);
                n = self.assignments(1:2,i);
                d = self.assignments(3,i);
                
                qq = p - (n'*p+d)*n;
                
                v = [n(2); -n(1)];
                plot([qq(1)-self.tangent_width*v(1) qq(1)+self.tangent_width*v(1)],  ...
                    [qq(2)-self.tangent_width*v(2) qq(2)+self.tangent_width*v(2)], ...
                    'color', self.tangent_color,'linewidth', self.tangent_thickness);
            end
        end
        
        % Generate proposals
        function proposals = generate_global_proposals(self, number_of_proposals)
            
            proposals = {};
            
            n = self.assignments(1:end-1,:);
            d = self.assignments(end,:);
            
            qq = self.points - (ones(2,1)*(sum(n.*self.points)+d)).*n;
            
            
            while size(proposals,2) < number_of_proposals/2;
                randpt = ceil(self.num_points*rand(2,1));
                p1 = self.points(:,randpt(1));
                p2 = self.points(:,randpt(2));
                
                v = null([p1 p2; ones(1,2)]');
                
                if size(v,2) == 1
                    v = v./norm(v(1:end-1));
                    proposals{end+1} = v;
                end
            end
            
            
            while size(proposals,2) < number_of_proposals;
                randpt = ceil(self.num_points*rand(3,1));
                p1 = qq(:,randpt(1));
                p2 = qq(:,randpt(2));
                
                v = null([p1 p2; ones(1,2)]');
                
                if size(v,2) == 1
                    v = v./norm(v(1:end-1));
                    proposals{end+1} = v;
                end
            end
        end
        
	% Can be used to produce high-quality figures for publication using pgfplots+latex.
        function save_for_pgfplots(self, name)
            % Pgfplots use:
            % \addplot table[col sep = comma] {filename.tsv}
     
           fname_points = sprintf('%s_points.tsv', name);
           csvwrite(fname_points, self.points()');
           
 
           fname_points = sprintf('%s_tangetlines.tsv', name);
           fid = fopen(fname_points, 'w');

           for i = 1:self.num_points
               p = self.points(:,i);
               n = self.assignments(1:2,i);
               d = self.assignments(3,i);
               qq = p - (n'*p+d)*n;
               v = [n(2); -n(1)];
               
               sz = self.tangent_width;
               fwrite(fid,sprintf('%g,%g\n', qq(1) - sz*v(1), qq(2) - sz*v(2)));
               fwrite(fid,sprintf('%g,%g\n', qq(1) + sz*v(1), qq(2) + sz*v(2)));
               fwrite(fid,sprintf('\n'));
               
           end
           fclose(fid);

        end

        function plot(self, text)
            clf; hold on;

            if (nargin < 2)
                text = '';
            end
            
            self.plot_input_points();           
            self.plot_tangents();
            legend('Initial points', 'Estimated tanget lines','Location','northoutside');


            self.plot_energy();
            xlabel(text);
            axis equal;
        end
    end
end
