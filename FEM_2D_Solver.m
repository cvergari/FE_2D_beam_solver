classdef FEM_2D_Solver
    
	properties (SetAccess = public)
        coordinates  = [];      % Nodes coordinates
        connectivity = [];      % Element connectivity matrix
        
        loads = [];             % Nodal loads
        loading_steps = 1;      % Number of loading steps
        
        E   = [];               % Young's modulus
        CSA = [];               % Cross-sectional area
        I   = [];               % Moment of inertia
        BC  = [];               % Border conditions
    end
    
    
    properties (SetAccess = private)
        stiffness_matrix = [];   % Global stiffness matrix
        dof = [];               % Degrees of freedom
        n_nodes = [];           % Number of nodes
        results = struct('deformed'      , [],...
                         'displacements' , [], ...
                         'rotations'     , [] , ...
                         'forces'        , []);
    end
    
    properties (SetAccess = private, GetAccess = private)
        intermediary_solution = [];
        intermediary_bc = [];
        intermediary_F = [];
    end
    
    methods
        function this = FEM_2D_Solver()
        end
        
        
        function this = demo(~)
            % Demo application. Plots the results and provides quantitative
            % results in the "results" property.
            % A simple beam with imposed displacement an point loads.
            % Comparison with ANSYS solution is provided
            
            if nargout == 0
                msg = {'Error! This function should be called as follows' ; ...
                       'obj = FEM_2D_Solver;' ; 
                       'obj = obj.demo()';
                       'disp(obj.results.displacements'};
                error(strjoin(msg , '\n'))
            end
            
            % Reset object to input data
            this = FEM_2D_Solver;
            
            this.E   = 1e7;     % Young's modulus
            this.CSA = 1;     % Cross-sectionala area
            this.I   = 1;    % Moment of inertia
            
            % Geometry of horizontal beam
            this.coordinates  = [linspace(0,10,11) ; zeros(1,11)]';
            this.connectivity = [1:this.n_nodes-1; 2:this.n_nodes]'; % Each element connected to the next
            
            % Boundary conditions
            this.BC = [0 0 0 ; NaN(this.n_nodes-1 , 3)];  % Base constrained in displacements and rotation
            this.BC(end,:) = [NaN 1 NaN];                 % Imposed displacement at the top
            
            % Loads
            this.loading_steps = 10;                      % Apply BC and loads in 10 steps
            this.loads = zeros(this.n_nodes , 3);               % Empty load matrix
            this.loads(round(this.n_nodes/2) , :) = [0 1e6 0];  % Add load mid-shaft
            
            this = this.solveFEM();
            
            % Plot output
            figure
            hold on;
            h = this.plot_structure();
            h = h(1);
            h2 = this.plot_deformed();
            h(2) = h2(1);
            
            % Plot border conditions
            BC_ = this.BC;
            BC_(isnan(BC_)) = 0;
            h(3) = quiver(this.coordinates(end,1) , this.coordinates(end,2) , BC_(end,1) , BC_(end,2),...
                'linewidth',1,'AutoScale','on', 'AutoScaleFactor', 1);
            
            
            % Plot force
            pos = find(this.loads ~= 0);
            for j = 1 : length(pos)
                [row , ~]= ind2sub(size(this.loads) , pos(j));
                force_direction = [sign(this.loads(row , 1)) , sign(this.loads(row , 2))];
                
                h(3+j) = quiver(this.coordinates(row,1) , this.coordinates(row,2) , ...
                                0.3*force_direction(1), ...
                                0.3*force_direction(2),... 
                                'linewidth',1,'AutoScale','on', 'AutoScaleFactor', 1);
            end
            
            ansys_displacements = [zeros(11,1) , [0; 0.096; 0.3393; 0.65587 ; 0.97467; 1.2240;...
                1.3487; 1.3604; 1.2873; 1.1578; 1.0000 ]];           
            ansys_solution = ansys_displacements + this.coordinates;
            
            h(end+1) = plot(ansys_solution(:,1) , ansys_solution(:,2) , 'g--');
            
            legend(h , {'Original structure' , 'Deformed structure' , 'Displacement' , 'Forces' , 'Ansys solution'} , 'Location' , 'NorthWest');
            
            axis([-0.1 , max(this.coordinates(:,1)) + 1 , -1 , 4.5])
        end
        
        
        %%
        function this = demo2(~)
            % Demo application 2. Plots the results and provides quantitative
            % results in the "results" property.
            % This is a slightly more complex structure than demo1, with 
            % imposed displacment and point load. Comparison with
            % ANSYS solution is provided
            
            if nargout == 0
                msg = {'Error! This function should be called as follows' ; ...
                       'obj = FEM_2D_Solver;' ; 
                       'obj = obj.demo()';
                       'disp(obj.results.displacements'};
                error(strjoin(msg , '\n'))
            end
            
            % Reset object to input data
            this = FEM_2D_Solver;
            
            this.E   = 1e7;  % Young's modulus
            this.CSA = 1;    % Cross-sectionala area
            this.I   = 1;    % Moment of inertia
            
            % Geometry of horizontal beam
            x = [0:10 , 9:-1:1]';
            y = [zeros(11,1) , ; (1:3)' ; zeros(3,1)+4 ;(3:-1:1)'];
            N = size(x,1);
            this.coordinates  = [x , y];
            this.connectivity = [1:this.n_nodes-1; 2:this.n_nodes]'; 
            this.connectivity = [this.connectivity ; [N , 1]];
                
            % Boundary conditions
            this.BC = [0 0 0 ; NaN(this.n_nodes-1 , 3)];  % Base constrained in displacements and rotation
            this.BC(11,:) = [NaN 1 NaN];                 % Imposed displacement at the far side
            
            % Loads
            this.loading_steps = 10;                      % Apply BC and loads in 10 steps
            this.loads = zeros(this.n_nodes , 3);               % Empty load matrix
            this.loads(16 , :) = [0 -1e6 0];  % Add load at the to
            
            this = this.solveFEM();
            
            
            % Plot output
            figure
            hold on;
            h = this.plot_structure();
            h = h(1);
            h2 = this.plot_deformed();
            h(2) = h2(1);
            
            % Plot border conditions
            BC_ = this.BC;
            BC_(isnan(BC_)) = 0;
            h(3) = quiver(this.coordinates(end,1) , this.coordinates(end,2) , BC_(end,1) , BC_(end,2),...
                'linewidth',1,'AutoScale','on', 'AutoScaleFactor', 1);
            
            
            % Plot force
            pos = find(this.loads ~= 0);
            for j = 1 : length(pos)
                [row , ~]= ind2sub(size(this.loads) , pos(j));
                force_direction = [sign(this.loads(row , 1)) , sign(this.loads(row , 2))];
                
                h(3+j) = quiver(this.coordinates(row,1) , this.coordinates(row,2) , ...
                                0.3*force_direction(1), ...
                                0.3*force_direction(2),... 
                                'linewidth',1,'AutoScale','on', 'AutoScaleFactor', 1);
            end
            
            
            ansys_displacements = [[0;0.03;0.07;0.1;0.14;0.17;0.21;0.24;0.28;0.31;0.35;0.19;0.02;-0.14;-0.26;-0.23;-0.19;-0.12;-0.079;-0.046] , ...
                                   [0;0.0095;0.0381;0.0862;0.1543;0.2426;0.3514;0.4812;0.6324;0.8052;1;0.7221;0.4339;0.1567;-0.088;-0.2223;-0.2856;-0.2387;-0.1608;-0.0739]];           
            ansys_solution = ansys_displacements + this.coordinates;
            
            h(end+1) = plot(ansys_solution(:,1) , ansys_solution(:,2) , 'g--');
            
            legend(h , {'Original structure' , 'Deformed structure' , 'Displacement' , 'Forces' , 'Ansys solution'} , 'Location' , 'NorthEast');
            
            axis([-0.1 , max(this.coordinates(:,1)) + 1 , -1 , 4.5])
            
        end % End demo2
        
        
        
        
        function this = solveFEM(this)
            % Solver.
            % Use as this = solveFEM() to store the solution.
            
            if nargout == 0
                msg = {'Error! This function should be called as follows' ; ...
                       'obj = FEM_2D_Solver;' ; 
                       'obj = obj.solveFEM()';
                       'disp(obj.results.displacements'};
                       
                error(strjoin(msg , '\n'))
            end
            
            % First step of intermediary solution: initial coordinates and
            % zero rotations.
            this.intermediary_solution = this.coordinates;
            cumulative_rotations  = zeros(this.n_nodes,1);

            % Apply loads and displacements in steps, and at each step
            % recalculate the stiffness matrix to account for beam
            % reorientations
            this.intermediary_bc = this.BC / this.loading_steps;
            this.intermediary_F  = this.loads / this.loading_steps;
            for j = 1 : this.loading_steps

                % Assemble martix
                this.stiffness_matrix = this.assemble_matrix();

                F = this.assemble_loads();
                % Add imposed displacements as modified forces; see
                % reference within the function
                [this, F] = this.add_displacement_conditions(F);

                % Solve
                delta = this.stiffness_matrix \ F;

                % Compute displacements
                [node_disp , node_rot] = this.extract_nodal_displacements(delta);

                this.intermediary_solution = this.intermediary_solution + node_disp;
                cumulative_rotations = cumulative_rotations + node_rot;
            end

            % Store results
            this.results.deformed      = this.intermediary_solution;
            this.results.displacements = this.intermediary_solution - this.coordinates;
            this.results.rotations  = cumulative_rotations;
            
        end
        
   
        function K = assemble_matrix(this)
            % Assemble matrix
            K = zeros(max(this.dof(:)));

            for k = 1 : size(this.connectivity , 1)
                % Two nodes
                p1 = this.coordinates(this.connectivity(k,1),:);
                p2 = this.coordinates(this.connectivity(k,2),:);
                % Length
                L = norm(p1 - p2);

                % Constants
                ka  = this.CSA  * this.E/L;
                k12 = 12*this.E * this.I / L^3;
                k6  = 6*this.E  * this.I / L^2;
                k4  = 4*this.E  * this.I / L;
                k2  = 2*this.E  * this.I / L;

                % Matrix in local coordinates
                kl = [ka  , 0    , 0   , -ka , 0    , 0;...
                      0   , k12  , k6  , 0   , -k12 , k6;...
                      0   , k6   , k4  , 0   , -k6  , k2;...
                      -ka , 0    , 0   , ka  , 0    , 0;...
                      0   , -k12 , -k6 , 0   , k12  , -k6;...
                      0   , k6   , k2  , 0   , -k6  , k4];

                % Matrix in global coordinates
                kg = this.local_to_global_stiffness(p1 , p2 , kl);      
                steering_vector = [this.dof(this.connectivity(k,1),:) , this.dof(this.connectivity(k,2),:)]';

                % Add to global stiffness matrix
                K = this.assemble_K(K, kg , steering_vector);
            end
            
        end

        
        %% These are the properties set functions
        function this = set.coordinates(this, coordinates)
            this.coordinates = coordinates;
            this.n_nodes = size(coordinates , 1);
            this.dof = this.count_DoF();
            if ~isempty(coordinates)
                this = checkData(this);
            end
        end
        function this = set.BC(this, BC)
            % Stores values and checks for consistency
            this.BC  = BC;
            this.dof = this.count_DoF();   % Update Degrees of freedom
            if ~isempty(BC)
                this = checkData(this);
            end
        end
        function this = set.connectivity(this, connectivity)
            % Stores values and checks for consistency
            this.connectivity  = connectivity;
            if ~isempty(connectivity)
                this = checkData(this);
            end
        end
        function this = set.loads(this, loads)
            % Stores values and checks for consistency
            this.loads = loads;
            if ~isempty(loads)
                this = checkData(this);
            end
        end
        
         function h = plot_deformed(this)
            % Plots the structure given by coordinates and connectivity.
            % Returns handle to plot
            %   h = obj.plot_deformed();
            c = this.connectivity;
            
            h = zeros(size(c, 1) , 1);
            for k = 1 : size(c, 1)
                h(k) = plot([this.results.deformed(c(k,1) , 1) , this.results.deformed(c(k,2) , 1)] , ...
                            [this.results.deformed(c(k,1) , 2) , this.results.deformed(c(k,2) , 2)] , 'r--');
                hold on;            
           end       
           axis equal;
         end
        
        function h = plot_structure(this)
            % Plots the structure given by coordinates and connectivity.
            % Returns handle to plot
            %   h = obj.plot_structure();
            c = this.connectivity;
            
            h = zeros(size(c, 1) , 1);
            for k = 1 : size(c , 1)
                h(k) = plot([this.coordinates(c(k,1),1) ; this.coordinates(c(k,2),1)] , ...
                            [this.coordinates(c(k,1),2) ; this.coordinates(c(k,2),2)] , 'b-*');
                hold on;
            end
            axis equal;
            
        end
    end
    
    %% Private methods
    methods (Access = private)
        
        function this = checkData(this)
            % Checks for data consistency and raises errors
            
            errors = cell(0);
            if size(this.coordinates , 2) ~= 2
                errors{end+1,1} = 'Coordinates should be an Nx2 matrix';
                this.coordinates = [];
            end
            
            if ~isempty(this.BC) && size(this.BC , 2) ~= 3
                errors{end+1,1} = 'BC should be an Nx3 matrix; number of columns is wrong';
                this.BC = [];
            end
            
            if ~isempty(this.n_nodes) && this.n_nodes ~= 0 && ...
                    ~isempty(this.BC) && size(this.BC,1) ~= this.n_nodes
                errors{end+1,1} = 'BC should be an Nx3 matrix; number of rows is wrong';
                this.BC = [];
            end
                
            if ~isempty(this.loads) && size(this.loads , 2) ~= 3
                errors{end+1,1} = 'loads should be an Nx3 matrix; number of columns is wrong';
                this.loads = [];
            end
            if ~isempty(this.n_nodes) && this.n_nodes ~= 0 && ...
                    ~isempty(this.loads) && size(this.loads , 1) ~= this.n_nodes
                errors{end+1,1} = 'loads should be an Nx3 matrix; number of rows is wrong';
                this.loads = [];
            end
            
            if ~isempty(this.n_nodes) && this.n_nodes ~= 0 && ...
                    ~isempty(this.connectivity) && size(this.connectivity , 2) ~= 2
                errors{end+1,1} = 'connectivity should have two columns';
                this.connectivity = [];
            end
            
            
            if ~isempty(errors)
                message = strjoin([{'Errors were found:'} ; errors] , '\n');
                error(message)
            end
            
        end
        
        
        
        function F = assemble_loads(this)
            % Assemble the loads in the right order and positions
            F = zeros(numel(this.loads),1);
            F(1:3:end) = this.intermediary_F(:,1);   % X forces
            F(2:3:end) = this.intermediary_F(:,2);   % Y forces
            F(3:3:end) = this.intermediary_F(:,3);   % Moments
            
            % Remove forces corresponding to constraints
            ind = reshape(this.dof' , [] , 1);  % DOF in the same order as forces
            F(ind == 0) = [];
        end
        
        
        function [this, F] = add_displacement_conditions(this, F)
            % This function embeds the imposed displacements (different
            % from zero) into the force vector.
            % Procedure from:
            % A note on imposing displacement boundary conditions in finite element analysis
            % Baisheng Wu, Zhonghai Xu, Zhengguang Li
            % Communications in Numerical Methods in Engineering
            %
            % We replace the forces corresponding to the imposed displacement with an
            % equivalent force to obtain the displacement
            % Furthermore, we modify the matrix K to only have that displacement in the
            % corrisponding row & column

            where = find(~isnan(this.intermediary_bc) & this.intermediary_bc ~= 0);  % extract imposed displacement
            n = this.dof(where);   % Degree of freedom corresponding to displacement bs


            for k = 1 : length(where)
                pos = n(k);

                % equivalement force = kn * displacement
                F(pos) = this.intermediary_bc(where(k));

                % Correct other forces to account for the modification
                % f_hat_i = f_i - ki*ui
                for i = 1 : length(F)
                    if i == pos;  continue;  end
                    F(i) = F(i) - this.stiffness_matrix(i,pos) * this.intermediary_bc(where(k));
                end

                % K has a row and column of zeros, and one in the diagonal
                this.stiffness_matrix(pos,:) = 0;
                this.stiffness_matrix(:,pos) = 0;
                this.stiffness_matrix(pos,pos) = 1;
            end
        end
        
        
        function K = assemble_K(~, K, kg , g)    
            % Assembles the global stiffness matrix from one local matrix

            eldof = 2*3; % Number of degrees of freedom per element 	
            for i=1:eldof
               if g(i) ~= 0
                  for j=1: eldof
                      if g(j) ~= 0
                         K(g(i),g(j))= K(g(i),g(j)) + kg(i,j);
                      end
                   end
               end
            end 

        end
        
        
        function kg = local_to_global_stiffness(~, p1 , p2 , kl)
            % Rotates the local stiffness matrix according to the
            % orientation of the beam in the global coordinate system
            
            % Orientation of the beam
            theta = atan2(p2(2) - p1(2) , p2(1) - p1(1)); 

            % Rotation matrix
            c = cos(theta);
            s = sin(theta);
            C = [c , -s , 0 , 0  , 0  , 0 ; ...
                 s , c , 0  , 0  , 0  , 0 ; ...
                 0 , 0 , 1  , 0  , 0  , 0 ; ...
                 0 , 0 , 0  , c  , -s , 0 ;
                 0 , 0 , 0  , s  , c  , 0 ;
                 0 , 0 , 0  , 0  , 0  , 1];
            
            % Rotate local stiffness matris
            kg = C * kl * C';

        end

   

        function dof = count_DoF(this)
            % Counts and labels the degrees of freedom

            if isempty(this.n_nodes) || isempty(this.BC)
                dof = [];
                return
            end
            
            dof = ones(this.n_nodes, 3);  % Initialise the matrix nf to 1
            dof(this.BC == 0) = 0;        % Prescribed nodal freedom
            
            cnt = sum(dof(:) == 1);   % Total number of degrees of freedom
            
            % Label the DOF in the right order
            nftmp = dof';
            nftmp(nftmp == 1) = 1:cnt;  
            dof = nftmp';

        end
        
        function [node_disp , node_rot] = extract_nodal_displacements(this, delta)
            % Extract nodal displacements
            
            % Dispatch results in three columns, in the right order
            output = this.dof';             % Fixed DOF remain zero
            output(output ~=0) = delta;     % Replace other values
            output = output';
            node_disp = [output(:,1) , output(:,2)];    % X and Z displacements
            node_rot = output(:,3);                     % Rotations

        end        
        
        
    end
 
end
