classdef FEM_3D_Solver
    
	properties
        coordinates  = [];      % Nodes coordinates
        connectivity = [];      % Element connectivity matrix
        
        loads = [];             % Nodal loads [Fx, Fy, Fz, Mx, My, Mz]
        loading_steps = 1;      % Number of loading steps
        
        CSA = [];               % Cross-sectional area
        G   = [];               % Shear modulus
        J   = [];               % Torsional constant (Polar Moment of Inertia)
        E   = [];               % Young's modulus
        Iy  = [];               % Moment of inertia y
        Iz  = [];               % Moment of inertia z        
        BC  = [];               % Boundary conditions [ux1; uy1; uz1; rotx1; roty1; rotz1; ux2; uy2; uz2; rotx2; roty2; rotz2;]
    end
    
    
    properties (SetAccess = private)
        stiffness_matrix = [];  % Global stiffness matrix
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
        function this = FEM_3D_Solver()

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
            this = FEM_3D_Solver;
            
            radius   = 5.642;     % Radius for a CSA of 100
            this.E   = 1e7;     % Young's modulus
            this.CSA = pi*radius^2;     % Cross-sectionala area
            this.Iy  = pi*radius^4/4;    % Moment of inertia
            this.Iz  = pi*radius^4/4;    % Moment of inertia
            this.G   = this.E / (2*(1+0.3));  % Shear modulus assuming Poisson ratio of 0.3
            this.J   = pi*radius^4/2;

            % Geometry of horizontal beam
            this.coordinates  = [linspace(0,500,11) ; zeros(1,11); zeros(1,11)]';
            this.connectivity = [1:this.n_nodes-1; 2:this.n_nodes]'; % Each element connected to the next
            
            % Boundary conditions
            this.BC = [0 0 0 0 0 0 ; NaN(this.n_nodes-1 , 6)];    % Base constrained in displacements and rotation    
            
            % Loads
            this.loading_steps = 20;                      % Apply BC and loads in 10 steps
            this.loads = zeros(this.n_nodes , 6);         % Empty load matrix
            this.loads(ceil(this.n_nodes/2) , :) = [0, 5e4, 0, 0, 0, 0];  % Add load mid-shaft
            this.loads(end , :) = [0, 0, 0, 0, 5e6, 0];   % Add moment at beam end
            
            
            this = this.solveFEM();
            
            % Plot output
            figure
            hold on;
            h = this.plot_structure();
            h = h(1);
            h2 = this.plot_deformed();
            h(2) = h2(1);

            % Plot force
            % pos = find(this.loads ~= 0);
            % for j = 1 : length(pos)
            %     [row , ~]= ind2sub(size(this.loads) , pos(j));
            %     force_direction = [sign(this.loads(row , 1)) , sign(this.loads(row , 2)), sign(this.loads(row , 3))];
            % 
            %     h(3+j) = quiver3(this.coordinates(row,1) , this.coordinates(row,2) , this.coordinates(row,3) , ...
            %                     100*force_direction(1), ...
            %                     100*force_direction(2),... 
            %                     100*force_direction(3),... 
            %                     'linewidth',1,'AutoScale','on', 'AutoScaleFactor', 1);
            % 
            %     % Moment
            %     moment_direction = [sign(this.loads(row , 4)) , sign(this.loads(row , 5)), sign(this.loads(row , 6))];
            %     h(3+j) = quiver3(this.coordinates(row,1) , this.coordinates(row,2) , this.coordinates(row,3) , ...
            %                     100*moment_direction(1), ...
            %                     100*moment_direction(2),... 
            %                     100*moment_direction(3),... 
            %                     'linewidth',1,'AutoScale','on', 'AutoScaleFactor', 1);                
            %     quiver3(this.coordinates(row,1) , this.coordinates(row,2) , this.coordinates(row,3) , ...
            %             80*moment_direction(1), ...
            %             80*moment_direction(2),... 
            %             80*moment_direction(3),... 
            %             'linewidth',1,'AutoScale','on', 'AutoScaleFactor', 1,...
            %             'Color', get(h(3+j), 'Color'));                   
            % end    
            
            % Displacements obtained from ansys simulation
            ansys_displacements = [[0, -0.037232	-0.32859	-0.99579	-2.0523	-3.4412	-5.1045	-7.0568	-9.3458	-12.019	-15.122]' ,...
                                   [0, 1.7645	6.6248	13.793	22.489	31.947	41.576	51.205	60.834	70.465	80.096]', ...
                                   [0, -0.78501	-3.1199	-6.9836	-12.364	-19.26	-27.678	-37.613	-49.055	-61.993	-76.414]'];

            ansys_solution = ansys_displacements + this.coordinates;
            
            h(end+1) = plot3(ansys_solution(:,1) , ansys_solution(:,2),  ansys_solution(:,3),  'g--');
            
            legend(h , {'Original structure' , 'Deformed structure' , 'Ansys solution'} , 'Location' , 'NorthWest');
            xlabel('x');  ylabel('y');  zlabel('z');
            axis([0 , max(this.coordinates(:,1)) + 20 , -50 , 100, -100, 20])
            grid on;  view([-39, 7]);

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
                       'disp(obj.results.displacements)'};
                error(strjoin(msg , '\n'))
            end
            
            % Reset object to input data
            this = FEM_3D_Solver;
            
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
            this.BC(11,:) = [NaN 1 NaN];                  % Imposed displacement at the far side
            
            % Loads
            this.loading_steps = 10;                      % Apply BC and loads in 10 steps
            this.loads = zeros(this.n_nodes , 3);         % Empty load matrix
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
            
            % Displacements obtained from ansys simulation
            ansys_displacements = [0,0;0.0389,0.010;0.077,0.040;0.115,0.087;0.15,0.153;0.188,0.237;0.222,0.341;0.253,0.466;...
                                   0.281,0.615;0.305,0.791;0.322,1;0.179,0.692;0.0225,0.361;-0.125,0.046;-0.238,-0.217;...
                                   -0.196,-0.346;-0.153,-0.393;-0.0935,0 ; -0.06762, 0; -0.046, 0; 0,0];


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
            
            % Check that constants are set
            this.checkConstants();

            % First step of intermediary solution: initial coordinates and
            % zero rotations.
            this.intermediary_solution = this.coordinates;
            cumulative_rotations  = zeros(this.n_nodes,1);
            cumdisp = zeros(this.n_nodes,3);

            % Apply loads and displacements in steps, and at each step
            % recalculate the stiffness matrix to account for beam
            % reorientations
            this.intermediary_bc = this.BC / this.loading_steps;
            this.intermediary_F  = this.loads / this.loading_steps;
            
            for j = 1 : this.loading_steps

                % Assemble matrix - use sparse for efficiency
                % this.stiffness_matrix = sparse(this.assemble_matrix());
                this.stiffness_matrix = this.assemble_matrix();

                % Loads must be recalculated at each step to account for
                % changed orientation of the elements
                F = this.assemble_loads();

                % Add imposed displacements as modified forces; see
                % reference within the function
                [this, F] = this.add_displacement_conditions(F);

                % Solve
                delta = this.stiffness_matrix \ F;


                % Compute displacements
                [node_disp , node_rot] = this.extract_nodal_displacements(delta);
                
                cumdisp = cumdisp + node_disp;
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
                p1 = this.intermediary_solution(this.connectivity(k,1),:);
                p2 = this.intermediary_solution(this.connectivity(k,2),:);

                % Constants
                L = norm(p1 - p2);
                E = this.E;
                Iy = this.Iy;
                Iz = this.Iz;
                A = this.CSA;
                G = this.G;
                J = this.J;
                

                % Matrix in local coordinates
                kl= [(E*A)/L        0             0              0        0             0          (-E*A)/L            0                  0              0             0              0;... 
                        0     (12*E*Iz)/L^3       0              0        0         (6*E*Iz)/L^2       0             (-12*E*Iz)/L^3       0              0             0         (6*E*Iz)/L^2;...
                        0           0       (12*E*Iy)/L^3        0   (-6*E*Iy)/L^2       0             0                0           (-12*E*Iy)/L^3       0        (-6*E*Iy)/L^2       0;...
                        0           0             0           (G*J)/L     0              0             0                0                 0           (-G*J)/L         0              0;...
                        0           0        (-6*E*Iy)/L^2       0     (4*E*Iy)/L        0             0                0            (6*E*Iy)/L^2        0         (2*E*Iy)/L         0;...
                        0      (6*E*Iz)/L^2       0              0        0         (4*E*Iz)/L         0           (-6*E*Iz)/L^2          0              0             0         (2*E*Iz)/L;...
                    (-E*A)/L        0             0              0        0              0           (E*A)/L            0                 0              0             0              0;... 
	                    0     (-12*E*Iz)/L^3      0              0        0        (-6*E*Iz)/L^2       0            (12*E*Iz)/L^3         0              0             0         (-6*E*Iz)/L^2;... 
                        0           0       (-12*E*Iy)/L^3       0   (6*E*Iy)/L^2        0             0                0            (12*E*Iy)/L^3       0         (6*E*Iy)/L^2       0;... 
                        0           0             0          (-G*J)/L     0              0             0                0                 0           (G*J)/L          0              0;... 
                        0           0       (-6*E*Iy)/L^2        0     (2*E*Iy)/L        0             0                0             (6*E*Iy)/L^2       0         (4*E*Iy)/L         0;...
                        0      (6*E*Iz)/L^2       0              0        0         (2*E*Iz)/L         0          (-6*E*Iz)/L^2            0             0             0         (4*E*Iz)/L;...
                    ];

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
                this = this.checkData();
            end
        end
        function this = set.BC(this, BC)
            % Stores values and checks for consistency
            % Order of boundary condition is 
            % [ux1; uy1; uz1; rotx1; roty1; rotz1; ux2; uy2; uz2; rotx2; roty2; rotz2;]

            this.BC  = BC;
            this.dof = this.count_DoF();   % Update Degrees of freedom
            if ~isempty(BC)
                this = this.checkData();
            end
        end
        function this = set.connectivity(this, connectivity)
            % Stores values and checks for consistency
            this.connectivity  = connectivity;
            if ~isempty(connectivity)
                this = this.checkData();
            end
        end
        function this = set.loads(this, loads)
            % Stores values and checks for consistency
            this.loads = loads;
            if ~isempty(loads)
                this = this.checkData();
            end
        end
        function loads = get.loads(this)
            % Returns the loads or an empty array for no loads
            
            loads = this.loads;
            % If no loads are assigned, assumed unloaded            
            if numel(loads) == 0
                loads = zeros(this.n_nodes , 6); 
            end
            

        end
        
         function h = plot_deformed(this)
            % Plots the structure given by coordinates and connectivity.
            % Returns handle to plot
            %   h = obj.plot_deformed();
            c = this.connectivity;
            
            h = zeros(size(c, 1) , 1);
            for k = 1 : size(c , 1)
                h(k) = plot3([this.results.deformed(c(k,1),1) ; this.results.deformed(c(k,2),1)] , ...
                             [this.results.deformed(c(k,1),2) ; this.results.deformed(c(k,2),2)] , ...
                             [this.results.deformed(c(k,1),3) ; this.results.deformed(c(k,2),3)], 'r-*');
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
                h(k) = plot3([this.coordinates(c(k,1),1) ; this.coordinates(c(k,2),1)] , ...
                             [this.coordinates(c(k,1),2) ; this.coordinates(c(k,2),2)] , ...
                             [this.coordinates(c(k,1),3) ; this.coordinates(c(k,2),3)], 'b-*');
                hold on;
            end
            axis equal;
            
        end
    end
    
    %% Private methods
    methods (Access = private)
        
        function checkConstants(this)
            
            constants = {'E', 'Young''s modulus';
                         'Iy', 'Moment of intertia';
                         'Iz', 'Moment of intertia';
                         'CSA', 'Cross-sectional area';
                         'G', 'Shear modulus';
                          'J', 'Torsional constant';};

            for k = 1 : size(constants, 1)
                if isempty(this.(constants{k,1})) || this.(constants{k,1}) < 0
                    error([constants{k,2} ' "' constants{k,1} '" should be set and higher than zero!'])
                end
            end
        end

        function this = checkData(this)
            % Checks for data consistency and raises errors
            
            errors = cell(0);
            if size(this.coordinates , 2) ~= 3
                errors{end+1,1} = 'Coordinates should be an Nx3 matrix';
                this.coordinates = [];
            end
            
            if ~isempty(this.BC) && size(this.BC , 2) ~= 6
                errors{end+1,1} = 'BC should be an Nx6 matrix; number of columns is wrong';
                this.BC = [];
            end
            
            if ~isempty(this.n_nodes) && this.n_nodes ~= 0 && ...
                    ~isempty(this.BC) && size(this.BC,1) ~= this.n_nodes
                errors{end+1,1} = 'BC should be an Nx6 matrix; number of rows is wrong';
                this.BC = [];
            end
                
            if ~isempty(this.loads) && size(this.loads , 2) ~= 6
                errors{end+1,1} = 'loads should be an Nx6 matrix; number of columns is wrong';
                this.loads = [];
            end
            if ~isempty(this.n_nodes) && this.n_nodes ~= 0 && ...
                    ~isempty(this.loads) && size(this.loads , 1) ~= this.n_nodes
                errors{end+1,1} = 'loads should be an Nx6 matrix; number of rows is wrong';
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

            N = numel(this.loads);

            F = zeros(N,1);
            for k = 1 : 6
                F(k:6:end) = this.intermediary_F(:,k);
            end

            
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
            % I tried a vectorized version of this loop, but it took longer
            % than the loop!
            %
            % Input:
            %   K: global stiffness matrix
            %   kg: local stiffess matrix in global coordinates
            %   g: steering vector (element numbers)

            eldof = 2*6; % Number of degrees of freedom per element 	
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
            
            
            % Rotation matrix
            n3 = (p2-p1)/norm((p2-p1));
            n1 = [1 0 0];    % Orientation of the coordinate system

            % n2(1) = n3(2)*n1(3)-n3(3)*n1(2);
            % n2(2) =-n1(3)*n3(1)+n1(1)*n3(3);
            % n2(3) = n3(1)*n1(2)-n1(1)*n3(2);
            % 
            % R=[n1; n2; n3];

            C = cross(n1, n3) ; 
            D = dot(n1, n3) ;
            NP0 = norm(n1);

            if all(C==0) % check for colinearity    

                R = sign(D) * (norm(n3) / NP0) ; % orientation and scaling

            else
                Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ; 
                R = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
            
                % Assemble in 12*12
                R=[  R     zeros(3) zeros(3) zeros(3);
                     zeros(3)   R     zeros(3) zeros(3);
                     zeros(3) zeros(3)   R     zeros(3);
                     zeros(3) zeros(3) zeros(3)   R    ];            
            end

            % Rotate local stiffness matris
            kg = R * kl * R';
        end

   

        function dof = count_DoF(this)
            % Counts and labels the degrees of freedom

            if isempty(this.n_nodes) || isempty(this.BC)
                dof = [];
                return
            end
            
            dof = ones(this.n_nodes, 6);  % Initialise the matrix nf to 1
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
            node_disp = output(:,1:3);    % X, Y and Z displacements
            node_rot  = output(:,4:6);                     % Rotations

        end        
        
        
    end
 
end
