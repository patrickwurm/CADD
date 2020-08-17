classdef Tensile_Test_rotatedBCHandler < nsModel.BCHandler
    
    properties
        bottom_node_nr_list_L
        left_node_nr_list_L
        top_node_nr_list_L
        bottom_node_nr_list_R
        left_node_nr_list_R
        top_node_nr_list_R
        empty_boundary_cond_matrix
        u0L
        u0R
    end
    
    methods
        function self = Tensile_Test_rotatedBCHandler
            self.bottom_node_nr_list_L = [];
            self.left_node_nr_list_L = [];
            self.top_node_nr_list_L = [];
            self.bottom_node_nr_list_R = [];
            self.left_node_nr_list_R = [];
            self.top_node_nr_list_R = [];
            self.empty_boundary_cond_matrix = [];
            self.u0L = 0;
            self.u0R = 0;
        end
        
        function boundary_cond_matrix = incorporateBC(self, model, stateflag, timestep, settings)
            
            if nargin<5
                tstart = 0;
                tend = 0;
                disp = 0;
            else
                tstart = settings.continuum_external_change(1);
                tend = settings.continuum_external_change(2);
                disp = settings.continuum_external_displacement;  
            end
            
            %===========================
            %DISPLACEMENT BC
            %=========================== 
            
            if stateflag==0 % this is the first call!
                boundary_cond_matrix = [];
                
                if isempty( self.bottom_node_nr_list_L)
                    bottom_boundary = 3;
                    for node = model.boundary_dict(bottom_boundary).node_list
                         self.bottom_node_nr_list_L(end+1) =  node.number;
                    end
                end
                if isempty( self.bottom_node_nr_list_R)
                    bottom_boundary = 7;
                    for node = model.boundary_dict(bottom_boundary).node_list
                         self.bottom_node_nr_list_R(end+1) =  node.number;
                    end
                end
                
                if isempty( self.left_node_nr_list_L)
                    left_boundary = 2;
                    for node = model.boundary_dict(left_boundary).node_list
                         self.left_node_nr_list_L(end+1) = node.number;
                    end
                end
                if isempty( self.left_node_nr_list_R)
                    left_boundary = 6;
                    for node = model.boundary_dict(left_boundary).node_list
                         self.left_node_nr_list_R(end+1) = node.number;
                    end
                end
                
                if isempty( self.top_node_nr_list_L)
                    top_boundary = 4;
                    for node = model.boundary_dict(top_boundary).node_list
                         self.top_node_nr_list_L(end+1) = node.number;
                    end
                end
                if isempty( self.top_node_nr_list_R)
                    top_boundary = 8;
                    for node = model.boundary_dict(top_boundary).node_list
                         self.top_node_nr_list_R(end+1) = node.number;
                    end
                end
                self.empty_boundary_cond_matrix = zeros(length( self.bottom_node_nr_list_L),4);
            else
                boundary_cond_matrix = zeros(size(self.empty_boundary_cond_matrix));
            end
            
            i_boundary_cond = 0;
            for i = 1:length(  self.left_node_nr_list_L )
                for j=1:model.dimension
                    global_col_index = ( self.left_node_nr_list_L(i)-1)*model.dimension + j;
                    i_boundary_cond = i_boundary_cond+1;
                    boundary_cond_matrix(i_boundary_cond, 1) = global_col_index;
                    if j==1 && timestep >= tstart && timestep < tend
                        boundary_cond_matrix(i_boundary_cond, 2) = self.u0L -disp/(tend-tstart)*(timestep-tstart); %displacement
                        boundary_cond_matrix(i_boundary_cond, 3) = -disp/(tend-tstart); %velocity
                    elseif j==1 && timestep >= tend
                        boundary_cond_matrix(i_boundary_cond, 2) = self.u0L -disp; %displacement
                        boundary_cond_matrix(i_boundary_cond, 3) = 0; %velocity
                    elseif j==1
                        boundary_cond_matrix(i_boundary_cond, 2) = self.u0L; %displacement
                        boundary_cond_matrix(i_boundary_cond, 3) = 0; %velocity
                    else
                        boundary_cond_matrix(i_boundary_cond, 2) = 0; %displacement
                        boundary_cond_matrix(i_boundary_cond, 3) = 0; %velocity
                    end
                    boundary_cond_matrix(i_boundary_cond, 4) = 0; %force
                    boundary_cond_matrix(i_boundary_cond, 5) =  self.left_node_nr_list_L(i);
                    boundary_cond_matrix(i_boundary_cond,6) = j;
                end
            end
            
            for i = 1:length(  self.left_node_nr_list_R )
                for j=1:model.dimension
                    global_col_index = ( self.left_node_nr_list_R(i)-1)*model.dimension + j;
                    i_boundary_cond = i_boundary_cond+1;
                    boundary_cond_matrix(i_boundary_cond, 1) = global_col_index;
                    if j==1 && timestep >= tstart && timestep < tend
                        boundary_cond_matrix(i_boundary_cond, 2) = self.u0R + disp/(tend-tstart)*(timestep-tstart); %displacement
                        boundary_cond_matrix(i_boundary_cond, 3) = disp/(tend-tstart); %velocity
                    elseif j==1 && timestep >= tend
                        boundary_cond_matrix(i_boundary_cond, 2) = self.u0R + disp; %displacement
                        boundary_cond_matrix(i_boundary_cond, 3) = 0; %velocity
                    elseif j==1
                        boundary_cond_matrix(i_boundary_cond, 2) = self.u0R; %displacement
                        boundary_cond_matrix(i_boundary_cond, 3) = 0; %velocity
                    else
                        boundary_cond_matrix(i_boundary_cond, 2) = 0; %displacement
                        boundary_cond_matrix(i_boundary_cond, 3) = 0; %velocity
                    end
                    boundary_cond_matrix(i_boundary_cond, 4) = 0; %force
                    boundary_cond_matrix(i_boundary_cond, 5) =  self.left_node_nr_list_R(i);
                    boundary_cond_matrix(i_boundary_cond,6) = j;
                end
            end
            
            
            for i = 1:length(  self.bottom_node_nr_list_L )
                global_col_index = ( self.bottom_node_nr_list_L(i)-1)*model.dimension + 2;
                i_boundary_cond = i_boundary_cond+1;
                boundary_cond_matrix(i_boundary_cond, 1) = global_col_index;
                boundary_cond_matrix(i_boundary_cond, 2) = 0;
                boundary_cond_matrix(i_boundary_cond, 3) = 0;
                boundary_cond_matrix(i_boundary_cond, 4) = 0;
                boundary_cond_matrix(i_boundary_cond, 5) =  self.bottom_node_nr_list_L(i);
                boundary_cond_matrix(i_boundary_cond,6) = 2;
            end
            for i = 1:length(  self.bottom_node_nr_list_R )
                global_col_index = ( self.bottom_node_nr_list_R(i)-1)*model.dimension + 2;
                i_boundary_cond = i_boundary_cond+1;
                boundary_cond_matrix(i_boundary_cond, 1) = global_col_index;
                boundary_cond_matrix(i_boundary_cond, 2) = 0;
                boundary_cond_matrix(i_boundary_cond, 3) = 0;
                boundary_cond_matrix(i_boundary_cond, 4) = 0;
                boundary_cond_matrix(i_boundary_cond, 5) =  self.bottom_node_nr_list_R(i);
                boundary_cond_matrix(i_boundary_cond,6) = 2;
            end
            
            for i = 1:length(  self.top_node_nr_list_L )
                global_col_index = ( self.top_node_nr_list_L(i)-1)*model.dimension + 2;
                i_boundary_cond = i_boundary_cond+1;
                boundary_cond_matrix(i_boundary_cond, 1) = global_col_index;
                boundary_cond_matrix(i_boundary_cond, 2) = 0;
                boundary_cond_matrix(i_boundary_cond, 3) = 0;
                boundary_cond_matrix(i_boundary_cond, 4) = 0;
                boundary_cond_matrix(i_boundary_cond, 5) =  self.top_node_nr_list_L(i);
                boundary_cond_matrix(i_boundary_cond,6) = 2;
            end
            for i = 1:length(  self.top_node_nr_list_R )
                global_col_index = ( self.top_node_nr_list_R(i)-1)*model.dimension + 2;
                i_boundary_cond = i_boundary_cond+1;
                boundary_cond_matrix(i_boundary_cond, 1) = global_col_index;
                boundary_cond_matrix(i_boundary_cond, 2) = 0;
                boundary_cond_matrix(i_boundary_cond, 3) = 0;
                boundary_cond_matrix(i_boundary_cond, 4) = 0;
                boundary_cond_matrix(i_boundary_cond, 5) =  self.top_node_nr_list_R(i);
                boundary_cond_matrix(i_boundary_cond,6) = 2;
            end
            
            %===========================
            %NODAL FORCES
            %===========================
            
            % ...
            
            %===========================
            %SURFACE FORCES
            %===========================
            
            % ...
            
            %===========================
            %BODY FORCES
            %===========================
            
            % ...
        end
    end
end
