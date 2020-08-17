classdef Analysis < handle
    properties( Access = protected )
        model;
        model_progression;
        bc_handler;
    end
    
    methods
        function self = Analysis( model, model_progression )
            self.model = model;
            self.model_progression = model_progression;
            self.bc_handler = analyze.BCHandler(model,model_progression);
        end
    end
    
    methods( Abstract = true )
        run(self)
    end
    
    methods( Access = protected )
        function [global_stiffness, load_vector] = ...
                assemble_k_fint__( self, global_stiffness, load_vector, iter )
            % disp('> Assemble global stiffness and internal forces');
            global global_data;
            dim = global_data.dimension;
            
            load_vector = zeros(self.model.nodeAmount()*dim,1);
            global_stiffness = sparse(self.model.nodeAmount()*dim, ...
                self.model.nodeAmount()*dim);

            %
            for element_counter = 1:1:self.model.elementAmount()
                
                element = self.model.getElement(element_counter);
                element_impl = element.element_type.implementation;
                [k_element, f_int] = element_impl.calcElementMatrices( element );
                
                for counter = 1:length(self.model.getElement(1).node_array)
                    node_number = ...
                        element.node_array(counter).number;
                    for d=1:dim
                        index = (node_number - 1)*dim + d;
                        load_vector(index) = ...
                            load_vector(index) - f_int(d,counter);
                    end
                end
                
                for row_counter = 1:length(k_element)
                    row_node_number = ...
                        element.node_array(ceil(row_counter/dim)).number;
                    row_coord_offset = ...
                        mod(row_counter-1,dim)+1;
                    row_index = (row_node_number - 1)*dim + row_coord_offset;
                    for col_counter = 1:length(k_element)
                        col_node_number = ...
                            element.node_array(ceil(col_counter/dim)).number;
                        col_coord_offset = ...
                            mod(col_counter-1,dim)+1;
                        col_index = (col_node_number - 1)*dim + col_coord_offset;
                        
                        node = element.node_array(ceil(col_counter/dim));
                        if( node.getDisplacementDOF(). ...
                                getDOFValue(col_coord_offset).constraint )
                            constr_value = node.getDisplacementDOF().getDOFValue(col_coord_offset).dof_value;
                            load_vector(row_index) = ...
                                load_vector(row_index) - constr_value * k_element(row_counter,col_counter);
                            
                            if( row_counter == col_counter )
                                global_stiffness(row_index,col_index) = ...
                                    global_stiffness(row_index,col_index) + k_element(row_counter,col_counter);
                            end
                        else
                            global_stiffness(row_index,col_index) = ...
                                global_stiffness(row_index,col_index) + k_element(row_counter,col_counter);
                        end
                    end
                end

            end
        end
        
        function [global_stiffness, load_vector] = ...
                assemble_fext__( self, global_stiffness, load_vector )
            % disp('> Assemble load vector');
            global global_data;
            dim = global_data.dimension;
            
            for node_counter = 1:self.model.nodeAmount()
                force_load = self.model.getNode(node_counter).force_load;
                if( ~ isempty(force_load) )
                    for coord_offset = 1:dim
                        index = (node_counter-1)*dim + coord_offset;
                        load_vector( index ) = ...
                            load_vector( index ) + force_load.getLoadValue(coord_offset);
                    end
                end
            end
        end   
        
        function updateDOF__(self, solution_vector)
            global global_data;
            dim = global_data.dimension;
            
            for node_counter = 1:self.model.nodeAmount()
                for dim_counter = 1:dim
                    index = (node_counter-1)*dim + dim_counter;
                    dof_value = self.model.getNode(node_counter).getDisplacementDOF().getDOFValue(dim_counter);
                    if( ~ dof_value.constraint )
                        dof_value.setValue( solution_vector(index) );
                    end
                end
            end
        end
        
        function unbalanced_energy = unbalanced_energy__(self, solution_vector, load_vector)
            global global_data;
            dim = global_data.dimension;
            unbalanced_energy = 0.0;
            for node_no = 1:self.model.nodeAmount()
                for dim_no = 1:dim
                    node = self.model.getNode(node_no);
                    if( ~node.getDisplacementDOF(). ...
                            getDOFValue(dim_no).constraint )
                        unbalanced_energy = unbalanced_energy + ...
                            abs( solution_vector( (node_no-1)*dim + dim_no ) * load_vector((node_no-1)*dim + dim_no ) );
                    end
                end
            end
        end
        
    end
end












