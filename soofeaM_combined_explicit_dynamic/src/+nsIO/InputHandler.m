classdef InputHandler < handle
    %INPUTHANDLER 
    
    properties ( SetAccess = protected )
        file_name
        dimension
    end
    
    methods
        function self = InputHandler( file_name, dimension )
            self.file_name = file_name;
            self.dimension = dimension;
        end
        
        function read( self, model )
            % open file
            file_id = fopen(self.file_name);
            
            % read mesh
            self.readMesh( model, file_id );
            
            %close file
            fclose(file_id);
        end
    end
    
    methods ( Abstract )
        readMesh ( self, model, file_id )
    end
end

