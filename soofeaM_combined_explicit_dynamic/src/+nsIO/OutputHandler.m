classdef ( Abstract ) OutputHandler < handle
    %OUTPUTHANDLER 
    
    properties ( SetAccess = protected )
        file_name
        folder_name
        dimension
        force_overwrite
    end
    
    methods
        function self = OutputHandler( file_name, folder_name, dimension, force_overwrite )
            self.file_name = file_name;
            self.folder_name = folder_name;
            self.dimension = dimension;
            self.force_overwrite = force_overwrite;
        end
        
        function write( self, model, vtks , atomrefpos, atompos, mxatm, mxpatm, step, energies, atom_velocities, nodal_velocities )
            
            current_file_path = [self.file_name,'.vtk.',num2str(vtks)];
            current_folder_path = fileparts(current_file_path);
            
            if exist(current_folder_path,'dir') == 7
            else
                mkdir (current_folder_path);
            end
                      
            if exist(current_file_path,'file')==2
                disp(['-> The file ',current_file_path, ' already exists.'])
                if self.force_overwrite == true
                    disp(['--> I am overwriting all vtk files with that name.'])
                    delnames = [self.file_name,'.vtk.*'];
                    delete (delnames)
                else
                    %Here goes some kind of dialog like "Do you want to
                    %delete the file?"
                    %Needs to be implemented
                end
            end
            
            standard_file_path = [self.file_name,'_standard'];
            
            if ~exist(standard_file_path,'file')
                write_standard( self, model, vtks , atomrefpos, atompos, mxatm, mxpatm, step )
            else
                copyfile(standard_file_path, current_file_path)
            end

            file_id = fopen(current_file_path,'a');

            % write output to file
            self.writeResults( model, file_id, atomrefpos, atompos, mxatm, mxpatm, step, energies, atom_velocities, nodal_velocities );
            
            fclose(file_id);
        end
        
        function write_standard( self, model, vtks, atomrefpos, atompos, mxatm, mxpatm, step )
            
            current_file_path = [self.file_name,'_standard'];
            current_folder_path = fileparts(current_file_path);
            
            if exist(current_folder_path,'dir') == 7
            else
                mkdir (current_folder_path);
            end
            
            if exist(current_file_path,'file')==2
                disp(['-> The file ',current_file_path, ' already exists.'])
                if self.force_overwrite == true
                    disp(['--> I am overwriting all vtk files with that name.'])
                    delnames = [self.file_name,'.vtk.*'];
                    delete (delnames)
                else
                    %Here goes some kind of dialog like "Do you want to
                    %delete the file?"
                    %Needs to be implemented
                end
            end
            
            % write output to file
            file_id = fopen(current_file_path,'w');
            
            self.writeResults_standard( model, file_id, atomrefpos, atompos, mxatm, mxpatm, step );
            
            fclose(file_id);
        end
        
    end
    
    methods ( Abstract )
        writeResults ( self, model, file_id )
    end
end

