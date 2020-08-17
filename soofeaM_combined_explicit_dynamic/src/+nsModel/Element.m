classdef Element < nsModel.NodeContainer

    properties ( SetAccess = private )
        material
        pad_atoms = nsModel.PadAtom.empty
    end
    
    methods
        function self = Element(number, node_number_list)
            self@nsModel.NodeContainer(number, node_number_list);
        end
        
        function setMaterial( self, material )
            self.material = material;
        end
        
        function appendPadAtom(self, pad_atom_id, coordinates)
             self.pad_atoms(end+1) = nsModel.PadAtom(pad_atom_id, coordinates);
        end
        
        function resetPadAtoms(self)
            self.pad_atoms = nsModel.PadAtom.empty;
        end
    end
end

