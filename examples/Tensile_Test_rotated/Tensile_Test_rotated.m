function [mymodel, analysis] = Tensile_Test_rotated(outfile_path, results_path, cont_type, settings)

if isfield(settings.additional_input,'temperature')
    mesh_file_name = ['Tensile_Test_rotated_',num2str(settings.additional_input.temperature),'K.msh'];
else
    mesh_file_name = 'Tensile_Test_rotated_0.001K.msh';
end
mesh_file_path = ['examples/Tensile_Test_rotated/',mesh_file_name];
output_file_name = 'Tensile_Test_rotated';
output_file_path = outfile_path;

dimension = 2;
force_overwrite = 1; % has nothing to do with forces :p but with overwriting the output file.
input_handler = nsIO.GmshInputHandler(mesh_file_path, dimension);
output_handler = nsIO.VTKCADDOutputHandler(output_file_path, results_path, dimension, force_overwrite);

mymodel = nsModel.Model( dimension );
mymodel.setTimeStep(1e-15);
mymodel.addType(nsModel.nsType.ElementType(1,1,'triangle',1)); 
mymodel.type_dict{1}.setImplementation(nsAnalyzer.nsImplementation.LinearElementImpl());

mymodel.addType( nsModel.nsType.EdgeType(2,1,2) );
mymodel.type_dict{2}.setImplementation(nsAnalyzer.nsImplementation.BoundaryImpl());

% Add a material
% mymodel.addMaterial( nsModel.nsMaterial.LinearStVenantKirchhoffMaterial(1,2.1e11,0.3,'plane_stress',7900));

if isfield(settings.additional_input,'temperature')
    if settings.additional_input.temperature == 0.001
            c11 = 3.680606775836401*1e11*0.82; 
            c12 = 2.752730267678228*1e11*0.82;
            rho = 5.531505669578787e-07; %area density: mass/area
    elseif settings.additional_input.temperature == 50
            c11 = 3.693724287725043*1e11; 
            c12 = 2.830776388236943*1e11;
            rho = 5.529039934005879e-07;
    elseif settings.additional_input.temperature == 100
            c11 = 3.639419756641961*1e11;
            c12 = 2.988973201386569*1e11;
            rho = 5.526248399519090e-07;
    elseif settings.additional_input.temperature == 150
            c11 = 3.994240063602212*1e11;
            c12 = 2.711083926109033*1e11;
            rho = 5.522367135145907e-07;
    elseif settings.additional_input.temperature == 250
            c11 = 3.918409088922712*1e11; 
            c12 = 2.831420111106603*1e11;
            rho = 5.512796129720994e-07;
    else
        error('No material parameters saved for target temperature.')
    end
else
    c11 = 3.680606775836401*1e11*0.82; 
    c12 = 2.752730267678228*1e11*0.82;
    rho = 5.531505669578787e-07; %area density: mass/area
end

c11 = c11*1e-10; % -> in the computation of the virial stress, a depth of 1e-10 is assumed
c12 = c12*1e-10; % -> in the computation of the virial stress, a depth of 1e-10 is assumed

mymodel.addMaterial( nsModel.nsMaterial.LinearHexagonalMaterial(1,c11,c12,rho));

input_handler.read( mymodel )

mymodel.setBCHandler( Tensile_Test_rotatedBCHandler )
mymodel.bc_handler.incorporateBC(mymodel, 0, 0); %this is just a first call to initiate some quantities in the bc handler
if strcmp(cont_type,'dynamic')
    analysis = nsAnalyzer.nsAnalysis.Linear_Dynamic_Analysis(mymodel, output_handler, output_file_name);
elseif strcmp(cont_type,'static')
    analysis = nsAnalyzer.nsAnalysis.Linear_Analysis(mymodel, output_handler, output_file_name);
elseif strcmp(cont_type,'hybrid')
    analysis = nsAnalyzer.nsAnalysis.Linear_Hybrid_Analysis(mymodel, output_handler, output_file_name);
else
    error('Type of continuum not implemented (supported types: static, dynamic)')
end

end
