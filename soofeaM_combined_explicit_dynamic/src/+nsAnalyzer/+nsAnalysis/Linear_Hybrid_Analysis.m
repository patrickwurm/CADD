classdef Linear_Hybrid_Analysis < handle

    properties
        dynamic_analysis
        static_analysis
        example
        output_handler
        u
        u1
    end
    
    methods
        function self = Linear_Hybrid_Analysis(mymodel, output_handler, output_file_name)
            self.static_analysis = nsAnalyzer.nsAnalysis.Linear_Analysis(mymodel, output_handler, output_file_name);
            self.dynamic_analysis = nsAnalyzer.nsAnalysis.Linear_Dynamic_Analysis(mymodel, output_handler, output_file_name, 0);
            self.example = output_file_name;
            self.output_handler = self.dynamic_analysis.output_handler; %It doesnt matter which one we choose, output handler is the same for static and dynamic analysis
            self.u1 = self.dynamic_analysis.u1;
            self.u = self.dynamic_analysis.u; 
        end
        
        function setModel(self,model)
            self.dynamic_analysis.setModel(model);
            self.static_analysis.setModel(model);
        end
        
        function ComputeHMatrix(self, mxpatm, mxatm)
            self.dynamic_analysis.ComputeHMatrix(mxpatm, mxatm);
            self.static_analysis.ComputeHMatrix(mxpatm, mxatm);
        end
    end
end

