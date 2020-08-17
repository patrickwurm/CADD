classdef NumInt < handle
    
    properties
    end
    
    methods(Static)    
        % Diese Methode erstellt ein array vom Typ
        % nsNumeric.NaturalIntegrationPoint, das die Integrationspunkte des
        % gewünschten shape-types mit der gewünschten Anzahl an
        % Integrationspunkten (pro Koordinatenrichtung) enthält.
        function natural_ips = getNaturalIntegrationPoints(shape_type, number_of_int_points)
            ip_data = {};
            natural_ips = nsNumeric.NaturalIntegrationPoint.empty;
            
            %Be careful with number_of_int_points:
            %For shape_types linear, quad and hex, number_of_int_points is the number
            %of integration points PER DIRECTION
            %For shape_types triangle and teta, number_of_int_points is the
            %total number of integration points
            %Question: Is there a better/more natural way to handle this?
            
            if strcmp(shape_type,'linear') || strcmp(shape_type,'quad')  || strcmp(shape_type,'hex')
                if strcmp(shape_type,'linear')
                    ip_data = nsNumeric.IPData.getLinearIPData(number_of_int_points);
                    
                    for i=1:numel(ip_data)
                        ip=ip_data{i};
                        natural_ips(i) = nsNumeric.NaturalIntegrationPoint(ip{2},ip{1});
                    end 
                elseif strcmp(shape_type,'quad')
                    ip_data = nsNumeric.IPData.getLinearIPData(number_of_int_points);
                   
                    counter=0;
                    for i=1:numel(ip_data)
                        for j=1:numel(ip_data)
                            counter = counter+1;
                            natural_ips(counter) = nsNumeric.NaturalIntegrationPoint(ip_data{i}{2}*ip_data{j}{2},[ip_data{i}{1},ip_data{j}{1}]);
                        end
                    end
                elseif strcmp(shape_type,'hex')
                    ip_data = nsNumeric.IPData.getLinearIPData(number_of_int_points);
                    
                    counter=0;
                    for i=1:numel(ip_data)
                        for j=1:numel(ip_data)
                            for k=1:numel(ip_data)
                            counter = counter+1;
                            natural_ips(counter) = nsNumeric.NaturalIntegrationPoint(ip_data{i}{2}*ip_data{j}{2}*ip_data{k}{2},[ip_data{i}{1},ip_data{j}{1},ip_data{k}{1}]);
                            end
                        end
                    end
                end
            elseif strcmp(shape_type,'triangle') || strcmp(shape_type,'tetra')
                if strcmp(shape_type,'triangle')
                    ip_data = nsNumeric.IPData.getTriangleIPData(number_of_int_points);
                elseif strcmp(shape_type,'tetra')
                    ip_data = nsNumeric.IPData.getTetraIPData(number_of_int_points);
                end
                for i=1:numel(ip_data)
                    ip=ip_data{i};
                    natural_ips(i) = nsNumeric.NaturalIntegrationPoint(ip{2},ip{1});
                end
            else
                error(['No integration points for the given shape type ', shape_type, ' implemented'])
            end
        end

        % Diese Methode wird ausschließlich zu Übungszwecken verwendet.
        function ival = integrate(func, int_point_array)
            ival = 0;
            for ip = int_point_array 
                ival = ival + func(ip.natural_coordinates) * ip.weight();
            end
        end
        
        % Mit Hilfe dieser Methode werden die Integrationen in soofeaM
        % durchgeführt.
        function ival = methodIntegrate(func, int_point_array, parameters)
            ival = 0;
            for ip = int_point_array 
                ival = ival + func(ip, parameters) * ip.getWeight();
            end
        end
    end    
end
