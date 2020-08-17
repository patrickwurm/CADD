classdef LinearIPData
    
    properties
    end
    
    methods(Static)
        function ip_data = getIpData(number_of_int_points)
            if number_of_int_points == 1
                ip_data = {{0, 2}};
            elseif number_of_int_points == 2
                ip_data = {{-0.577350269189626, 1},...
                    {0.577350269189626, 1}};
            elseif number_of_int_points == 3
                ip_data = {{-0.774596669241483, 0.555555555555556},...
                    {0, 0.888888888888889},...
                    {0.774596669241483, 0.555555555555556}};
            elseif number_of_int_points == 4
                ip_data = {{-0.861136311594953, 0.347854845137454},...
                    {-0.339981043584856, 0.652145154862546},...
                    {0.339981043584856, 0.652145154862546},...
                    {0.861136311594953, 0.347854845137454}};
            else
                error('number of integration points not implemented')
            end
        end
    end
end

