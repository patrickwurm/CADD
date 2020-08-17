function [ trigger, filter_output ] = Cont_trigger_DSP(trigger_case, state, pos_matrix, kernel_matrix, interface_pos, run_in)

switch trigger_case
    case 7
        filter_output = dot(kernel_matrix,pos_matrix,3);
        if state ~= 0
            %the ",1" is to improve performance.
            %it only returns the first relevant index and then stops
            
             a = length(find(filter_output > run_in.max_distance));
             b = length(find(filter_output < run_in.min_distance));
             term_out = sprintf('Exceeded EQ fluc (max) in %d pair(s), exceeded EQ fluc (min) in %d pair(s)',a,b);
             disp(term_out)
            
            if length(find(filter_output > run_in.max_distance,1))==1 || length(find(filter_output < run_in.min_distance,1))==1
                trigger = 1;
            else
                trigger = 0;
            end
        else
            trigger = 1;
        end
    case 0 % Always invoke continuum
        trigger = 1;
        filter_output = [];
    case 100 % Never invoke continuum
        trigger = 0;
        filter_output = [];
end
end
