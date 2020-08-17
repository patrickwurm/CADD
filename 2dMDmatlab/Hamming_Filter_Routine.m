function [ filter_output_hist ] = Hamming_Filter_Routine( input_signal, M_window, cutoff_freq, sample_freq )

signal_length = length(input_signal);

disp_vector = zeros(1, M_window+1);

fc = cutoff_freq/sample_freq;

kernel = zeros(size(disp_vector));

for i = 1:length(kernel)
    if i==M_window/2
        kernel(i) = 2*pi*fc;
    else
        kernel(i) = sin(2*pi*fc*(i-M_window/2))/(i-M_window/2);
    end
end

kernel = kernel'.*hamming(length(kernel));

kernel = kernel/sum(kernel);

filter_output_hist = zeros(1,signal_length);

for i = 1:signal_length
    
    disp_vector = cat(2, input_signal(i), disp_vector(1:end-1));
    
    output_signal = dot(kernel,disp_vector);
    
    filter_output_hist(i) = output_signal;
end

end

