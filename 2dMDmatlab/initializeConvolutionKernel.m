function [kernel_matrix, pos_matrix] = initializeConvolutionKernel(M_window, band_atoms, band_atom_ass_atoms, fc, D, pos )
    
    M_window = round(1500/D); %Window size
    pos_matrix = zeros(length(band_atoms) ,2, M_window+1);
    kernel_matrix = zeros(size(pos_matrix));
    fc = 1.885e12/(2*pi)/(1e15/D);
    
    kernel = zeros(1,M_window+1);
    
    for i = 1:length(kernel)
        if i==M_window/2
            kernel(i) = 2*pi*fc;
        else
            kernel(i) = sin(2*pi*fc*(i-M_window/2))/(i-M_window/2);
        end
    end
    
    kernel = kernel'.*hamming(length(kernel));
    
    kernel = kernel/sum(kernel);
    
    for i = 1:length(band_atoms)
        kernel_matrix(i,:,:) = [kernel'; kernel'];
    end
    for i = 1:length(kernel)
        pos_matrix(:,:,i) = [pos(band_atoms,1)-pos(band_atom_ass_atoms,1), pos(band_atoms,2)-pos(band_atom_ass_atoms,2)];
    end

end

