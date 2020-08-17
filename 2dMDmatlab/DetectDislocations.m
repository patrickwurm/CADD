function [ detected_in_el, burgers ] = DetectDislocations( reference_vectors, pos, elements, a, correction)

%"elements" is a matrix which contains the elements in its rows
%and the atoms of each element in its columns (the atoms should be sorted
%clockwise)

%"reference_vectors" is a matrix which contains the reference_vectors as
%its columns

%"detected" should be something like the element number in which a
%dislocation was detected

% reference_vectors = [ 0 ,sqrt(3)/2*a, sqrt(3)/2*a, 0, -sqrt(3)/2*a, -sqrt(3)/2*a; ...
%                      a, 1/2*a, -1/2*a, -a, -1/2*a];

number_of_elements = length(elements(:,1));

test_matrix = zeros(size(reference_vectors));

detected_in_el = [];
burgers = [];

best_vec = zeros(2,3);

for i = 1:number_of_elements
    p1 = pos(elements(i,1),:);
    p2 = pos(elements(i,2),:);
    p3 = pos(elements(i,3),:);
    
    if correction(i,1)~=0
        if correction(i,1)==1
            p1 = p1 + correction(i,2:3);
        elseif correction(i,1)==2
            p2 = p2 + correction(i,2:3);
        else
            p3 = p3 + correction(i,2:3);
        end
    end
    
    for j = 1:3
        
        if j==1
            test_vec = p2 - p1;
        elseif j==2
            test_vec = p3 - p2;
        else
            test_vec = p1 - p3;
        end
        
        test_vectors = repmat(test_vec', 1, length(reference_vectors(1,:)));
        norms = sqrt(sum((test_vectors-reference_vectors).^2,1));
        [~, min_index] = min(norms);
        best_vec(:,j) = reference_vectors(:,min_index);
    end
    sum_vec = sum(best_vec,2);
    if norm(sum_vec) >= a*10^-3
        detected_in_el = [detected_in_el, i];
        burgers = [burgers, sum_vec];
    end
end

end

