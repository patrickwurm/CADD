function [elements] = CreateDetectionBand(distance, a)


%This function returns the node numbers of the triangles in the detection
%band
%elements has the triangles in the rows and the nodes of every triangle in
%the columns
%This function is very problem specific, so use it with care

global inactive_atoms
global refpos

postol = 10^-3*a;

atoms = [];

left_end = min(refpos(inactive_atoms,1))-2*sqrt(3)*a-postol;
right_end = max(refpos(inactive_atoms,1))+2*sqrt(3)*a+postol;
top_end = distance+postol;
bottom_end = distance-a-postol;

for i=1:length(refpos(:,1))
    if refpos(i,1) >= left_end && refpos(i,1) <= right_end && refpos(i,2) >= bottom_end && refpos(i,2) <= top_end
        atoms(end+1) = i;
    end
end

%group atoms into top/middle/bottom layer
atoms_top=[];
atoms_middle=[];
atoms_bottom=[];

for i=atoms
    if refpos(i,2)>=distance-postol
        atoms_top(end+1) = i;
    elseif refpos(i,2)>=distance-a/2-postol
        atoms_middle(end+1) = i;
    else
        atoms_bottom(end+1) = i;
    end
end

%sort atoms in each group by their x position

[~, order] = sort(refpos(atoms_top,1));
atoms_top = atoms_top(order);
[~, order] = sort(refpos(atoms_middle,1));
atoms_middle = atoms_middle(order);
[~, order] = sort(refpos(atoms_bottom,1));
atoms_bottom = atoms_bottom(order);

%determine type of left end of the band
if refpos(atoms_middle(1),1) < refpos(atoms_top(1),1)
    left_end = 1;
else
    left_end = 2;
end

%determine type of right end of the band
if refpos(atoms_middle(end),1) > refpos(atoms_top(end),1)
    right_end = 1;
else
    right_end = 2;
end

if left_end == 1
    if right_end == 2
        last_index_right = length(atoms_top);
        for i = 1:last_index_right
            elements(i,:) = [atoms_middle(i), atoms_top(i), atoms_bottom(i)];
        end
        for i = 1:last_index_right-1
            elements(end+1,:) = [atoms_top(i), atoms_middle(i+1), atoms_bottom(i)];
        end
    else
        last_index_right = length(atoms_top);
        for i = 1:last_index_right
            elements(i,:) = [atoms_middle(i), atoms_top(i), atoms_bottom(i)];
        end
        for i = 1:last_index_right
            elements(end+1,:) = [atoms_top(i), atoms_middle(i+1), atoms_bottom(i)];
        end
    end
%     for i = 1:last_index_right
%         elements(i,:) = [atoms_middle(i), atoms_top(i), atoms_bottom(i)];
%     end
%     for i = 1:last_index_right-1
%         elements(end+1,:) = [atoms_top(i), atoms_middle(i+1), atoms_bottom(i)];
%     end
else
    if right_end == 2
        last_index_right = length(atoms_top)-1;
        for i = 1:last_index_right
            elements(i,:) = [atoms_top(i), atoms_middle(i), atoms_bottom(i)];
        end
        for i = 1:last_index_right
            elements(end+1,:) = [atoms_middle(i), atoms_top(i+1), atoms_bottom(i+1)];
        end
    else
        last_index_right = length(atoms_top);
        for i = 1:last_index_right
            elements(i,:) = [atoms_top(i), atoms_middle(i), atoms_bottom(i)];
        end
        for i = 1:last_index_right-1
            elements(end+1,:) = [atoms_middle(i), atoms_top(i+1), atoms_bottom(i+1)];
        end
    end

end

%TEST THE FUNCTION
figure(10)
for i=1:length(elements(:,1))
    scatter(refpos(elements(i,1),1),refpos(elements(i,1),2),'r')
    hold on
    scatter(refpos(elements(i,2),1),refpos(elements(i,2),2),'g')
    scatter(refpos(elements(i,3),1),refpos(elements(i,3),2),'b')
end

end

