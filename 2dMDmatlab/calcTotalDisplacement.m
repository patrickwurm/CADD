function [total_displacement, max_displacement_x] = calcTotalDisplacement(pos,refpos)

u = pos-refpos; 

total_displacement = sum(sum(abs(u)));

max_displacement_x = max(u(:,1))+abs(min(u(:,1)));


end

