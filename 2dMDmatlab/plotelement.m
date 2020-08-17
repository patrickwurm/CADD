function plotelement( element, pos, natural_pos, pad_atom_id )

close all;

a = element.node_list(1).coordinates;
b = element.node_list(2).coordinates;
c = element.node_list(3).coordinates;

if isempty(element.pad_atoms)
    disp('No pad atoms associated')
else
    p = pos(pad_atom_id,:);
    np = natural_pos;
end

figure(72)
line([a(1),b(1)], [a(2), b(2)], 'Color', 'r');
hold on;
line([b(1),c(1)], [b(2), c(2)], 'Color', 'r');
line([a(1),c(1)], [a(2), c(2)], 'Color', 'r');
plot(p(1),p(2),'k^');

figure(73)
line([0,1], [0, 0], 'Color', 'r');
hold on;
line([1,0], [0, 1], 'Color', 'r');
line([0,0], [0, 1], 'Color', 'r');
plot(np(1),np(2),'k^');

end

