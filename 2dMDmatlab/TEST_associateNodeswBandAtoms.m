function [ ] = TEST_associateNodeswBandAtoms(band_atom_ass_atoms, band_atoms, pos, interface_atom_list)

figure(30)
hold on

j=0;
for i = band_atoms
    j=j+1
    scatter(pos(i,1),pos(i,2))
    scatter(pos(band_atom_ass_atoms(j),1), pos(band_atom_ass_atoms(j),2))
end

end

