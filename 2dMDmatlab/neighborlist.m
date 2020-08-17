function [ nl, numneigh ] = neighborlist(pos, mxatm, mxpatm, rcut, pbc, xbox, ybox, a, example, settings)
%Note that there are also neighborlist entries for the pad atoms.
%This is needed for the force evaluation.
%See "getforce".

counter =0;
nl=[];
numneigh = zeros(mxatm+mxpatm,1);

global refpos
global refpadpos

rneigh = settings.effective_cutoff; %"effective" cutoff radius

if strcmp(example, 'Example1')

for ia=1:length(pos)-1
    for ja=ia+1:length(pos)
        sx = (pos(ia,1)-pos(ja,1));
        sy = (pos(ia,2)-pos(ja,2));
        %PBC in y-direction (only valid in Example2)
        if pbc
            sy = sy-ybox * round(sy/ybox);
        end
        rrr2 = sx^2 + sy^2;
        if rrr2<=rneigh^2
            numneigh(ia) = numneigh(ia)+1;
            nl(ia,numneigh(ia))=ja;
            numneigh(ja) = numneigh(ja)+1;
            nl(ja,numneigh(ja))=ia;
        end
    end
end

elseif strcmp(example, 'Example2') || strcmp(example, 'Example2_rotated') || strcmp(example, 'Tensile_Test') || strcmp(example, 'Tensile_Test_rotated')
    
    global exclude_atoms
    
    atoms = 1:mxatm+mxpatm;
    nl=zeros(mxatm+mxpatm,40);
    numneigh = zeros(mxatm+mxpatm,1);
    
    atoms(find(ismember(atoms,exclude_atoms)))=[];
    for ia=1:length(atoms)-1
            for ja=ia+1:length(atoms)
                    sx = (pos(atoms(ia),1)-pos(atoms(ja),1));
                    sy = (pos(atoms(ia),2)-pos(atoms(ja),2));
                    %PBC in y-direction (only valid in Example2)
                    if pbc
                        sy = sy-ybox * round(sy/ybox);
                    end
                    rrr2 = sx^2 + sy^2;
                    if rrr2<=rneigh^2
                        numneigh(atoms(ia)) = numneigh(atoms(ia))+1;
                        nl(atoms(ia),numneigh(atoms(ia)))=atoms(ja);
                        numneigh(atoms(ja)) = numneigh(atoms(ja))+1;
                        nl(atoms(ja),numneigh(atoms(ja)))=atoms(ia);
                    end
            end
    end

elseif strcmp(example, 'Dislocation') || strcmp(example, 'Radial_Pulse')
    
    for ia=1:length(pos)-1
        for ja=ia+1:length(pos)
            sx = (pos(ia,1)-pos(ja,1));
            sy = (pos(ia,2)-pos(ja,2));
            rrr2 = sx^2 + sy^2;
            if rrr2<=rneigh^2
                numneigh(ia) = numneigh(ia)+1;
                nl(ia,numneigh(ia))=ja;
                numneigh(ja) = numneigh(ja)+1;
                nl(ja,numneigh(ja))=ia;
            end
        end
    end

else
    error('Example not found')
end

end