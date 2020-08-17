function associatePadAtomswElements( fe_model, refpadpos, mxatm )

%loop over all pad atoms (their index starts with mxatm+1)
for i=1:length(refpadpos)
    %loop over all elements and find the element(s) to which the pad atoms
    %belongs to
    associated = false;
    for element = fe_model.element_dict
        if triangleTest(element.node_list, refpadpos(i,:))
            element.appendPadAtom(mxatm+i, refpadpos(i,:));
            %disp(['I associated Pad atom ',num2str(i), ' with element ',num2str(element.number),'.'])
            
            %why break? 1 association is enough to describe the motion of
            %the pad atom.
            %it makes no difference if the pad atom e.g. lies at a node and
            %is shared by 7 elements; if the motion is described by 1
            %element, it must be the same for all elements
            associated = true;
            break
        end
    end
    if associated == false
        error(['Pad atom ', num2str(i),' is not associated with any element.'])
    end
end


end

