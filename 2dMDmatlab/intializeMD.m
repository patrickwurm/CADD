function vel = intializeMD(vel, tset, kb, amass, mxatm, refpos, aymax, a)

%refpos, aymax and a can be deleted from the input parameters
%we only need these variables to set the velocities of some atoms to zero


disp(' ')
disp('Initializing MD Model ...')

sumvx=0;
sumvy=0;
sumv2x=0;
sumv2y=0;
sumv2=0;
sumvxcheck=0;
sumvycheck=0;

for ia=1:mxatm
    vel(ia,1) = (rand-0.5);
    vel(ia,2) = (rand-0.5);
    sumvx = sumvx + vel(ia,1);
    sumvy = sumvy + vel(ia,2);
end
sumvx=sumvx/mxatm;
sumvy=sumvy/mxatm;

for ia=1:mxatm
    vel(ia,1)=(vel(ia,1)-sumvx);
    vel(ia,2)=(vel(ia,2)-sumvy);
    sumvxcheck = sumvxcheck + vel(ia,1);
    sumvycheck = sumvycheck + vel(ia,2);
end

if 0>1
    sumv2 = 0;
    %set top atoms velocity to zero
    disp('-> I am setting the velocity of the top atoms to zero...')
    for ia=1:mxatm
        if refpos(ia,2)>aymax-3*a
            vel(ia,1) = 0;
            vel(ia,2) = 0;
        end
    end
end

for ia=1:mxatm
    sumv2 = sumv2 + vel(ia,1)^2+vel(ia,2)^2;
end

temp = amass/(kb*(2*mxatm-2))*sumv2;

disp(['-> Temperature before rescaling: ',num2str(temp),' K'])

fs = sqrt(2*kb*(2*mxatm-2)*tset/(amass*sumv2));
%the first factor 2 in the root is due to the fact that the
%kinetic energy will drop to about 1/2 of its initial value

sumv2=0;
for ia=1:mxatm
    vel(ia,1) = vel(ia,1)*fs;
    vel(ia,2) = vel(ia,2)*fs;
    sumv2 = sumv2 + vel(ia,1)^2+vel(ia,2)^2;
end

temp = amass/(kb*(2*mxatm-2))*sumv2;

disp(['-> Temperature after rescaling: ',num2str(temp),' K. (This should be 2x tset, as ~1/2 of the kinetic energy will transform to potential energy.)'])

end
