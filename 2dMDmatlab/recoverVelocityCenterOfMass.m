function vel = recoverVelocityCenterOfMass(vel, mxatm)

sumvx=0;
sumvy=0;
sumvxcheck=0;
sumvycheck=0;
for ia=1:mxatm
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

end

