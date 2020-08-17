close all
clear all
a = 2.929295*1e-10;

axmin = 0;
aymin = 0;
axmax = 17*sqrt(3)*a;
aymax = 40*a;

delta = 15*a;

pos = Volterra_Cut_Test(a, axmin, aymin, axmax, aymax, 0);

disl_burgers = -a;
original_disl_pos = [(axmax-axmin)/2-a*sqrt(3)/2*1/3, 10*a];
new_disl_pos = original_disl_pos + [0, -delta];

scatter(pos(:,1), pos(:,2))
axis equal
hold on
scatter(original_disl_pos(1), original_disl_pos(2))
scatter(new_disl_pos(1), new_disl_pos(2))

displacement = zeros(length(pos(:,1)),2);

incl_original = 0;
incl_ghost = 1;
incl_new = 1;

needleman = 1;
hirth = 0;

if incl_original
for i = 1:length(pos(:,1))
    %original dislocation
    
    b = -disl_burgers;
    % we gotta get Poissons number from the FE model
    % (somehow)
    nu = 0.3243; %Poissons number
    
    dx2 = pos(i,1) - original_disl_pos(1);
    dx1 = pos(i,2) - original_disl_pos(2);
    %dx1 = -dx1;
    %dx2 = -dx2;
    
    if abs(dx1) < a/100 && abs(dx2) < a/100;
        displacement = 0;
    else
        if needleman
        % Needleman: wrong??
        displacement_y = b/(2*pi*(1-nu)) * ( 1/2*(dx1.*dx2)/(dx1.^2+dx2.^2) - (1-nu)*atan(dx1./dx2) );
        displacement_x = b/(2*pi*(1-nu)) * ( 1/2*(dx2.^2)/(dx1.^2+dx2.^2) - 1/4*(1-2*nu)*log((dx1.^2+dx2.^2)/b^2) );
        else
        displacement_y = b/(2*pi) * ( 1/2*(dx1.*dx2)/((1-nu)*(dx1.^2+dx2.^2)) + atan(dx2./dx1) );
        displacement_x = -b/(2*pi) * ( 1/(4*(1-nu))*(1-2*nu)*log((dx1.^2+dx2.^2)) + (dx1^2-dx2^2)/(4*(1-nu)*(dx1.^2+dx2.^2)) );
        end
        displacement(i,:) = displacement(i,:) + [displacement_x, displacement_y];
    end
    
end
end

if incl_ghost
for i = 1:length(pos(:,1))
    %ghost dislocation
    
    b = -disl_burgers;
    % we gotta get Poissons number from the FE model
    % (somehow)
    nu = 0.3243; %Poissons number
    
    dx2 = pos(i,1) - original_disl_pos(1);
    dx1 = pos(i,2) - original_disl_pos(2);
    %dx1 = -dx1;
    %dx2 = -dx2;
    
    if abs(dx1) < a/100 && abs(dx2) < a/100;
        displacement = 0;
    else
        if needleman
        % Needleman: wrong??
        displacement_y = b/(2*pi*(1-nu)) * ( 1/2*(dx1.*dx2)/(dx1.^2+dx2.^2) - (1-nu)*atan(dx1./dx2) );
        displacement_x = b/(2*pi*(1-nu)) * ( 1/2*(dx2.^2)/(dx1.^2+dx2.^2) - 1/4*(1-2*nu)*log((dx1.^2+dx2.^2)/b^2) );
        else
        displacement_y = b/(2*pi) * ( 1/2*(dx1.*dx2)/((1-nu)*(dx1.^2+dx2.^2)) + atan(dx2./dx1) );
        displacement_x = -b/(2*pi) * ( 1/(4*(1-nu))*(1-2*nu)*log((dx1.^2+dx2.^2)) + (dx1^2-dx2^2)/(4*(1-nu)*(dx1.^2+dx2.^2)) );
        end
        displacement(i,:) = displacement(i,:) + [displacement_x, displacement_y];
    end
    
end
end

if incl_new
for i = 1:length(pos(:,1))
    %new dislocation
    
    b = disl_burgers;
    % we gotta get Poissons number from the FE model
    % (somehow)
    nu = 0.3243; %Poissons number
    
    dx2 = pos(i,1) - new_disl_pos(1);
    dx1 = pos(i,2) - new_disl_pos(2);
    %dx1 = -dx1;
    %dx2 = -dx2;
    
    if abs(dx1) < a/100 && abs(dx2) < a/100;
        displacement = 0;
    else
        if needleman
        % Needleman: wrong??
        displacement_y = b/(2*pi*(1-nu)) * ( 1/2*(dx1.*dx2)/(dx1.^2+dx2.^2) - (1-nu)*atan(dx1./dx2) );
        displacement_x = b/(2*pi*(1-nu)) * ( 1/2*(dx2.^2)/(dx1.^2+dx2.^2) - 1/4*(1-2*nu)*log((dx1.^2+dx2.^2)/b^2) );
        else
        displacement_y = b/(2*pi) * ( 1/2*(dx1.*dx2)/((1-nu)*(dx1.^2+dx2.^2)) + atan(dx2./dx1) );
        displacement_x = -b/(2*pi) * ( 1/(4*(1-nu))*(1-2*nu)*log((dx1.^2+dx2.^2)) + (dx1^2-dx2^2)/(4*(1-nu)*(dx1.^2+dx2.^2)) );
        end
        displacement(i,:) = displacement(i,:) + [displacement_x, displacement_y];
    end
    
end
end



quiver(pos(:,1),pos(:,2),displacement(:,1),displacement(:,2));

pos_new =  pos +displacement;

figure(33)
scatter(pos(:,1),pos(:,2))
hold on
axis equal
scatter(pos_new(:,1), pos_new(:,2))


