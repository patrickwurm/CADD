function padpos = genpadatoms(a,axmin,aymin,axmax,aymax, rcut, example, dx, dy)

inboundsx=true;
inboundsy=true;

tolerance = 10e-12;

%set rcut larger, to to make sure, that 

padpos=[];
i=1;

if strcmp(example, 'Example1')
    xmin = axmin - floor(rcut/dx)*dx;
    ymin = aymin - floor(rcut/(dy/2))*dy/2;
    xmax = axmax + floor(rcut/dx)*dx;
    ymax = aymax;
    
    
    if mod(floor(rcut/(dy/2)),2)~=0
        ymin = ymin - dy/2;
    end
    
    
    for k=1:2
        if k==1
            xstart = xmin;
            ystart = ymin;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        else
            xstart = xmin+dx/2;
            ystart = ymin+dy/2;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        end
        while inboundsy
            while inboundsx
                if x<(axmin-tolerance) || y<(aymin-tolerance) || x>(axmax+tolerance)
                    padpos(i,1)=x;
                    padpos(i,2)=y;
                    i=i+1;
                end
                x = x+dx;
                if x>(xmax+tolerance)
                    inboundsx=false;
                end
            end
            x=xstart;
            y=y+dy;
            inboundsx=true;
            if y>ymax+tolerance
                inboundsy = false;
            end
        end
    end
elseif strcmp(example, 'Example2')
    %xmin = axmin - floor(rcut/a)*a;
    xmin = axmin - floor((2*rcut+1.5*dy)/(dx/2))*dx/2;
    ymin = aymin;
    xmax = axmin;
    ymax = aymax;
    
    if mod(floor(rcut/(dx/2)),2)~=0
        %ymin = ymin - sqrt(3)/2*a;
    end
    
    for k=1:2
        if k==1
            xstart = xmin;
            ystart = ymin+dy/2;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        else
            xstart = xmin+dx/2;
            ystart = ymin;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        end
        while inboundsy
            while inboundsx
                if x<(axmin-tolerance) || y<(aymin-tolerance) || x>(axmax+tolerance)
                    padpos(i,1)=x;
                    padpos(i,2)=y;
                    i=i+1;
                end
                x = x+dx;
                if x>(xmax+tolerance)
                    inboundsx=false;
                end
            end
            x=xstart;
            y=y+dy;
            inboundsx=true;
            if y>ymax+tolerance
                inboundsy = false;
            end
        end
    end
    
    
elseif strcmp(example, 'Dislocation')
    xmin = axmin - floor(rcut/(dx/2))*dx/2;
    ymin = aymin - floor(rcut/dy)*dy;
    xmax = axmax + floor(rcut/(dx/2))*dx/2;
    ymax = aymax;
    
    %do i need this?
    %     if mod(floor(rcut/(sqrt(3)/2*a)),2)~=0
    %         ymin = ymin - sqrt(3)/2*a;
    %     end
    
    
    for k=1:2
        if k==1
            xstart = xmin;
            ystart = ymin+dy/2;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        else
            xstart = xmin+dx/2;
            ystart = ymin;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        end
        while inboundsy
            while inboundsx
                if x<(axmin-tolerance) || y<(aymin-tolerance) || x>(axmax+tolerance)
                    padpos(i,1)=x;
                    padpos(i,2)=y;
                    i=i+1;
                end
                x = x+dx;
                if x>(xmax+tolerance)
                    inboundsx=false;
                end
            end
            x=xstart;
            y=y+dy;
            inboundsx=true;
            if y>ymax+tolerance
                inboundsy = false;
            end
        end
    end    

elseif strcmp(example, 'Radial_Pulse')
    xmin = axmin - floor(rcut/(dx/2))*dx/2;
    ymin = aymin - floor(rcut/dy)*dy;
    xmax = axmax + floor(rcut/(dx/2))*dx/2;
    ymax = aymax + floor(rcut/dy)*dy;
    
    for k=1:2
        if k==1
            xstart = xmin;
            ystart = ymin+dy/2;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        else
            xstart = xmin+dx/2;
            ystart = ymin;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        end
        while inboundsy
            while inboundsx
                if x<(axmin-tolerance) || y<(aymin-tolerance) || x>(axmax+tolerance) || y>(aymax+tolerance)
                    padpos(i,1)=x;
                    padpos(i,2)=y;
                    i=i+1;
                end
                x = x+dx;
                if x>(xmax+tolerance)
                    inboundsx=false;
                end
            end
            x=xstart;
            y=y+dy;
            inboundsx=true;
            if y>ymax+tolerance
                inboundsy = false;
            end
        end
    end    
    
elseif strcmp(example, 'Example2_rotated')
    xmin = axmin - floor((2*rcut+1.5*dx)/dx)*dx;
    ymin = aymin;
    xmax = axmin;
    ymax = aymax;
    
    if mod(floor(rcut/(dy/2)),2)~=0
        %ymin = ymin - sqrt(3)/2*a;
    end
    
        for k=1:2
            if k==1
                xstart = xmin;
                ystart = ymin;
                x = xstart;
                y = ystart;
                inboundsx=true;
                inboundsy=true;
            else
                xstart = xmin+dx/2;
                ystart = ymin+dy/2;
                x = xstart;
                y = ystart;
                inboundsx=true;
                inboundsy=true;
            end
            while inboundsy
                while inboundsx
                    if x<(axmin-tolerance) || y<(aymin-tolerance) || x>(axmax+tolerance)
                        padpos(i,1)=x;
                        padpos(i,2)=y;
                        i=i+1;
                    end
                    x = x+dx;
                    if x>(xmax+tolerance)
                        inboundsx=false;
                    end
                end
                x=xstart;
                y=y+dy;
                inboundsx=true;
                if y>ymax+tolerance
                    inboundsy = false;
                end
            end
        end
       
        
elseif strcmp(example, 'Tensile_Test')
    xmin = axmin - floor((2*rcut+1.5*dy)/(dx/2))*dx/2;
    ymin = aymin-5*dy;
    xmax = axmax + floor((2*rcut+1.5*dy)/(dx/2))*dx/2;
    ymax = aymax+5*dy;
    
    %do i need this?
    %     if mod(floor(rcut/(sqrt(3)/2*a)),2)~=0
    %         ymin = ymin - sqrt(3)/2*a;
    %     end
    
    
    for k=1:2
        if k==1
            xstart = xmin;
            ystart = ymin+dy/2;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        else
            xstart = xmin+dx/2;
            ystart = ymin;
            x = xstart;
            y = ystart;
            inboundsx=true;
            inboundsy=true;
        end
        while inboundsy
            while inboundsx
                if x<(axmin-tolerance) || x>(axmax+tolerance)
                    padpos(i,1)=x;
                    padpos(i,2)=y;
                    i=i+1;
                end
                x = x+dx;
                if x>(xmax+tolerance)
                    inboundsx=false;
                end
            end
            x=xstart;
            y=y+dy;
            inboundsx=true;
            if y>ymax+tolerance
                inboundsy = false;
            end
        end
    end
    
    elseif strcmp(example, 'Tensile_Test_rotated')
    xmin = axmin - floor((2*rcut+1.5*dx)/dx)*dx;
    ymin = aymin - 3*dy;
    xmax = axmax + floor((2*rcut+1.5*dx)/dx)*dx;
    ymax = aymax + 3*dy;
    
    if mod(floor(rcut/(dy/2)),2)~=0
        %ymin = ymin - sqrt(3)/2*a;
    end
    
        for k=1:2
            if k==1
                xstart = xmin;
                ystart = ymin;
                x = xstart;
                y = ystart;
                inboundsx=true;
                inboundsy=true;
            else
                xstart = xmin+dx/2;
                ystart = ymin+dy/2;
                x = xstart;
                y = ystart;
                inboundsx=true;
                inboundsy=true;
            end
            while inboundsy
                while inboundsx
                    if x<(axmin-tolerance) || x>(axmax+tolerance)
                        padpos(i,1)=x;
                        padpos(i,2)=y;
                        i=i+1;
                    end
                    x = x+dx;
                    if x>(xmax+tolerance)
                        inboundsx=false;
                    end
                end
                x=xstart;
                y=y+dy;
                inboundsx=true;
                if y>ymax+tolerance
                    inboundsy = false;
                end
            end
        end
       
    
else
    error('Example not implemented')
    
end

disp(['-> I generated ',num2str(length(padpos(:,1))),' pad atoms.'])


end