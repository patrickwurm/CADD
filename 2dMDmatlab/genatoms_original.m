function pos = genatoms_original(a,axmin,aymin,axmax,aymax, dx, dy)

inboundsx=true;
inboundsy=true;

tolerance = 10e-12;

pos=[];
i=1;


for k=1:2
    if k==1
        xstart = axmin;
        ystart = aymin;
        x = xstart;
        y = ystart;
        inboundsx=true;
        inboundsy=true;
    else
        xstart = axmin+dx/2;
        ystart = aymin+dy/2;
        x = xstart;
        y = ystart;
        inboundsx=true;
        inboundsy=true;
    end
    while inboundsy
        while inboundsx
            pos(i,1)=x;
            pos(i,2)=y;
            x = x+dx;
            i=i+1;
            if x>axmax+tolerance
                inboundsx=false;
            end
        end
        x=xstart;
        y=y+dy;
        inboundsx=true;
        if y>aymax+tolerance
            inboundsy = false;
        end
    end
end
    disp(['-> I generated ',num2str(length(pos(:,1))),' atoms.'])
end