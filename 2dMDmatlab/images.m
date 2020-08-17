function [ sx, sy ] = images( sx, sy, pbc, xbox, ybox )

if pbc==true
    %if sx>xbox/2 sx=sx-xbox; end
    %if sx<-xbox/2 sx=sx+xbox; end
    if sy>ybox/2 sy=sy-ybox; end
    if sy<-ybox/2 sy=sy+ybox; end
end

end

