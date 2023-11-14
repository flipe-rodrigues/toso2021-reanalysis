function Ry = rotyd(theta,dim)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    if nargin == 1
        dim = 3;
    end

    Ry = [...
         cosd(theta), 0, sind(theta); ...
                   0, 1,           0; ...
        -sind(theta), 0, cosd(theta); ...
    ];

    Ry = Ry(1:dim,1:dim);
end

