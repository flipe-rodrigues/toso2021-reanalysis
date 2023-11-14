function Rx = rotxd(theta,dim)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    if nargin == 1
        dim = 3;
    end
    
    Rx = [...
        1,           0,            0; ...
        0, cosd(theta), -sind(theta); ...
        0, sind(theta),  cosd(theta); ...
    ];

    Rx = Rx(1:dim,1:dim);
end

