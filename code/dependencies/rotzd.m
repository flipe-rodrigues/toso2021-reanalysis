function Rz = rotzd(theta,dim)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    if nargin == 1
        dim = 3;
    end
    
    Rz = [...
        cosd(theta), -sind(theta), 0; ...
        sind(theta),  cosd(theta), 0; ...
                  0,            0, 1; ...
    ];

    Rz = Rz(1:dim,1:dim);
end

