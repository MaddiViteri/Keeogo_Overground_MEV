%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/21/2025 (MATLAB R2024a)

%this program converts string of height in feet to meters

%-------------------------------------------------------------------------%

function [height] = fToM (participant)
    
    height = 0;

    hfeet_String = participant.height;

    height_split = split(hfeet_String,"'");

    height_inches = (str2num(height_split(1))*12)+str2num(height_split(2));

    height = height_inches/39.37; %convert to meters

end%end function