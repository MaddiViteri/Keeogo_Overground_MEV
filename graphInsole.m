%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/20/2025 (MATLAB R2024a)

%where foot is left or right

%-------------------------------------------------------------------------%

function graphInsole (foot,trialInsole)
    
foot = caseInsensitivePattern(foot);

    foot = string(foot); % Ensure foot is a string

    if contains(foot, 'l')
        I = imread("INSOLEL.png");
        I = imrotate(I, 180); % Fix orientation for left insole
        desiredRows = trialInsole.forefootWidth;
        desiredColumns = trialInsole.length;
        K = imresize(I, [desiredColumns, desiredRows]); 
        imshow(K)

    elseif contains(foot, 'r')
        I = imread("INSOLER.png");
        desiredRows = trialInsole.forefootWidth;
        desiredColumns = trialInsole.length;
        J = imresize(I, [desiredColumns, desiredRows]);
        imshow(J)
    
    else 
        disp("Error: image cannot be added to graph");
    end 

end%end function