%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/21/2025 (MATLAB R2024a)

%this program calculates demographic data from age, height, weight 

%-------------------------------------------------------------------------%

function [age, height, weight] = calcDemographics(participants)

    age_Mean = mean([participants.age]);
    age_SD = std([participants.age]);

    age = struct('mean',age_Mean,'std',age_SD);

    for i = 1:length([participants.height])
        height_Meters = fToM(participants(i));
        participants(i).height = height_Meters;
    end %convert height to meters

    height_Mean = mean([participants.height]);
    height_SD = std([participants.height]);

    height = struct('mean',height_Mean,'std',height_SD);

    for i = 1:length([participants.weight])
        weight_Kg = participants(i).weight / 2.205;
        participants(i).weight = weight_Kg;
    end %convert weight to kg
    
    weight_Mean = mean([participants.weight]);
    weight_SD = std([participants.weight]);

    weight = struct('mean',weight_Mean,'std',weight_SD);

end%end function