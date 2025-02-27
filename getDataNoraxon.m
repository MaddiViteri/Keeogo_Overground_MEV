%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/23/2025

%This program will import data collected from Noraxon and convert it to a
%matrix readable by matlab

%This function accepts a string as path, retrives data from that folder and
% outputs the data as matrirx data_Vicon

%assumes that each subject folder has a folder labled Noraxon holding the
%motion capture data

%also assumes that the data was exported as separate .csv
%i don't know why i chose this but that's what we're working with

%most og this is similar to the getDataVicon function, just more involved
%to concatenate files into tables

%Edit History
%02/23/2025 - changed tables/cells to structures

%-------------------------------------------------------------------------%

function [data_Noraxon] = getDataNoraxon (path)
    
    path_Noraxon = path + "/Noraxon";

    fprintf(path); %troubleshooting
    fprintf ('\n'); %troubleshooting

    %matrix to rename files
    trialnames_TRUE = ["NONKeeogo_SS01","NONKeeogo_SS02","NONKeeogo_FS01","NONKeeogo_FS02","_Keeogo_SS01",...
        "_Keeogo_SS02","_Keeogo_FS01","_Keeogo_FS02","MVC_LGas","MVC_RGas","MVC_LHam","MVC_RHam","MVC_LQuad","MVC_RQuad"];

    %initialize empty data array
    numStructs = length(trialnames_TRUE);
    data_Noraxon = repmat(struct('trial', [], 'InsoleL', [], 'LegL', [], 'InsoleR', [], 'LegR', []), numStructs, 1);

    %start for loop to go through all of trialnames_TRUE
    %look for FOLDER that matches 

    %for k = 1:length(trialnames_TRUE)
    %choosing a folder
    %path = uigetdir(path); %chooses folder from /Noraxon folder
    %^replaced with automatic compare/contains. these foldernames don't
    %change so no need for user control/decision making

    %compare/find file names from trialnames_true to folder contents
    %test = trialnames_TRUE(1);

    %folderContents = dir(path);
    
    
    fprintf(path_Noraxon + '\n'); %troubleshooting (should print once per folder if this works)

    path_Noraxon = path_Noraxon + "/"; %otherwise it wont read it (very annoying)

    folderContents = dir(path_Noraxon); 

    trialnames = {folderContents.name};

    trialnames = string(trialnames);

    loc = length(trialnames);

    
    found = false;
        %compare fileNames to trialnames_True and run the through for loop
        for j = 1:length(trialnames_TRUE)
            
        
            pat = caseInsensitivePattern(trialnames_TRUE(j)); %accounts for noncase sensitive patterns (an array)

            for k = 1:length(trialnames)
                if (contains(trialnames(k),pat,'IgnoreCase',false)) %also accounts for noncase sensitive patterns (you need both)
                    %disp(trialnames(i) + " file matched to trialnames_TRUE: " + trialnames_TRUE(j)); %troubleshooting
                    path = strcat(path_Noraxon,trialnames(k));
                    [data_LeftInsole, data_LeftLeg, data_RightInsole, data_RightLeg] = DataNoraxonConcat (path);
                    temp = struct('trial',trialnames_TRUE(j),'InsoleL',data_LeftInsole,'LegL',data_LeftLeg,'InsoleR',data_RightInsole, 'LegR', data_RightLeg);
                    data_Noraxon(j) = temp;
                    disp(trialnames_TRUE(j) + " added to data_Noraxon")
                    loc = i; found = 1;
                end%end of if

            end %trialnames for loop
    
            %catch all -- if loc = length of trialnames array, fileisnot found b/c last file is info.csv or something
            if (found == 0) 
                    disp("%-------------------------------------------------------------------------%");
                    disp("!file not found! no file matched to trialnames_TRUE: " + trialnames_TRUE(j));
                    disp("%-------------------------------------------------------------------------%");
            end %end if catch all

            found = false;

        end%trialnames_TRUE for loop
       
    
end%end of function