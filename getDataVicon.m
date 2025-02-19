%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/07/2025

%This function will import data collected from vicon and convert it to a
%matrix readable by matlab

%This function accepts a string as path, retrives data from that folder and
% outputs the data as matrirx data_Vicon

%assumes that each subject folder has a folder labled Vicon holding the
%motion capture data

%also assumes that the data was exported as .xlsx files

%-------------------------------------------------------------------------%

function [data_Vicon] = getDataVicon (path)

    data_Vicon = {}; %troubleshooting
    
    path = path + "/Vicon"; %Changes path to correct folder

    %fprintf(path); %troubleshooting
    
    %matrix to rename files
    trialnames_TRUE = ["NON_SS01","NON_SS02","NON_FS01","NON_FS02","KON_SS01","KON_SS02","KON_FS01","KON_FS02"];

    %do this 8 times

    for j = 1:length(trialnames_TRUE)

    %choose file
    notice = msgbox("Choose file for trial: " + trialnames_TRUE(j),"Choose file"); %dialog box - notifies which trial the user is selecting the file for

    [f,p]=uigetfile('*.*','File Selector',path); %chooses file from /vicon folder. I'm allowing the user to choose cause the trial names are not always the same

    if isvalid(notice); delete(notice); end %automatically deletes message box

    data_temp = readtable(strcat(p,f));
    %reformat data

    markers = ["SacralUL","SacralMid","SacralUR", "LASIS","RASIS","LPSIS","RPSIS",...
        "LAntSupThigh","LAntInfLatThigh","LAntInfMedThigh","LPosThigh","LLatEpi","LMedEpi",...
        "RAntSupThigh","RAntInfLatThigh","RAntInfMedThigh","RPosThigh","RLatEpi","RMedEpi",...
        "LSupLatShank","LSupMedShank","LInfLatShank","LInfMedShank","LLatMal","LMedMal",...
        "RSupLatShank","RSupMedShank","RInfLatShank","RInfMedShank","RLatMal","RMedMal",...
        "LHee","L5Meta","LToe","L2Meta","RHee","R5Meta","RToe","R2Meta",...
        "Clav","Stern","C7","LShould","RShould"];

    h = convertCharsToStrings(data_temp.Properties.VariableNames); %for the contain function to work, .VariableNames is in char by default

    %parse through markers and h arrays and find column where marker is in

    disp("Reformating " + trialnames_TRUE(j)+ " data from " + f + " ... ");

    loc = 0; 

        %data_Vicon table
        for i=1:length(markers)
            for k=1:length(h)
                if(contains(h(k),markers(i))) 
                 loc = k;
                end
            end
            %disp(markers(i) + " is at index "+ loc); %troubleshootig
            %rename columns 
            data_temp = renamevars(data_temp, [loc loc+1 loc+2], [markers(i)+"_X" markers(i)+"_Y" markers(i)+"_Z"]);
            %disp("columns relabled"); %troubleshooting

            %loc = 0; %resets location
        end
    
    %now no matter how many markers are or what theyre called, they can all
    %be referenced the same in the code

    %put reformated data into correct matrix

    data_Vicon {j} = data_temp;
    
    end %should occur 8 times (one per trial)
    
end %end of function


