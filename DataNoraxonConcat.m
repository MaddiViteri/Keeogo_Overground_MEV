%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/07/2025

%this program concatonates all the noraxon data and outputs it as a table
%back to getDataNoraxon (it was just getting long and hard to troubleshoot
%so i split it up)

%-------------------------------------------------------------------------%

function [data_Small] = DataNoraxonConcat (path)

    %Files needed from each trial folder (taken from noraxon output):
    %LT_MED, LT_RECT, LT_SEMI, RT_MED, RT_RECT, RT_SEMI, 
    %LT_Insole-Arch, LT_Insole-Center_of_pressure, LT_Insole-Heel,
    %LT_Insole-Met, LT_Insole-Total
    %RT_Insole-Arch, RT_Insole-Center_of_pressure, RT_Insole-Heel,
    %RT_Insole-Met, RT_Insole-Total
    fileNames_TRUE = ["LT_MED","LT_RECT","LT_SEMI", "RT_MED", "RT_RECT", "RT_SEMI",...
        "LT_Insole-Arch", "LT_Insole-Center_of_pressure","LT_Insole-Heel", "LT_Insole-Met", "LT_Insole-Total"...
        "RT_Insole-Arch", "RT_Insole-Center_of_pressure","RT_Insole-Heel", "RT_Insole-Met", "RT_Insole-Total"];

    %these need to be found and concatenated into a readable table with
    %variable names

    path = path + "/"; %otherwise it wont read it (very annoying)

    % disp(path); %troubleshooting

    folderContents = dir(path);

    fileNames = {folderContents.name};

    fileNames = string(fileNames);

    %same thing, we're going to look for the files through the h table

    loc = 0;

    data_LeftLeg = {}; data_RightLeg = {}; data_LeftInsole = {}; data_RightInsole = {};

    for j = 1:length(fileNames_TRUE)

        for i = 1:length(fileNames)
            if (contains(fileNames(i),fileNames_TRUE(j)))
                loc = i;
            end % end if find file name
        end %end after length of fileNames

        %catch all
        if (loc == 0)

            disp("no file matched to fileNames_TRUE: " + fileNames_TRUE(j));

        else

            %disp("file matched to fileNames_TRUE: " + fileNames_TRUE(j)); %troubleshooting

            %add data to data_Noraxon - different 
    
            data_temp = readtable(strcat(path,fileNames(loc)));
    
            %rename header variables

            %for fileNames_TRUE(1) - add both columns to data_LeftLeg, rename 'value' to
            %FileNames_True(1)
            %for fileNames_TRUE(2,3) - add value column to data_LeftLeg, rename 'value' to
            %FileNames_True(2,3)
    
            if(j==1)
                data_LeftLeg = data_temp;
                data_LeftLeg = renamevars(data_LeftLeg,'value',fileNames_TRUE(j));
            end
         
    
            %for fileNames_TRUE(2,3) - add value column to data_LeftLeg, rename 'value' to
            %FileNames_True(2,3)
            if(j==2 || j ==3)
                %data_temp = table(data_temp.value);
                data_LeftLeg = join(data_LeftLeg,data_temp); 
                data_LeftLeg = renamevars(data_LeftLeg,'value',fileNames_TRUE(j));
            end
    
            %for fileNames_TRUE(4) - add both columns to data_RightLeg, rename 'value' to
            %FileNames_True(4)
            
            if(j==4)
                data_RightLeg = data_temp;
                data_RightLeg = renamevars(data_RightLeg,'value',fileNames_TRUE(j));
            end
    
            %for fileNames_TRUE(5,6) - add value column to data_RightLeg, rename 'value' to
            %FileNames_True(5,6)

            if(j==5 || j ==6)
                data_RightLeg = join(data_RightLeg,data_temp); 
                data_RightLeg = renamevars(data_RightLeg,'value',fileNames_TRUE(j));
            end
    
            %for fileNames_TRUE(7) - add both columns to data_LeftInsole,
            %rename value to fileNames_TRUE(7)
    
            if(j == 7)
                data_LeftInsole = data_temp;
                data_LeftInsole = renamevars(data_LeftInsole,'value',fileNames_TRUE(j));
            end
    
            %for fileNames_TRUE(8) - add both colums to data_LeftInsole, rename
            %x to COP_X and y to COP_Y
    
            if(j == 8)
                %data_temp = table(data_temp.x,data_temp.y);
                data_LeftInsole = join(data_LeftInsole,data_temp);
                data_LeftInsole = renamevars(data_LeftInsole,{'x', 'y'},["COP_X", "COP_Y"]);
            end
            
            %for fileNames_TRUE(9, 10, 11) - add value column to
            %data_LeftInsole, rename value to fileNames_TRUE(9, 10, 11)
    
            if(j == 9 || j == 10 || j == 11)
                data_LeftInsole = join(data_LeftInsole,data_temp); 
                data_LeftInsole = renamevars(data_LeftInsole,'value',fileNames_TRUE(j));
            end

            %for fileNames_TRUE(12) - add both columns to data_RightInsole,
            %rename value to fileNames_TRUE(12)
    
            if(j == 12)
                data_RightInsole = data_temp;
                data_RightInsole = renamevars(data_RightInsole,'value',fileNames_TRUE(j));
            end
    
            %for fileNames_TRUE(13) - add both colums to data_RightInsole, rename
            %x to COP_X and y to COP_Y
    
            if(j == 13)
                data_RightInsole = join(data_RightInsole,data_temp);
                data_RightInsole = renamevars(data_RightInsole,{'x', 'y'},["COP_X" "COP_Y"]);
            end
    
            %for fileNames_TRUE(14, 15, 16) - add value column to
            %data_RightInsole, rename value to fileNames_TRUE(14, 15, 16)
    
            if(j == 14 || j == 15 || j == 16)
                %data_temp = table(data_temp.Var1);
                data_RightInsole = join(data_RightInsole,data_temp); 
                data_RightInsole = renamevars(data_RightInsole,'value',fileNames_TRUE(j));
            end
        
            data_NoraxonSmall = {data_LeftLeg data_RightLeg data_LeftInsole data_RightInsole}; 
            %disp("data concat sucessful"); %troubleshooting

            data_Small = data_NoraxonSmall;
            %disp("data exported to getDataNoraxon");%troubleshooting

        end%elseif

    end %end after all fileNames_True have been found and added to the table
    
end%funtion
