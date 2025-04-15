%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 01/29/2025

%this function reformats the data from the read table function 
%accepts table as argument, returns table to getDataVicon

function data = format_Data (data)
    header_Titles = data.Properties.VariableNames;
    fprintf(header_Titles)
    data = 0;
end