function marker = get_marker(markerStruct, oldName, nameMap)
    if isfield(markerStruct, oldName)
        marker = markerStruct.(oldName);
    elseif isfield(nameMap, oldName)
        mappedName = nameMap.(oldName);
        if isfield(markerStruct, mappedName)
            marker = markerStruct.(mappedName);
        else
            warning("Mapped field '%s' also not found. Available fields:\n%s", ...
                mappedName, strjoin(fieldnames(markerStruct), ', '));
            error("Marker field not found: %s", oldName);
        end
    else
        warning("No mapping exists for field: %s", oldName);
        error("Marker field not found: %s", oldName);
    end
end
