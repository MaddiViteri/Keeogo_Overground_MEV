function nameMap = get_marker_name_map()
    % Returns a map of expected marker names to actual field names in the dataset.
    % Useful for renaming or resolving marker name mismatches.

    nameMap = containers.Map();

    % Core pelvis markers
    nameMap('ULSacral') = 'ULSacral';
    nameMap('LLSacral') = 'LLSacral';
    nameMap('URSacral') = 'URSacral';
    nameMap('LRSacral') = 'LRSacral';
    nameMap('RASIS') = 'RASIS';
    nameMap('LASIS') = 'LASIS';
    nameMap('RPSIS') = 'RPSIS';
    nameMap('LPSIS') = 'LPSIS';

    % Thigh cluster - Left
    nameMap('LAntLatThigh') = 'LAntLatThigh';
    nameMap('LAntMedThigh') = 'LAntMedThigh';
    nameMap('LPosLatThigh') = 'LPosLatThigh';
    nameMap('LPosMedThigh') = 'LPosMedThigh';

    % Thigh bony landmarks - Left
    nameMap('LLatEpi') = 'LLatEpi';
    nameMap('LMedEpi') = 'LMedEpi';

    % Thigh cluster - Right
    nameMap('RAntLatThigh') = 'RAntLatThigh';
    nameMap('RAntMedThigh') = 'RAntMedThigh';
    nameMap('RPosLatThigh') = 'RPosLatThigh';
    nameMap('RPosMedThigh') = 'RPosMedThigh';

    % Thigh bony landmarks - Right
    nameMap('RLatEpi') = 'RLatEpi';
    nameMap('RMedEpi') = 'RMedEpi';

    % Shank cluster - Left
    nameMap('LSupLatShank') = 'LSupLatShank';
    nameMap('LSupMedShank') = 'LSupMedShank';
    nameMap('LInfLatShank') = 'LInfLatShank';
    nameMap('LInfMedShank') = 'LInfMedShank';

    % Ankle - Left
    nameMap('LLatMal') = 'LLatMal';
    nameMap('LMedMal') = 'LMedMal';

    % Shank cluster - Right
    nameMap('RSupLatShank') = 'RSupLatShank';
    nameMap('RSupMedShank') = 'RSupMedShank';
    nameMap('RInfLatShank') = 'RInfLatShank';
    nameMap('RInfMedShank') = 'RInfMedShank';

    % Ankle - Right
    nameMap('RLatMal') = 'RLatMal';
    nameMap('RMedMal') = 'RMedMal';

    % Foot - Left
    nameMap('LHEE') = 'LHEE';
    nameMap('L5Meta') = 'L5Meta';
    nameMap('LTOE') = 'LTOE';
    nameMap('L2Meta') = 'L2Meta';

    % Foot - Right
    nameMap('RHEE') = 'RHEE';
    nameMap('R5Meta') = 'R5Meta';
    nameMap('RTOE') = 'RTOE';
    nameMap('R2Meta') = 'R2Meta';

    % Upper body - Left
    nameMap('LShould') = 'LShould';
    nameMap('LElb') = 'LElb';
    nameMap('LRad') = 'LRad';
    nameMap('LUln') = 'LUln';

    % Upper body - Right
    nameMap('RShould') = 'RShould';
    nameMap('RElb') = 'RElb';
    nameMap('RRad') = 'RRad';
    nameMap('RUln') = 'RUln';

    % Trunk and head
    nameMap('Stern') = 'Stern';
    nameMap('C7') = 'C7';
    nameMap('LFHead') = 'LFHead';
    nameMap('RFHead') = 'RFHead';
    nameMap('LBHead') = 'LBHead';
    nameMap('RBHead') = 'RBHead';

    % Example alias mapping (uncomment/add as needed)
    % nameMap('RAntInfLatThigh') = 'RAntLatThigh';
    % nameMap('RAntInfMedThigh') = 'RAntMedThigh';
end
