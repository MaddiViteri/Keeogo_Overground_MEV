%% Example code to parse data recorded in the Keeogo Data Logging Format (KDLF).
%% THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.
% This function takes the name of a file containing KDLF data as input parameter,
% performs basic validation of the file, parses the data and returns a data
% structure containing the parsed data in physical units. Fields in the
% data structure are named as in TEC-00446, e.g. Data.LK_Angle contains the
% angle data of the left knee.
% Version 1.0
% Compatible with Keeogo firmware version 5.14.0.0
function Data = ParseKDLF(FileName)

    %Open file with read-only permission.
    f = fopen(FileName, 'r');
    
    %% Parse file header
    %Declarations of expected values for fields in the header.
    SIGNATURE = hex2dec("0XA245");
    VERSION = 1;
    HEADER_SIZE = 248;
    NB_DATA_FIELDS = 47;
    FRAME_SIZE = 98;

    signature = fread(f, 1, '*uint16');
    
    %Exit prematurely if passed file has the wrong signature.
    if signature ~= SIGNATURE
        fclose(f);
        error("Header has wrong signature");
    end

    version = fread(f, 1, '*uint8');

    %Exit prematurely if passed file has the wrong version.
    if version ~= VERSION
        fclose(f);
        error("Header has wrong version");
    end

    headerSize = fread(f, 1, '*uint16');

    %Exit prematurely if passed file has the wrong header size.
    if headerSize ~= HEADER_SIZE
        fclose(f);
        error("Header has wrong size");
    end

    numberOfDataFields = fread(f, 1, '*uint16');

    %Exit prematurely if passed file has the wrong number of data fields.
    if numberOfDataFields ~= NB_DATA_FIELDS
        fclose(f);
        error("Header has wrong number of data fields");
    end

    frameSize = fread(f, 1, '*uint16', 2);

    %Exit prematurely if passed file has the wrong frame size.
    if frameSize ~= FRAME_SIZE
        fclose(f);
        error("Header has wrong frame size");
    end

    %% Parse data frames
    %Declarations of conversion factors to obtain physical values from 
    %fixed-point values.
    ANGLE_GAIN_DEGREES = 180 / 32768;       %Gain to obtain angle in degrees.
    ACCEL_GAIN_MG = 4000 / 32768;           %Gain to obtain acceleration in milli-g.
    VEL_GAIN_DPS_KNEES = 573.44 / 32768;    %Gain to obtain angular velocity in degrees per second;
    VEL_GAIN_DPS_HIPS = 2000 / 32768;       %Gain to obtain angular velocity in degrees per second;
    EULER_GAIN_DEGREES = 1 / 128;           %Gain to obtain Euler angles in degrees.
    TORQUE_GAIN_NM = 1 / 152;               %Gain to obtain torque in Nm.
    TEMP_GAIN_CELSIUS = 1 / 256;            %Gain to obtain temperature in Celsius.

    %Set offset to the first data item of the first data frame.
    offset = HEADER_SIZE + 2 + 2;
    fseek(f, offset, 'bof');

    SKIP = FRAME_SIZE - 2;

    Data.RK_Angle = fread(f, inf, 'int16=>double', SKIP) .* ANGLE_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_Angle = fread(f, inf, 'int16=>double', SKIP) .* ANGLE_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_AccFrontal = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_AccVertical = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_AccLateral = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_AccFrontal = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_AccVertical = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_AccLateral = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_GyroFrontal = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_KNEES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_GyroVertical = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_KNEES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_GyroLateral = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_KNEES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_GyroFrontal = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_KNEES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_GyroVertical = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_KNEES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_GyroLateral = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_KNEES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_Roll = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_Pitch = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_Yaw = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_Roll = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_Pitch = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_Yaw = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_Torque = fread(f, inf, 'int16=>double', SKIP) .* TORQUE_GAIN_NM;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_Torque = fread(f, inf, 'int16=>double', SKIP) .* TORQUE_GAIN_NM;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RK_Temperature = fread(f, inf, 'int16=>double', SKIP) .* TEMP_GAIN_CELSIUS;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LK_Temperature = fread(f, inf, 'int16=>double', SKIP) .* TEMP_GAIN_CELSIUS;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_Angle = fread(f, inf, 'int16=>double', SKIP) .* ANGLE_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_Angle = fread(f, inf, 'int16=>double', SKIP) .* ANGLE_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_AccFrontal = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_AccVertical = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_AccLateral = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_AccFrontal = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_AccVertical = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_AccLateral = fread(f, inf, 'int16=>double', SKIP) .* ACCEL_GAIN_MG;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_GyroFrontal = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_HIPS;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_GyroVertical = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_HIPS;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_GyroLateral = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_HIPS;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_GyroFrontal = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_KNEES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_GyroVertical = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_KNEES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_GyroLateral = fread(f, inf, 'int16=>double', SKIP) .* VEL_GAIN_DPS_KNEES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_Roll = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_Pitch = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.RH_Yaw = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_Roll = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_Pitch = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.LH_Yaw = fread(f, inf, 'int16=>double', SKIP) .* EULER_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.ModeSwitchAngle = fread(f, inf, 'int16=>double', SKIP) .* ANGLE_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.AssistanceSwitchAngle = fread(f, inf, 'int16=>double', SKIP) .* ANGLE_GAIN_DEGREES;

    offset = offset + 2;
    fseek(f, offset, 'bof');

    Data.OnOffSwitchAngle = fread(f, inf, 'int16=>double', SKIP) .* ANGLE_GAIN_DEGREES;

    fclose(f);

end