%%% This code was adapted from Ahlad Neti's emgAnalysis code from the CATT
%%% study. The code filters the data using both RMS and traditional
%%% filtering (bandpass filter, rectify, low-pass filter).
%%% The code splits the MVC trials and the regular trials
function [all_data,trial_names] = emgAnalysis(emg_path,muscle_names)
% Create filter
fs = 2000;
    fs_nyq = fs/2;
    win_smooth = 0.1;
    win_size = win_smooth*fs;
    bp_low = 10; % prev 10
    bp_high = 500; % prev 500
    fc_low = 2;
    Wn1 = bp_low/fs_nyq;
    Wn2 = bp_high/fs_nyq;
    Wn = fc_low/fs_nyq;
    
    [bb,ba] = butter(4,[Wn1,Wn2],'bandpass');
    [d,c] = butter(4,Wn,'low');
    movRMS = dsp.MovingRMS(win_size,"Method","Sliding window","SpecifyWindowLength",true);

    c1 = 1;
    emg_filenames = {dir(emg_path).name}; % Extract filenames
    for n = 3:length(emg_filenames) % Ignore empty cells
        true_name = emg_filenames{n};
        name = true_name(18:end-4); % name of the trial
        % Load raw data
        raw_data = load(strcat(emg_path,"\",true_name));
        f = fieldnames(raw_data);
        signal_data = raw_data.(f{:}).movements.sources(1).signals;
        sigs = fieldnames(signal_data);

        if ~contains(upper(name), "MVC") % Determine actual trials
            c2 = 1;
            for s = 1:length(sigs)
                sig_name = sigs{s};
                % muscle_names = raw_data.(f{:}).movements.sources(1).signals.(sigs{s}).name;
                signal = signal_data.(sig_name).data;
                data = filtfilt(bb,ba,signal);
                % data_rms = filtfilt(d,c,data);
                recdata = abs(data - mean(data));
                smoothed = filtfilt(d,c,recdata);
                rms = movRMS(data);
                all_data.(name).emg.(muscle_names{c2}) = smoothed;
                all_data.(name).rms.(muscle_names{c2}) = rms; %fix
                all_data.(name).time = 1/fs:1/fs:length(all_data.(name).emg.(muscle_names{1}))/fs;
                c2 = c2 + 1;
            end
            trial_names{c1} = name;
            c1 = c1 + 1;
        else
            c2 = 1;
            for s = 1:length(sigs)
                sig_name = sigs{s};
                signal = signal_data.(sig_name).data;
                data = filtfilt(bb,ba,signal);
                % data_rms = filtfilt(d,c,data);
                recdata = abs(data - mean(data));
                % smoothed = smoothdata(recdata,'sgolay',win_size);
                smoothed = filtfilt(d,c,recdata);
                rms = movRMS(data);
                
                all_data.MVC.(name).(muscle_names{c2}) = smoothed;
                all_data.MVC.(name).rms.(muscle_names{c2}) = rms;
                all_data.MVC.(name).time = 1/fs:1/fs:(length(all_data.MVC.(name).(muscle_names{1})))/fs;

                % all_data.MVC.(name).(muscle_names{c2}) = smoothed;
                c2 = c2 + 1;
            end
        end
    end
end