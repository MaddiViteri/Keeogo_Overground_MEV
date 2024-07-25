% Function is adapted from ktk function detect cycles
% https://kineticstoolkit.uqam.ca/doc/api/ktk.cycles.detect_cycles.html
% [docs]def detect_cycles(
%     ts: TimeSeries,
%     data_key: str,
%     *,
%     event_names: list[str] = ["phase1", "phase2"],
%     thresholds: list[float] = [0.0, 1.0],
%     directions: list[str] = ["rising", "falling"],
%     min_durations: list[float] = [0.0, 0.0],
%     max_durations: list[float] = [np.Inf, np.Inf],
%     min_peak_heights: list[float] = [-np.Inf, -np.Inf],
%     max_peak_heights: list[float] = [np.Inf, np.Inf],
% ) -> TimeSeries:
    % """
    % Detect cycles in a TimeSeries based on a dual threshold approach.
    % 
    % This function detects biphasic cycles and identifies the transitions as
    % new events in the output TimeSeries. These new events are named:
    % 
    % - event_names[0]:
    %   corresponds to the start of phase 1
    % - event_names[1]:
    %   corresponds to the start of phase 2
    % - '_':
    %   corresponds to the end of the cycle.
    % 
    % Parameters
    % ----------
    % ts
    %     TimeSeries to analyze.
    % data_key
    %     Name of the data key to analyze in the TimeSeries. This data must be
    %     unidimensional.
    % event_names
    %     Optional. Event names to add in the output TimeSeries. Default is
    %     ['phase1', 'phase2'].
    % thresholds
    %     Optional. Values to cross to register phase changes. Default is
    %     [0., 1.].
    % directions
    %     Optional. Directions to cross thresholds to register phase changes.
    %     Either ['rising', 'falling'] or ['falling', 'rising']. Default is
    %     ['rising', 'falling'].
    % min_durations
    %     Optional. Minimal phase durations in seconds. Default is [0., 0.].
    % max_durations
    %     Optional. Maximal phase durations in seconds. Default is
    %     [np.Inf, np.Inf]
    % min_peak_heights
    %     Optional. Minimal peak values to be reached in both phases. Default is
    %     [-np.Inf, -np.Inf].
    % max_peak_heights
    %     Optional. Maximal peak values to be reached in both phases. Default is
    %     [np.Inf, np.Inf].
    % 
    % Returns
    % -------
    % TimeSeries
    %     A copy of `ts` with the events added.
    % 
    % """
    % check_types(detect_cycles, locals())
    % 
function [tsout1, heel_strike_data,toe_off_data] = detect_cycles(data,time,thresholds)

    %% Planning
    directions = [{'rising'}, {'falling'}];
    event_names = {'Heel Contact', 'Toe Off'};
    max_durations = [2,4]; % 2s max time between heel strike and toe off and 4s max time between heel strike to heelstrike
    min_peak_heights = -Inf;
    max_peak_heights = Inf;
    directions(1) = lower(directions(1));
    if directions{1} ~= "rising" && directions{1} ~="falling"
        error("""directions(1) must be 'rising' or 'falling'""")
    end
    min_durations =[0.2,0.1];

    % % Testing Code
    % time = trial.downhill_keeogo02.force.time;
    % thresholds = 50;
    % ForceLT = trial.downhill_keeogo02.force.f1;
    % data =ForceLT(:,2);

    events = []; % Initialize events variable

    is_phase1 = true; % Set as phase is true

    for i = 1:length(time)
        % Determine whether data is on or off
        if directions(1) == "rising"
            crossing1 = data(i) >= thresholds;
            crossing2 = data(i) <= thresholds;
        else
            crossing1 = data(i) <= thresholds;
            crossing2 = data(i) >= thresholds;
        end

        if is_phase1 && crossing1
            is_phase1 = false; % Reset to false for next step
            % Record events for heel strike
            events = cat(1,events,[{time(i)}, {event_names{1}}]); 

        elseif (~is_phase1) && crossing2
            is_phase1 = true; % Reset to true for next step
            % Record events for toe off
            events = cat(1,events,[{time(i)}, {event_names{2}}]);
            % events(end+1,:) = [{time(i)},{event_names{2}}];
        end
    end

    % Ensure that we start with event_name1 and that it's not on time0
    while ~contains(string(events{1,2}), event_names{1}) || (events{1,1} == time(1))
        events = events(2:end,:);
    end
%% Remove cycles where criteria are not reached.
    valid_events = [];

    for i_event= 1:2:(length(events)-1)
        time1 = events(i_event,1);
        index_time1 = time1{1}*1000;
        time2 = events(i_event + 1,1);
        index_time2 = time2{1}*1000;
        try
            time3 = events(i_event + 2,1);
            index_time3 = time3{1}*1000;
        catch
            time3 = {Inf};
        end

        sub_ts1 = data(round(index_time1):round(index_time2));
        sub_ts2 = data(round(index_time1):round(index_time3));

        if string(directions(1)) == "rising"
            the_peak1 = max(sub_ts1);
            the_peak2 = min(sub_ts2);
        else
            the_peak1 = min(sub_ts1);
            the_peak2 = max(sub_ts2);
        end

        if (...
            time1{1}>20 ... % Only record after 20 seconds
            && time2{1}< 60 ... % Record only before 60 seconds
            && time2{1} - time1{1} >= min_durations(1)... % Heel strike to toe off time is > min duration
            && time2{1} - time1{1} <= max_durations(1)... % Heel strike to toe off time is < max duration
            && time3{1} - time2{1} <= max_durations(2)... % Heel strike to next heel strike time is > max duration
            && time3{1} - time2{1} >= min_durations(2)... % Toe off to next heel strike time is > max duration
            && the_peak1 >= min_peak_heights...
            && the_peak1 <= max_peak_heights...
            && the_peak2 >= min_peak_heights...
            && the_peak2 <= max_peak_heights...
        )    
            % Save it.
            valid_events= cat(1,valid_events,events(i_event,:));
            valid_events = cat(1,valid_events,events(i_event + 1,:));
        end
        
    end
    %%
    valid_events_str = string(valid_events);
    tsout = unique(valid_events_str,'rows','stable');
    tsout = cellstr(tsout);
    time_steps = str2double(tsout(:,1));
    tsout1 = [num2cell(time_steps),tsout(:,2)];
    %% Plot
    % figure;
    % plot(time,data);
    % xline(cell2mat(tsout1(:,1)),'-.',tsout1(:,2));
    % xlabel('Time (s)'); ylabel('Force (N)');
    %%
    numericvalues = [tsout1{:,1}];
    heel_index = strcmp(tsout(:,2),event_names{1});
    toe_off_index = strcmp(tsout1(:,2),event_names{2});
    heel_strike_data = [num2cell(numericvalues(heel_index))',tsout1(heel_index,2)];
    toe_off_data = [num2cell(numericvalues(toe_off_index))',tsout1(toe_off_index,2)];