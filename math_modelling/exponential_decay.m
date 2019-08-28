NUM_CASES = 139;
NUM_SEGMENTS = 4;

signal_per_segment = cell(NUM_SEGMENTS, 1);
for case_number=2:NUM_CASES
    echo_time_values = TEs;
    n_echotimes = size(echo_time_values, 2);
    for segment = 1:NUM_SEGMENTS
        if isempty(signal_per_segment{segment})
            segment_map = containers.Map;
        else
            segment_map = signal_per_segment{segment};
        end
        segmentation_masks = all_data{case_number, 2};
        segment_mask = segmentation_masks(:, segment);
        segment_indexes = find(segment_mask);
        case_segment_voxels = all_data{case_number, 1};
        case_signals = case_segment_voxels(segment_indexes, :).*segment_mask(segment_indexes);
        signal_in_segment = mean(case_signals);
        for echo_idx=1:n_echotimes
            echo_time = int2str(echo_time_values(echo_idx));
            signal_value = signal_in_segment(echo_idx);
            if segment_map.isKey(echo_time)
                t_value = segment_map(echo_time);
                segment_map(echo_time) = (t_value+signal_value)/2;
            else
                segment_map(echo_time) = signal_value;
            end
        end
        signal_per_segment{segment} = segment_map;
    end
end

echotimes_per_segment = cell(NUM_SEGMENTS, 1);
for segment = 1:NUM_SEGMENTS
     segment_map = signal_per_segment{segment};
     echo_times = keys(segment_map);
     echo_times_arr = [];
     signal_segment_arr = [];
     for echo_idx=1:size(echo_times, 2)
         et = echo_times{echo_idx};
         intensity = segment_map(et);
         echo_times_arr = [echo_times_arr; str2num(et)];
         signal_segment_arr = [signal_segment_arr; intensity];
     end
     echotimes_per_segment{segment} = sort(echo_times_arr);
     signal_per_segment{segment} = sort(signal_segment_arr, 'descend');
end

f = figure('Name','Time-intensity curves of individual voxels or regions of interest',...
           'NumberTitle','off','Visible','off');
plot(echotimes_per_segment{2}, signal_per_segment{2}, '-o', ...
     echotimes_per_segment{3}, signal_per_segment{3}, '-x', ...
     echotimes_per_segment{4}, signal_per_segment{4}, '-+');
xlabel('Echo Time');
ylabel('Intensity');
legend('CSF', 'Grey Matter', 'White Matter');
saveas(f,[pwd '/figures_q1/signal_segment_big_data.png']);