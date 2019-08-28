NUM_PARCELS = 12;

parcel_names=["Right occipital", "Left occipital", "Left temporal", "Right temporal", ...
              "Left frontal", "Right frontal", "Left pre-frontal lobe", ...
              "Right pre-frontal lobe", "Left parietal", "Right parietal", ...
              "Left limbic lobe", "Right limbic lobe"];

mwf_p_values = zeros(NUM_PARCELS, 1);
tf_p_values = zeros(NUM_PARCELS, 1);
for i=1:NUM_PARCELS
    parcel_ep = ep_mwf_per_parcel{i};
    parcel_ft = ft_mwf_per_parcel{i};
    y = vertcat(parcel_ep, parcel_ft);
    g = vertcat(repmat("Preterm", size(parcel_ep, 1), 1), repmat("Term", size(parcel_ft, 1), 1));
    p = anova1(y, g, 'off');
    mwf_p_values(i) = p;
    f1 = figure('Name','Myelin Water Fraction analysis of term and preterm subjects',...
               'NumberTitle','off','Visible','off');    
    boxplot(y, g);
    title(parcel_names(i));
    ylabel('liters');
    saveas(f1,[pwd strcat('/mwf_boxplot_', num2str(i), '.png')]);
    
    parcel_ep = ep_tf_per_parcel{i};
    parcel_ft = ft_tf_per_parcel{i};
    y = vertcat(parcel_ep, parcel_ft);
    g = vertcat(repmat("Preterm", size(parcel_ep, 1), 1), repmat("Term", size(parcel_ft, 1), 1));
    p = anova1(y, g, 'off');
    tf_p_values(i) = p;
    f2 = figure('Name','Tissue Fraction analysis of term and preterm subjects',...
                'NumberTitle','off','Visible','off');     
    boxplot(y, g);
    title(parcel_names(i));
    ylabel('liters');
    saveas(f2,[pwd strcat('/tf_boxplot_', num2str(i), '.png')]);
end

    %[R, p_value] = corrcoef(parcel_ep, parcel_ft);
    %f = figure('Name','Time-intensity curves of individual voxels or regions of interest',...
    %           'NumberTitle','off','Visible','off');
%     scatter(parcel_ep, parcel_ft, '.');
%     h1 = lsline;
%     h1.Color = 'r';
%     h1.LineWidth = 3;
%     title(parcel_names(i));
%     coord_x = max(parcel_ep)*0.75;
%     coord_y = max(parcel_ft)*0.93;
%     text(coord_x, coord_y, strcat("R^2=", num2str(R(2))),'FontSize',12,'FontWeight','bold');
%     text(coord_x, coord_y+0.05, strcat("p-value=", num2str(p_value(2))), 'FontSize',12,'FontWeight','bold');
%     xlabel('Extremly Preterm MWF');
%     ylabel('Full term MWF');
