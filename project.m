% Import data (assuming CSV format)
data = readtable('pwdfaults_dataset.csv');

% %Dataset checking
% %Check for missing values
% missing_values = sum(ismissing(data));
% disp('Missing values per column:');
% disp(missing_values);
% %Basic statistics
% % Get summary statistics for all variables
% summary(data);



% %eida thik ase Chi-squre
% %Chi-Squere test in fault type data VS Weather Condition data and Visualization
% data.WeatherCondition = categorical(data.WeatherCondition);
% data.FaultType = categorical(data.FaultType);
% 
% % Create contingency table and perform chi-square test
% [tbl, chi2, p] = crosstab(data.FaultType, data.WeatherCondition);
% fprintf('Chi-square test results:\nchi^2 = %.2f, p-value = %.4f\n', chi2, p);
% 
% % Create heatmap with title
% h = heatmap(categories(data.WeatherCondition), categories(data.FaultType), tbl');
% h.Title = 'Chi-Square Test: Fault Type vs Weather Condition';
% h.XLabel = 'Weather Condition';
% h.YLabel = 'Fault Type';
% h.ColorbarVisible = 'on';
% 
% % Add statistical results as subtitle
% annotation('textbox', [0.15, 0.85, 0.1, 0.1], 'String', ...
%     sprintf('?² = %.2f, p = %.4f', chi2, p), ...
%     'EdgeColor', 'none', 'FontSize', 10);




% %eida thik ase T-test er jonno
% data.FaultType = strtrim(lower(data.FaultType)); % Normalize text
% 
% % Define the fault types we want to compare
% fault_types = {'line breakage', 'transformer failure', 'overheating'};
% 
% % 2. Check sample counts for each fault type
% sample_counts = cellfun(@(x) sum(strcmp(data.FaultType, x)), fault_types);
% disp('Sample counts per fault type:');
% disp(table(fault_types', sample_counts', 'VariableNames', {'FaultType','Count'}));
% 
% % 3. Compare Voltage and Power Load Demand between groups
% if all(sample_counts >= 2) % Need at least 2 samples per group
%     fprintf('\n=== VOLTAGE COMPARISON (V) ===\n');
%     compare_metrics(data, 'Voltage_V', fault_types);
%     
%     fprintf('\n=== POWER LOAD DEMAND COMPARISON (MW) ===\n');
%     compare_metrics(data, 'PowerLoad_Demand_MW', fault_types);
%     
%     % 4. Create ECDF visualizations instead of boxplots
%     create_ecdf_plots(data, {'Voltage_V', 'PowerLoad_Demand_MW'}, fault_types);
% else
%     fprintf('\nCannot perform t-tests - insufficient samples in one or more groups\n');
%     fprintf('At least 2 samples per group are required\n');
% end
% 
% % ========== HELPER FUNCTIONS ==========
% 
% function compare_metrics(data, metric, groups)
%     % Perform pairwise t-tests between all groups for a given metric
%     for i = 1:length(groups)-1
%         for j = i+1:length(groups)
%             group1 = data.(metric)(strcmp(data.FaultType, groups{i}));
%             group2 = data.(metric)(strcmp(data.FaultType, groups{j}));
%             
%             % Remove missing values
%             group1 = group1(~isnan(group1));
%             group2 = group2(~isnan(group2));
%             
%             fprintf('\n%s vs %s (%s)\n', groups{i}, groups{j}, metric);
%             fprintf('Group sizes: %d vs %d\n', length(group1), length(group2));
%             
%             if length(group1) < 2 || length(group2) < 2
%                 fprintf('Insufficient samples for t-test\n');
%             else
%                 % Perform Welch's t-test (unequal variances)
%                 [h,p,ci,stats] = ttest2(group1, group2, 'Vartype', 'unequal');
%                 
%                 % Display results
%                 fprintf('t(%0.1f) = %.2f, p = %.4f\n', stats.df, stats.tstat, p);
%                 fprintf('Mean ± SD: %.2f ± %.2f vs %.2f ± %.2f\n', ...
%                     mean(group1), std(group1), mean(group2), std(group2));
%                 fprintf('95%% CI for difference: [%.2f, %.2f]\n', ci(1), ci(2));
%             end
%         end
%     end
% end
% 
% function create_ecdf_plots(data, metrics, groups)
%     % Create ECDF plots for each metric
%     colors = lines(length(groups)); % Different colors for each group
%     
%     for m = 1:length(metrics)
%         figure;
%         hold on;
%         legend_entries = cell(1, length(groups));
%         
%         for g = 1:length(groups)
%             vals = data.(metrics{m})(strcmp(data.FaultType, groups{g}));
%             vals = vals(~isnan(vals)); % Remove NaN values
%             
%             if ~isempty(vals)
%                 [f, x] = ecdf(vals);
%                 stairs(x, f, 'Color', colors(g,:), 'LineWidth', 2);
%                 legend_entries{g} = sprintf('%s (n=%d)', groups{g}, length(vals));
%             end
%         end
%         
%         title(['ECDF of ' metrics{m} ' by Fault Type']);
%         xlabel(metrics{m});
%         ylabel('Cumulative Probability');
%         legend(legend_entries(~cellfun(@isempty, legend_entries)), 'Location', 'southeast');
%         grid on;
%         
%         % Add KS test results to plot
%         if length(groups) == 2
%             group1 = data.(metrics{m})(strcmp(data.FaultType, groups{1}));
%             group2 = data.(metrics{m})(strcmp(data.FaultType, groups{2}));
%             [~, p_ks] = kstest2(group1, group2);
%             annotation('textbox', [0.2, 0.7, 0.1, 0.1], 'String', ...
%                 sprintf('KS-test p = %.3f', p_ks), 'EdgeColor', 'none');
%         end
%     end
% end

% 
% 
% % eida thik ase Monte Carlo 
% % Simulate fault probabilities under different weather conditions
% num_simulations = 100000;
% weather_types = categories(categorical(data.WeatherCondition));
% fault_probs = groupsummary(data, 'WeatherCondition', 'mean', 'DownTime_hrs');
% 
% % Get observed maximum downtime
% max_downtime = max(data.DownTime_hrs); % This will be 10 in your case
% 
% % Monte Carlo simulation with truncation
% simulated_downtimes = zeros(num_simulations, 1);
% for i = 1:num_simulations
%     % Randomly select a weather condition
%     rand_weather = weather_types(randi(length(weather_types)));
%     idx = strcmp(fault_probs.WeatherCondition, rand_weather);
%     mu = fault_probs.mean_DownTime_hrs(idx);
%     
%     % Generate truncated exponential values
%     while true
%         x = exprnd(mu);
%         if x <= max_downtime
%             simulated_downtimes(i) = x;
%             break;
%         end
%     end
% end
% 
% % Analyze results
% figure;
% histogram(simulated_downtimes, 'Normalization', 'probability', 'BinWidth', 0.5);
% title('Simulation of Downtime Distribution (Truncated at 10 hours)');
% xlabel('Downtime (hours)');
% ylabel('Probability');
% 
% xlim([0 max_downtime*1.2]); % Show up to 10% beyond max
% 
% % Calculate and display statistics
% fprintf('Simulation Statistics:\n');
% fprintf('Maximum observed downtime: %.2f hours\n', max_downtime);
% fprintf('Maximum simulated downtime: %.2f hours\n', max(simulated_downtimes));
% fprintf('95%% of simulated downtimes < %.2f hours\n', prctile(simulated_downtimes, 95));
% fprintf('Mean simulated downtime: %.2f hours\n', mean(simulated_downtimes));

% 
% 
%eida thik ase ks-test
% Define all weather conditions to compare
weather_conditions = {'rainy', 'clear', 'thunderstorm', 'windstorm'};
colors = {'r', 'b', 'g', 'm'}; % Colors for each condition
line_styles = {'-', '--', ':', '-.'}; % Different line styles

% Initialize variables
downtimes = cell(1, length(weather_conditions));
valid_conditions = {};

% Extract downtime data for each weather condition
for i = 1:length(weather_conditions)
    condition = weather_conditions{i};
    downtimes{i} = data.DownTime_hrs(strcmpi(data.WeatherCondition, condition) & ~isnan(data.DownTime_hrs));
    
    % Only keep conditions with at least 2 data points
    if length(downtimes{i}) >= 2
        valid_conditions{end+1} = condition;
    else
        fprintf('Warning: Insufficient data (%d points) for %s condition\n', length(downtimes{i}), condition);
    end
end

% Perform pairwise KS tests between all valid conditions
fprintf('\n=== Kolmogorov-Smirnov Test Results ===\n');
for i = 1:length(valid_conditions)-1
    for j = i+1:length(valid_conditions)
        [h,p,ksstat] = kstest2(downtimes{strcmp(weather_conditions, valid_conditions{i})}, ...
                               downtimes{strcmp(weather_conditions, valid_conditions{j})});
        fprintf('%s vs %s: D = %.3f, p = %.4f\n', ...
                valid_conditions{i}, valid_conditions{j}, ksstat, p);
    end
end

% Create ECDF plot for all valid conditions
figure;
hold on;
legend_entries = cell(1, length(valid_conditions));

for i = 1:length(valid_conditions)
    cond_idx = find(strcmp(weather_conditions, valid_conditions{i}));
    [f,x] = ecdf(downtimes{cond_idx});
    plot(x, f, 'Color', colors{cond_idx}, 'LineStyle', line_styles{i}, 'LineWidth', 2);
    legend_entries{i} = sprintf('%s (n=%d)', valid_conditions{i}, length(downtimes{cond_idx}));
end

% Format plot
title('ECDF Comparison of Downtime by Weather Condition');
xlabel('Downtime (hours)');
ylabel('Cumulative Probability');
legend(legend_entries, 'Location', 'southeast');
grid on;
