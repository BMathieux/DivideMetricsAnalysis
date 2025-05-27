function ws = SelectWs(HMetrics, interactiveWs, defaultWs)
    % Select window size for hillslope metric calculations.
    % This function plots mean and standard deviation of hilltop curvature to allow 
    % interactive selection of the window size (ws), or returns a default value if 
    % interactive selection is disabled. The selected window size is used for curvature 
    % analysis in divide metrics.
    %
    % Inputs:
    %   HMetrics      - Hillslope metrics from basin analysis.
    %   interactiveWs - Enable interactive plot for selection (logical, true/false).
    %   defaultWs     - Default window size if non-interactive (numeric, meters).
    %
    % Outputs:
    %   ws            - Selected window size (meters).

    % --- Use Default if Non-Interactive ---
    if ~interactiveWs
        ws = defaultWs;
        fprintf('Using default ws: %d m\n', ws);
        return;
    end

    % --- Compute Metrics ---
    [mC, sC, ws_vals] = computeMetrics(HMetrics);

    % --- Select Window Size ---
    ws = selectWindowSize(mC, sC, ws_vals);

    %% --- SUBFUNCTIONS ---

    % --- Compute Metrics ---
    function [mC, sC, ws_vals] = computeMetrics(HMetrics)
        % Calculate mean and standard deviation of curvature for each window size.
        ws_vals = HMetrics{1}.ws;
        mC = nan(numel(HMetrics), numel(ws_vals));
        sC = nan(numel(HMetrics), numel(ws_vals));

        for i = 1:numel(HMetrics)
            for j = 1:numel(ws_vals)
                ci = HMetrics{i}.Cht(:, j);
                ci = ci(~isnan(ci));
                mC(i, j) = mean(ci);
                sC(i, j) = std(ci);
            end
        end
    end

    % --- Select Window Size ---
    function ws = selectWindowSize(mC, sC, ws_vals)
        % Interactively select window size using curvature plots.
        f = figure('Name', 'Select Window Size', 'Position', [100, 100, 900, 600]);

        % Mean curvature plot
        ax1 = subplot(2, 2, 1); hold on; box on;
        plot(ax1, ws_vals, mC', 'Color', [0.5 0.5 0.5]);
        plot(ax1, ws_vals, mean(mC, 1, 'omitnan'), 'k', 'LineWidth', 3);
        xlabel('Window length (m)'); ylabel('Mean curvature (1/m)');
        grid on; grid minor;

        % Standard deviation plot
        ax2 = subplot(2, 2, 2); hold on; box on;
        plot(ax2, ws_vals, sC', 'Color', [0.5 0.5 0.5]);
        plot(ax2, ws_vals, mean(sC, 1, 'omitnan'), 'k', 'LineWidth', 3);
        xlabel('Window length (m)'); ylabel('Std curvature (1/m)');
        grid on; grid minor;

        % Gradient plot
        ax3 = subplot(2, 2, [3 4]); hold on; box on;
        gradC = abs(gradient(mean(mC, 1, 'omitnan')));
        plot(ax3, ws_vals, gradC, 'k', 'LineWidth', 3);
        xlabel('Window length (m)'); ylabel('gradient (mean curvature)');
        set(ax3, 'YScale', 'log');
        grid on; grid minor;

        % Initialize slider and vertical lines
        xL1 = xline(ax1, ws_vals(1), 'r-', 'LineWidth', 2);
        xL2 = xline(ax2, ws_vals(1), 'r-', 'LineWidth', 2);
        xL3 = xline(ax3, ws_vals(1), 'r-', 'LineWidth', 2);

        sld = uicontrol('Style', 'slider', 'Min', 1, 'Max', numel(ws_vals), ...
            'Value', 1, 'Units', 'normalized', 'Position', [0.15, 0.02, 0.6, 0.04], ...
            'SliderStep', [1/(numel(ws_vals)-1) 1/(numel(ws_vals)-1)]);

        txt = uicontrol('Style', 'text', 'Units', 'normalized', ...
            'Position', [0.76, 0.02, 0.1, 0.04], 'String', num2str(ws_vals(1)), ...
            'FontSize', 14, 'FontWeight', 'bold');

        uicontrol('Style', 'pushbutton', 'String', 'Confirm', ...
            'Units', 'normalized', 'Position', [0.86, 0.02, 0.1, 0.04], ...
            'Callback', @(~,~) uiresume(f));

        addlistener(sld, 'ContinuousValueChange', @(src, ~) updateSlider(round(src.Value)));

        updateSlider(1);
        uiwait(f);
        ws = ws_vals(round(sld.Value));
        fprintf('Selected ws: %d m\n', ws);
        close(f);

        % Update slider callback
        function updateSlider(idx)
            xL1.Value = ws_vals(idx);
            xL2.Value = ws_vals(idx);
            xL3.Value = ws_vals(idx);

            delete(findall(ax1, 'Tag', 'wsText'));
            delete(findall(ax2, 'Tag', 'wsText'));
            delete(findall(ax3, 'Tag', 'wsText'));

            text(ax1, ws_vals(idx), ax1.YLim(2), sprintf('%d m', ws_vals(idx)), ...
                'Color', 'k', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', 'Tag', 'wsText');
            text(ax2, ws_vals(idx), ax2.YLim(2), sprintf('%d m', ws_vals(idx)), ...
                'Color', 'k', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', 'Tag', 'wsText');
            text(ax3, ws_vals(idx), ax3.YLim(2), sprintf('%d m', ws_vals(idx)), ...
                'Color', 'k', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', 'Tag', 'wsText');

            txt.String = num2str(ws_vals(idx));
            drawnow;
        end
    end
end