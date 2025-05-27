function Basins = extractBasins(DEM1, params, varargin)
    % Select divide and extract neighboring basins.
    % This function processes a DEM to delineate basins on either side of a divide
    % and extracts basin characteristics. Outputs include basin data and a resampled DEM.
    % Note: Amin sets the reference drainage area for basin extraction.
    %
    % Inputs:
    %   DEM1     - High-resolution digital elevation model (GRIDobj).
    %   params   - Analysis parameters defined by the DefineParams function.
    %   varargin - Optional settings: 'Parallel' (false), 'Verbose' (true), 
    %              'InteractiveAmin' (false).
    %
    % Outputs:
    %   Basins    - Structure containing basin characteristics and divide data.
    %   DEM2     - Resampled DEM used for basin delineation.

    % --- Parse Inputs ---
    p = inputParser;
    addParameter(p, 'Verbose', true, @islogical);
    addParameter(p, 'InteractiveAmin', false, @islogical);
    parse(p, varargin{:});
    verbose = p.Results.Verbose;
    interactiveAmin = p.Results.InteractiveAmin;

    % --- Resample DEM ---
    DEM2 = resample(DEM1, params.res); % Resample for divide and basin extraction

    % --- Extract River Network ---
    FD2 = FLOWobj(DEM2, 'preprocess', 'carve');
    A2 = flowacc(FD2);
    S2 = STREAMobj(FD2, A2.*DEM2.cellsize^2 > params.A0);

    % --- Define Main Drainage Divide ---
    fprintf('Divide identification...\n');
    D2 = DIVIDEobj(FD2, S2, 'verbose', false); 
    D2 = shrink(D2, FD2, 'distance', params.Dmin);
    [xd, yd] = ind2coord(D2, vertcat(D2.IX));

    ixD = SelectDivide(D2, DEM2, 'verbose', verbose, 'mode', 'single');
    ixD = ixD(:,1); ixD(isnan(ixD)) = [];
    [xdd, ydd] = ind2coord(D2, ixD);
    ixdd = coord2ind(DEM2, xdd, ydd);

    [idx, ~] = rangesearch([xdd ydd], [xd yd], DEM2.cellsize * 2);
    rm = any(cell2mat(cellfun(@(x) ~isempty(x), idx, 'UniformOutput', false)), 2);
    xd(rm) = []; yd(rm) = []; % Delete hilltop nodes at divide

    % --- Set Drainage Area Threshold ---
    params = setThresholds(params, getnal(S2, A2), params.Amin, interactiveAmin, params.nb_ite, DEM2, FD2, A2, xdd, ydd);
    Si = STREAMobj(FD2, (params.A0/DEM2.cellsize^2 < A2) & (A2 < params.Amin));

    % --- Extract Basins ---
    Basins = extractBasins(FD2, Si, xdd, ydd, ixdd, DEM2, verbose);
    Basins.xdi = xd; Basins.ydi = yd;

    %% --- SUBFUNCTIONS ---

    % --- Set Thresholds ---
    function params = setThresholds(params, data, defaultValue, interactive, nb_ite, DEM2, FD2, A2, xd, yd)
        % Set a threshold value either to a default or interactively via ECDF plot.
        if ~interactive
            if verbose
                fprintf('Using default reference drainage area: %.2fkm²\n',round(defaultValue)*1e-6);
                params.Amin = params.Amin/DEM2.cellsize^2;
                return
            end
        end
    
        fprintf('Setting Amin interactively...\n');

        colors = lines(nb_ite); % Distinct colors for plot
    
        thresholds = 10.^(linspace(log10(params.A0/DEM2.cellsize^2), log10(max(data)), nb_ite+2));
        thresholds = thresholds(2:end-1); % Remove extremes
    
        % Create interactive threshold selection plot
        f = figure('Name', sprintf('Select reference drainage area for basin extraction'), 'Position', [100, 100, 800, 600]);

        sgtitle('Adjust slider to select suitable reference drainage area for basin extraction; confirm selection when satisfied', 'FontSize', 14);
    
        % ECDF subplot
        ax1 = axes('Parent', f, 'Position', [0.05, 0.25, 0.45, 0.65]); hold on;
        [x, y] = ecdf(data);
        plot(ax1, x, y, 'k-', 'LineWidth', 1.5);
        thresholdPlots = gobjects(1, nb_ite);
        thresholdText = gobjects(1, nb_ite);
        for i3 = 1:nb_ite
            idx = find(y >= thresholds(i3), 1);
            thresholdPlots(i3) = plot(ax1, x(idx), y(idx), 'o', 'Color', colors(i3,:), 'MarkerSize', 8, 'LineWidth', 1.5);
            thresholdText(i3) = text(ax1, x(idx), y(idx), sprintf(' %.2fkm²', thresholds(i3) * DEM2.cellsize^2 * 1e-6), ...
                'Color', colors(i3,:), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        end
        title(ax1, sprintf('ECDF with Possible Amin Thresholds'));
        ylabel(ax1, 'Cumulative Probability');
    
        % Network subplot with hillshade
        ax2 = axes('Parent', f, 'Position', [0.55, 0.25, 0.45, 0.65]);
        RGB = imageschs(DEM2, [], 'colormap', [.3 .3 .3], 'colorbar', false);
        [~, R] = GRIDobj2im(DEM2);
        imshow(flipud(RGB), R, 'Parent', ax2);
        hold(ax2, 'on');
        set(ax2, 'YDir', 'normal');
    
        % Pre-plot network options
        networkPlots = gobjects(1, nb_ite);
        try
            for i3 = 1:nb_ite
                Si = STREAMobj(FD2, (params.A0/DEM2.cellsize^2 < A2) & (A2 < thresholds(i3)));
                ixs = streampoi(Si, 'outlet', 'ix');
                D = drainagebasins(FD2, ixs); D.Z = double(D.Z);
                D.Z(D.Z==0) = nan;
                h = imagesc(D, 'Parent', ax2, 'Visible', 'off');
                h.AlphaData = 0.5;
                colormap(ax2, [0.3 0.3 0.3; lines()]);
                networkPlots(i3) = h;
            end
        catch
            error('Memory issue. Increase params.res in DefineParams function to handle coarser DEM instead.')
        end

        networkPlots(1).Visible = 'on';
        scatter(ax2, xd, yd, 'k.');
        title(ax2, sprintf('Basins outline for specific reference drainage area'));
    
        % Add slider
        sld = uicontrol('Style', 'slider', 'Min', 1, 'Max', nb_ite, 'Value', 1, ...
            'Units', 'normalized', 'Position', [0.2, 0.05, 0.6, 0.05], ...
            'SliderStep', [1/(nb_ite-1) 1/(nb_ite-1)]);
    
        % Display current threshold
        txt = uicontrol('Style', 'text', 'Units', 'normalized', ...
            'Position', [0.85, 0.05, 0.1, 0.05], 'String', num2str(thresholds(1)), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
        % Confirm button
        uicontrol('Style', 'pushbutton', 'String', 'Confirm', ...
            'Units', 'normalized', 'Position', [0.4, 0.01, 0.2, 0.05], ...
            'Callback', @(~,~) uiresume(f));
    
        % Slider listener
        addlistener(sld, 'ContinuousValueChange', @(src, evt) updatePlot(thresholdPlots, thresholdText, networkPlots, sld, txt));
        updatePlot(thresholdPlots, thresholdText, networkPlots, sld, txt);
        uiwait(f);
    
        % Set final threshold
        choiceNum = round(sld.Value);
        params.Amin = round(thresholds(choiceNum));
        fprintf('Set Amin to: %.2fkm²\n', round(thresholds(choiceNum) * DEM2.cellsize^2 *1e-6));
        close(f);
    
        % Update plot callback
        function updatePlot(thresholdPlots, thresholdText, networkPlots, sld, txt)
            choiceNum = round(sld.Value);
            selectedThreshold = thresholds(choiceNum);
            for i2 = 1:nb_ite
                if i2 == choiceNum
                    set(thresholdPlots(i2), 'MarkerSize', 12, 'LineWidth', 2);
                    set(thresholdText(i2), 'FontWeight', 'bold');
                    set(networkPlots(i2), 'Visible', 'on');
                else
                    set(thresholdPlots(i2), 'MarkerSize', 8, 'LineWidth', 1.5);
                    set(thresholdText(i2), 'FontWeight', 'normal');
                    set(networkPlots(i2), 'Visible', 'off');
                end
            end
            txt.String = sprintf('%.2f km²', selectedThreshold * DEM2.cellsize^2 * 1e-6);
            drawnow;
            pause(0.01);
        end
    end

    % --- Extract Basins ---
    function Basins = extractBasins(FD2, Si, xdd, ydd, ixdd, DEM2, verbose)
        % Extract drainage basins along the divide.
        di = [0; cumsum(sqrt(diff(xdd).^2 + diff(ydd).^2))];
        xys = streampoi(Si, 'outlet', 'xy');
        ixs = coord2ind(DEM2, xys(:, 1), xys(:, 2));
        [~, xDB, yDB] = drainagebasins(FD2, ixs);
        ixDB = coord2ind(DEM2, xDB, yDB);

        % Initialize basin data structure
        Basins = struct();
        Basins.DEM2 = DEM2;
        Basins.xd = xdd; Basins.yd = ydd;
        Basins.DEMi = cell(0);
        Basins.ixB = [];
        Basins.side = [];
        Basins.ixD = cell(0);
        Basins.di = cell(0);

        % Initialize progress bar
        if verbose
            PB = ProgressBar(numel(ixDB), 'taskname', 'Extract subcatchments...', 'ui', 'cli');
            figure; imageschs(DEM2, [], 'colormap', [0.7 0.7 0.7], 'colorbar', false); hold on;
            title('Drainage Basins Along Divide');
        end

        for i1 = 1:numel(ixDB)
            D = dependencemap(FD2, ixDB(i1));
            DEMii = crop(DEM2, D, NaN);
            [xb, yb] = getcoordinates(DEMii, 'matrix');
            xb(isnan(DEMii.Z)) = []; xb = xb(:);
            yb(isnan(DEMii.Z)) = []; yb = yb(:);

            r = DEM2.cellsize*2;
            idp = [];
            n = 0;
            m = 10;
            
            while 1
                ca = rangesearch([xb, yb], [xdd, ydd], r);
                idx = vertcat(find(cellfun(@(c) ~isempty(c), ca)));
                if isempty(idx), break, end
                if isequal(idx, idp)
                    n = n+1;
                else
                    n = 0; idp = idx;
                end
                xii = xdd(idx); yii = ydd(idx);
                h = all(abs(diff(yii))<1e-12);
                v = all(abs(diff(xii))<1e-12);
                if ~(h||v) || n>=m, break, end
                r = r*1.5;
            end
            
            if any(idx)
                Basins.DEMi{end+1} = DEMii;
                Basins.ixB(end+1) = ixDB(i1);
            
                ii = ixdd(idx);
                Basins.ixD{end+1} = ii;
                Basins.di{end+1} = di(idx);
                [xii, yii] = ind2coord(DEM2, ii);
                [xb, yb] = getoutline(DEMii, true);
                cx = mean(xb, 'omitnan'); cy = mean(yb, 'omitnan');
            
                isSide1 = zeros(length(xii),1);            
                if numel(xii) > 2
                    for i2 = 2:length(xii)-1
                        x1 = xii(i2-1); y1 = yii(i2-1);
                        x2 = xii(i2+1); y2 = yii(i2+1);
                        isSide1(i2) = ((x2-x1)*(cy-y1) - (y2-y1)*(cx-x1)) >= 0;
                    end
                    isSide1 = mode(isSide1(2:end-1));
                else
                    isSide1 = ((cx-mean(xii)) + (cy-mean(yii))) < 0;
                end
            
                if verbose
                    scatter(xii, yii, 'ko'); hold on;
                    if isSide1
                        Basins.side(end+1) = 1;
                        plot(xb, yb, 'b'); hold on;
                    else
                        Basins.side(end+1) = 2;
                        plot(xb, yb, 'r'); hold on;
                    end
                    drawnow;
                else
                    Basins.side(end+1) = isSide1 * 1 + ~isSide1 * 2;
                end
            end
            if verbose
                count(PB);
            end
        end


        % Sort basins by distance
        flat_di = cellfun(@(x) x(1), Basins.di);
        [~, sort_idx] = sort(flat_di);
        Basins.DEMi = Basins.DEMi(sort_idx);
        Basins.ixB = Basins.ixB(sort_idx);
        Basins.side = Basins.side(sort_idx);
    end
end