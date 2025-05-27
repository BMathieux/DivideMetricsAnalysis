function ixD = SelectDivide(D, DEM, varargin)
    % This function selects divide segments from a DIVIDEobj instance by allowing the
    % user to interactively specify divide node locations on a topographic map.
    % It converts the DIVIDEobj into a graph object to extract the shortest path between
    % the selected points, modified from the TopoToolbox function STREAMobj/modify
    % (Schwanghart, 2019).
    %
    % Inputs:
    % - D (DIVIDEobj): Instance of DIVIDEobj containing the divide network.
    % - DEM (GRIDobj): Digital elevation model used for visualization.
    % - varargin: Optional arguments:
    %   - 'verbose' (logical, default: true): Display progress messages and progress bar.
    %   - 'mode' (string, default: 'single'): Selection mode, either 'single'
    %     (stops after one selection) or 'multi' (allows multiple selections).
    %
    % Output:
    % - ixD (matrix): A two-column matrix where:
    %   - 1st column: Linear indices of the selected divide segments in the DIVIDEobj.
    %   - 2nd column: Iteration number identifying each selected segment group.
    %
    % Example:
    %   DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
    %   FD = FLOWobj(DEM, 'preprocess', 'c');
    %   ST = STREAMobj(FD, flowacc(FD) > 1000);
    %   D = DIVIDEobj(FD, ST);
    %   D = removeshortdivides(D, FD, 500);
    %   ixD = SelectDivide(D, DEM, 'verbose', true, 'mode', 'single');
    %
    % See also: DIVIDEobj, STREAMobj/modify, divorder
    %
    % Author: Bastien Mathieux (b.mathieux[at]gmail.com)
    % Date: 8 February 2023

    % --- Parse Inputs ---
    p = inputParser;
    p.FunctionName = 'SelectDivide';
    addRequired(p, 'D', @(x) isa(x, 'DIVIDEobj'));
    addRequired(p, 'DEM', @(x) isa(x, 'GRIDobj'));
    addParameter(p, 'verbose', true, @islogical);
    addParameter(p, 'mode', 'single', @(x) ischar(x) && any(strcmpi(x, {'single', 'multi'})));
    parse(p, D, DEM, varargin{:});
    useVerbose = p.Results.verbose;
    mode = p.Results.mode;

    % --- Prepare Divide Network ---
    if ~any(D.ordertype)
        D = divorder(D, 'topo'); % Order divide network for visualization
    end
    ixD = [];
    answ = 0;
    ite = 0;

    [dx, dy] = ind2coord(D, D.IX);

    % --- Create graph object ---

    ix = D.IX(:); ix = ix(~isnan(ix));
    [uix, ~, idxu] = unique(ix);
    
    % cearte graph object
    s = []; t = [];
    b = isnan(D.IX(:)); ix_full = D.IX(:);
    i1 = [1; find(b)+1];
    i2 = [find(b)-1; numel(ix_full)];
    for k = 1:numel(i1)
        seg = ix_full(i1(k):i2(k));
        seg = seg(~isnan(seg));
        if numel(seg)<2, continue; end
        [~, segmap] = ismember(seg, uix);
        s = [s; segmap(1:end-1); segmap(2:end)];
        t = [t; segmap(2:end);   segmap(1:end-1)];
    end
    G = graph(s, t);

    % --- Select Divide Segments ---
    if useVerbose
        disp('Starting interactive divide selection');
    end

    f = figure;
    imageschs(DEM, [], 'colormap', [0.7 0.7 0.7], 'colorbar', false); hold on;
    scatter(dx, dy, 10, '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
    plot(D, 'color', [0 0 0]);
    axis image; box on;

    ax = gca;
    lims = axis;
    xlim(ax, [lims(1) - (lims(2) - lims(1))/20, lims(2) + (lims(2) - lims(1))/20]);
    ylim(ax, [lims(3) - (lims(4) - lims(3))/20, lims(4) + (lims(4) - lims(3))/20]);

    while answ == 0
        ixpath = [];
        ixe = 0; ixi = [];
        hpath = [];

        % Prompt for upper reach (green)
        title('Set last (green) divide node location');
        hpstart = impoint('PositionConstraintFcn', @(pos) getnearest(pos));
        setColor(hpstart, [0 1 0]);
        addNewPositionCallback(hpstart, @(pos) updatePath(pos, 'start'));
        setPosition(hpstart, getPosition(hpstart));

        % Prompt for lower reach (red)
        title('Set first (red) divide node location');
        hpend = impoint('PositionConstraintFcn', @(pos) getnearestend(pos));
        setColor(hpend, [1 0 0]);
        addNewPositionCallback(hpend, @(pos) updatePath(pos, 'end'));
        setPosition(hpend, getPosition(hpend));

        % Wait for confirmation
        title(["Move divide nodes locations and press any key to extract reach" ...
               "If valley bottom connections, shrink DIVIDEobj"])

        set(gcf, 'WindowKeyPressFcn', @(k, l) uiresume);
        uiwait;

        delete(hpstart); delete(hpend);
        ite = ite + 1;

        % Store selected path
        if ~isempty(ixpath)
            ixD = [ixD; [ixpath; nan], ones(numel(ixpath)+1, 1) * ite];
        end

        if strcmpi(mode, 'single')
            if ~isempty(ixpath)
                [x, y] = ind2coord(D, ixpath);
                hpath = plot(x, y, 'color', [0.3 0.3 1], 'LineWidth', 1.51);
                hold on;
            end
            answ = 1;
        else
            answer = questdlg('Would you like to keep?', '', 'Yes', 'No', 'Yes');
            if strcmp(answer, 'Yes') && ~isempty(ixpath)
                [x, y] = ind2coord(D, ixpath);
                hpath = plot(x, y, 'color', [0.3 0.3 1], 'LineWidth', 1.51);
                hold on;
            else
                answ = 1;
            end
        end
    end

    close(f);

    %% --- SUBFUNCTIONS ---

    % --- Get Nearest ---
    function posn = getnearest(pos)
        % Snap to nearest divide point for start position.
        [~, ixi] = min((dx - pos(1)).^2 + (dy - pos(2)).^2);
        posn = [dx(ixi) dy(ixi)];
    end

    % --- Get Nearest End ---
    function posn = getnearestend(pos)
        % Snap to nearest divide point for end position.
        [~, ixe] = min((dx - pos(1)).^2 + (dy - pos(2)).^2);
        posn = [dx(ixe) dy(ixe)];
    end

    function updatePath(pos, pointType)
        % pos = [x, y] clicked
        [~, new_idx] = min((dx - pos(1)).^2 + (dy - pos(2)).^2); % dx, dy: coordinates of D.IX
        linear_idx = D.IX(new_idx); % This is a DEM index
    
        % Find the corresponding node number in the graph
        node_in_graph = find(uix == linear_idx);
    
        if strcmp(pointType, 'start')
            ixi = node_in_graph;
        else % 'end'
            ixe = node_in_graph;
        end
    
        % Only if both points are set:
        if ixi ~= 0 && ixe ~= 0
            P = shortestpath(G, ixi, ixe); % These are graph node numbers
            ixpath = uix(P);               % Convert graph node numbers back to DEM linear indices
            [x, y] = ind2coord(D, ixpath);
            if ishandle(hpath)
                set(hpath, 'XData', x, 'YData', y);
            else
                hold on;
                hpath = plot(x, y, 'r', 'LineWidth', 1.5);
                hold on;
            end
        end
    end
end