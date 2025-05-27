function [DTable, MTable] = CalcDivideTable(HMetrics, BData, DEMg, info, params, varargin)
    % Calculate divide metrics for paired basins across a drainage divide.
    % This function computes topographic and fluvial metrics (gradient, elevation, 
    % relief, ksn, chi) for each basin, then maps these to paired basins on either 
    % side of a divide. Optionally includes hilltop curvature, denudation, and 
    % migration rates if RunCHt is true. Outputs are structured as tables (DTable 
    % for detailed metrics, MTable for mean values).
    % Note: Uncertainties are calculated as the range of bootstrapped means. Modify 
    % the calcMetric subfunction to use a different method.
    %
    % Inputs:
    %   HMetrics - Hillslope metrics (cell array, empty or NaN if RunCHt is false).
    %   BData    - Basin characteristics, including DEMs and divide coordinates.
    %   DEMg     - Global DEM for computing fluvial metrics (GRIDobj).
    %   info     - Geotiff metadata for coordinate projection.
    %   params   - Analysis parameters defined with DefineParams function.
    %   varargin - Optional settings: 'Parallel' (true), 'Verbose' (true), 'RunCHt' (true).
    %
    % Outputs:
    %   DTable - Detailed table with metrics per basin pair.
    %   MTable - Summary table with mean metrics per unique pair.

    % --- Parse Inputs ---
    p = inputParser;
    addParameter(p, 'Parallel', true, @islogical);
    addParameter(p, 'Verbose', true, @islogical);
    addParameter(p, 'RunCHt', true, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Initialize Parallel Pool ---
    if opts.Parallel
        if isempty(gcp('nocreate'))
            parpool();
        end
    end

    % --- Initialize DTable Structure ---
    DTable = initDTable(BData, HMetrics, params, opts.RunCHt);

    % --- Pair Basins ---
    [DTable.id1, DTable.id2] = pairBasins(DTable.xd, DTable.yd, BData, opts.Verbose);

    % --- Compute Metrics ---
    metrics = computeMetrics(DTable, BData, DEMg, HMetrics, params, opts);

    % --- Create DTable ---
    DTable = createDTable(DTable, metrics, info, opts.RunCHt, params.zb, params.mn);

    % --- Create MTable ---
    MTable = createMTable(DTable, opts.RunCHt);
end

%% --- SUBFUNCTIONS ---

% --- Initialize DTable ---
function DTable = initDTable(BData, HMetrics, params, RunCHt)
    % Initialize DTable with basin and divide data.
    DTable.B.ixB = BData.ixB;
    DTable.B.side = BData.side;
    DTable.B.ixD = BData.ixD;
    DTable.B.di = BData.di;
    DTable.xd = BData.xd;
    DTable.yd = BData.yd;
    DTable.d = [0; cumsum(sqrt(diff(BData.xd).^2 + diff(BData.yd).^2))];
    DTable.params.nbstrp = params.nbstrp;

    if RunCHt
        DTable.B.ws = HMetrics{1}.ws;
        for i = 1:numel(BData.ixB)
            DTable.B.Cht{i} = HMetrics{i}.Cht;
            DTable.B.Sht{i} = HMetrics{i}.Sht;
            DTable.B.Rht{i} = HMetrics{i}.Rht;
            DTable.B.Lht{i} = HMetrics{i}.Lht;
            DTable.B.xht{i} = HMetrics{i}.xht;
            DTable.B.fpq{i} = HMetrics{i}.fpq;
            DTable.B.fpp{i} = HMetrics{i}.fpp;
        end
        DTable.ws = HMetrics{1}.ws;
        DTable.params.K = params.K;
        DTable.params.rp = params.rp;
        DTable.params.wH = params.wH;
    end
end

% --- Pair Basins ---
function [id1, id2] = pairBasins(xd, yd, BData, verbose)
    % Pair basins on either side of the divide by finding closest points.
    if verbose
        fprintf('Pairing divide nodes to neighboring basins...\n');
        tic
    end
    n = numel(xd);
    
    % Validate divide points
    valid = ~isnan(xd) & ~isnan(yd);
    xd = xd(valid); yd = yd(valid);
    n = numel(xd);
    if n == 0
        fprintf('No valid divide points found.\n');
        id1 = []; id2 = [];
        return;
    end

    % Initialize basin outline arrays
    B1 = []; B2 = [];
    for i = 1:numel(BData.DEMi)
        try
            [xb, yb] = getoutline(BData.DEMi{i}, true);
            if numel(xb) < 2 || numel(yb) < 2
                continue; % Skip basins with insufficient outline points
            end
            xb = interp1(1:length(xb), xb, linspace(1, length(xb), n), 'linear');
            yb = interp1(1:length(yb), yb, linspace(1, length(yb), n), 'linear');
            if BData.side(i) == 1
                B1 = [B1; [xb(:), yb(:), repmat(i, numel(xb), 1)]];
            else
                B2 = [B2; [xb(:), yb(:), repmat(i, numel(xb), 1)]];
            end
        catch
            continue; % Skip basins with outline errors
        end
    end

    % Check if B1 or B2 is empty
    if isempty(B1) || isempty(B2)
        id1 = []; id2 = [];
        return;
    end

    id1 = zeros(n, 1); id2 = zeros(n, 1);
    for i = 1:n
        % Find closest point in B1
        [~, idx1] = min(pdist2([B1(:,1), B1(:,2)], [xd(i), yd(i)]));
        if idx1 <= size(B1, 1)
            id1(i) = B1(idx1, 3);
        else
            id1(i) = 0; % Invalid index
        end

        % Find closest point in B2
        [~, idx2] = min(pdist2([B2(:,1), B2(:,2)], [xd(i), yd(i)]));
        if idx2 <= size(B2, 1)
            id2(i) = B2(idx2, 3);
        else
            id2(i) = 0; % Invalid index
        end
    end

    % Remove invalid pairs
    valid_pairs = id1 > 0 & id2 > 0;
    id1 = id1(valid_pairs);
    id2 = id2(valid_pairs);

    if verbose
        toc
    end
end

% --- Compute Metrics ---
function metrics = computeMetrics(DTable, BData, DEMg, HMetrics, params, opts)
    % Compute topographic and fluvial metrics for each basin.
    n = numel(DTable.xd); % Total number of paired divide nodes
    nb = numel(BData.ixB); % Number of basins

    % Pre-allocate basin_metrics as a cell array of structures
    basin_metrics = cell(nb, 1);
    for i = 1:nb
        basin_metrics{i} = struct('G', NaN, 'Gu', NaN, 'Z', NaN, 'H', NaN, 'Hu', NaN, ...
                                  'ksn', NaN, 'ksnu', NaN, 'chi', NaN(numel(params.zb), 1));
        if opts.RunCHt
            basin_metrics{i}.Ct = NaN; basin_metrics{i}.Ctu = NaN;
            basin_metrics{i}.Rt = NaN; basin_metrics{i}.Rtu = NaN;
            basin_metrics{i}.St = NaN; basin_metrics{i}.Stu = NaN;
            basin_metrics{i}.Lt = NaN; basin_metrics{i}.Ltu = NaN;
            basin_metrics{i}.E = NaN; basin_metrics{i}.Eu = NaN;
        end
    end

    % Compute fluvial metrics
    fluv = computeFluvialMetrics(DEMg, params.A0, params.mn, params.zb, opts.Verbose);


    % Compute metrics for each basin
    if opts.Verbose
        fprintf('Computing metrics for %d basins...\n', nb);
    end
    if opts.Parallel
        parfor i = 1:nb
            DEMi = BData.DEMi{i};
            HMetricsi = HMetrics{i};
            basin_metrics{i} = computeSingleBasinMetric(DEMi, DEMg, HMetricsi, fluv, params, opts);
        end
    else
        for i = 1:nb
            DEMi = BData.DEMi{i};
            if opts.RunCHt, HMetricsi = HMetrics{i}; else, HMetricsi = nan; end
                basin_metrics{i} = computeSingleBasinMetric(DEMi, DEMg, HMetricsi, fluv, params, opts);
        end
    end

    if opts.Verbose
        fprintf('Basin metrics computation completed.\n');
    end

    % Initialize pair-level metric arrays
    G1 = nan(1, n); G2 = nan(1, n); G1u = nan(1, n); G2u = nan(1, n);
    Z1 = nan(1, n); Z2 = nan(1, n);
    H1 = nan(1, n); H2 = nan(1, n); H1u = nan(1, n); H2u = nan(1, n);
    dG = nan(1, n); dGu = nan(1, n); dZ = nan(1, n); dH = nan(1, n); dHu = nan(1, n);
    ksn1 = nan(1, n); ksn2 = nan(1, n); ksn1u = nan(1, n); ksn2u = nan(1, n);
    dksn = nan(1, n); dksnu = nan(1, n);
    chi1 = nan(numel(params.zb), n); chi2 = nan(numel(params.zb), n); dchi = nan(numel(params.zb), n);

    if opts.RunCHt
        Ct1 = nan(1, n); Ct2 = nan(1, n); Ct1u = nan(1, n); Ct2u = nan(1, n);
        Rt1 = nan(1, n); Rt2 = nan(1, n); Rt1u = nan(1, n); Rt2u = nan(1, n);
        St1 = nan(1, n); St2 = nan(1, n); St1u = nan(1, n); St2u = nan(1, n);
        Lt1 = nan(1, n); Lt2 = nan(1, n); Lt1u = nan(1, n); Lt2u = nan(1, n);
        E1 = nan(1, n); E2 = nan(1, n); E1u = nan(1, n); E2u = nan(1, n);
        dE = nan(1, n); dEu = nan(1, n); MIG = nan(1, n); MIGu = nan(1, n);
    end

    % Map basin metrics to pairs
    for i = 1:n
        id1 = DTable.id1(i);
        id2 = DTable.id2(i);
        G1(i) = basin_metrics{id1}.G; G1u(i) = basin_metrics{id1}.Gu;
        G2(i) = basin_metrics{id2}.G; G2u(i) = basin_metrics{id2}.Gu;
        Z1(i) = basin_metrics{id1}.Z;
        Z2(i) = basin_metrics{id2}.Z;
        H1(i) = basin_metrics{id1}.H; H1u(i) = basin_metrics{id1}.Hu;
        H2(i) = basin_metrics{id2}.H; H2u(i) = basin_metrics{id2}.Hu;
        ksn1(i) = basin_metrics{id1}.ksn; ksn1u(i) = basin_metrics{id1}.ksnu;
        ksn2(i) = basin_metrics{id2}.ksn; ksn2u(i) = basin_metrics{id2}.ksnu;
        chi1(:,i) = basin_metrics{id1}.chi;
        chi2(:,i) = basin_metrics{id2}.chi;
        dG(i) = G2(i) - G1(i);
        dGu(i) = sqrt(G1u(i)^2 + G2u(i)^2);
        dZ(i) = Z2(i) - Z1(i);
        dH(i) = H2(i) - H1(i);
        dHu(i) = sqrt(H1u(i)^2 + H2u(i)^2);
        dksn(i) = ksn2(i) - ksn1(i);
        dksnu(i) = sqrt(ksn1u(i)^2 + ksn2u(i)^2);
        dchi(:,i) = chi2(:,i) - chi1(:,i);

        if opts.RunCHt
            Ct1(i) = basin_metrics{id1}.Ct; Ct1u(i) = basin_metrics{id1}.Ctu;
            Ct2(i) = basin_metrics{id2}.Ct; Ct2u(i) = basin_metrics{id2}.Ctu;
            Rt1(i) = basin_metrics{id1}.Rt; Rt1u(i) = basin_metrics{id1}.Rtu;
            Rt2(i) = basin_metrics{id2}.Rt; Rt2u(i) = basin_metrics{id2}.Rtu;
            St1(i) = basin_metrics{id1}.St; St1u(i) = basin_metrics{id1}.Stu;
            St2(i) = basin_metrics{id2}.St; St2u(i) = basin_metrics{id2}.Stu;
            Lt1(i) = basin_metrics{id1}.Lt; Lt1u(i) = basin_metrics{id1}.Ltu;
            Lt2(i) = basin_metrics{id2}.Lt; Lt2u(i) = basin_metrics{id2}.Ltu;
            E1(i) = basin_metrics{id1}.E; E1u(i) = basin_metrics{id1}.Eu;
            E2(i) = basin_metrics{id2}.E; E2u(i) = basin_metrics{id2}.Eu;
            dE(i) = E2(i) - E1(i);
            dEu(i) = sqrt(E1u(i)^2 + E2u(i)^2);
            MIG(i) = dE(i) / (G1(i) + G2(i));
            MIGu(i) = abs(MIG(i)) * sqrt((dEu(i)/dE(i))^2 + (dGu(i)/(G1(i) + G2(i)))^2);
        end
    end

    % Assemble metrics structure
    metrics = struct('G1', G1, 'G2', G2, 'G1u', G1u, 'G2u', G2u, ...
                     'Z1', Z1, 'Z2', Z2, ...
                     'H1', H1, 'H2', H2, 'H1u', H1u, 'H2u', H2u, ...
                     'dG', dG, 'dGu', dGu, 'dZ', dZ, 'dH', dH, 'dHu', dHu, ...
                     'ksn1', ksn1, 'ksn2', ksn2, 'ksn1u', ksn1u, 'ksn2u', ksn2u, ...
                     'dksn', dksn, 'dksnu', dksnu, ...
                     'chi1', chi1, 'chi2', chi2, 'dchi', dchi);
    if opts.RunCHt
        metrics.Ct1 = Ct1; metrics.Ct2 = Ct2; metrics.Ct1u = Ct1u; metrics.Ct2u = Ct2u;
        metrics.Rt1 = Rt1; metrics.Rt2 = Rt2; metrics.Rt1u = Rt1u; metrics.Rt2u = Rt2u;
        metrics.St1 = St1; metrics.St2 = St2; metrics.St1u = St1u; metrics.St2u = St2u;
        metrics.Lt1 = Lt1; metrics.Lt2 = Lt2; metrics.Lt1u = Lt1u; metrics.Lt2u = Lt2u;
        metrics.E1 = E1; metrics.E2 = E2; metrics.E1u = E1u; metrics.E2u = E2u;
        metrics.dE = dE; metrics.dEu = dEu; metrics.MIG = MIG; metrics.MIGu = MIGu;
    end
    

end

% --- Compute Single Basin Metric ---
function single_metrics = computeSingleBasinMetric(DEM, DEMg, HMetricsi, fluv, params, opts)
    % Compute metrics for a single basin.
    single_metrics = struct('G', NaN, 'Gu', NaN, 'Z', NaN, 'H', NaN, 'Hu', NaN, ...
                            'ksn', NaN, 'ksnu', NaN, 'chi', NaN(numel(params.zb), 1));
    if opts.RunCHt
        single_metrics.Ct = NaN; single_metrics.Ctu = NaN;
        single_metrics.Rt = NaN; single_metrics.Rtu = NaN;
        single_metrics.St = NaN; single_metrics.Stu = NaN;
        single_metrics.Lt = NaN; single_metrics.Ltu = NaN;
        single_metrics.E = NaN; single_metrics.Eu = NaN;
    end

    try
        % Compute topographic metrics
        x = gradient8(DEM); x = x.Z(:); x = x(x > 0);
        [single_metrics.G, single_metrics.Gu] = calcMetric(x, @mean, params.nbstrp);
        x = DEM.Z(:); x = x(x > 0);
        single_metrics.Z = mean(x);
        x = localtopography(DEM, params.wH); x = x.Z(:); x = x(x > 0);
        [single_metrics.H, single_metrics.Hu] = calcMetric(x, @mean, params.nbstrp);

        % Compute fluvial metrics
        [xdem,ydem] = getcoordinates(DEM,'matrix');
        xdem(isnan(DEM.Z)) = []; ydem(isnan(DEM.Z)) = [];
        [~,idx_fluv,~] = intersect(coord2ind(DEMg,fluv.xksn, fluv.yksn), ...
                                   coord2ind(DEMg,xdem,ydem));

        x = fluv.ksn(idx_fluv); x = x(~isnan(x));
        [single_metrics.ksn, single_metrics.ksnu] = calcMetric(x, @mean, params.nbstrp);
        for j = 1:numel(params.zb)
            [~,idx_chi,~] = intersect(coord2ind(DEMg,fluv.xchi{j}, fluv.ychi{j}), ...
                                   coord2ind(DEMg,xdem,ydem));
            valid = idx_chi(~isnan(fluv.chi{j}(idx_chi)));
            if ~isempty(valid)
                single_metrics.chi(j) = max(fluv.chi{j}(valid));
            end
        end

        % Compute hillslope metrics
        if opts.RunCHt
            [~, ws_idx] = min(abs(HMetricsi.ws - params.wl));
            x = HMetricsi.Cht(:, ws_idx); x = x(~isnan(x));
            [single_metrics.Ct, single_metrics.Ctu] = calcMetric(x, @mean, params.nbstrp);
            x = HMetricsi.Rht(:, ws_idx); x = x(~isnan(x));
            [single_metrics.Rt, single_metrics.Rtu] = calcMetric(x, @mean, params.nbstrp);
            x = HMetricsi.Sht(:, ws_idx); x = x(~isnan(x));
            [single_metrics.St, single_metrics.Stu] = calcMetric(x, @mean, params.nbstrp);
            x = HMetricsi.Lht(:, ws_idx); x = x(~isnan(x));
            [single_metrics.Lt, single_metrics.Ltu] = calcMetric(x, @mean, params.nbstrp);
            single_metrics.E = -params.K * single_metrics.Ct * params.rp;
            single_metrics.Eu = abs(single_metrics.E) * sqrt((params.dK/params.K)^2 + (single_metrics.Ctu/single_metrics.Ct)^2);
        end
    catch
    end
end

% --- Compute Fluvial Metrics ---
function fluv = computeFluvialMetrics(DEMg, A0, mn, zb, verbose)
    % Compute ksn and chi for fluvial analysis.
    if verbose
        fprintf('Extracting fluvial metrics...\n');
        tic
    end
    FDg = FLOWobj(DEMg, 'preprocess', 'c');
    Ag = flowacc(FDg);
    Sg = STREAMobj(FDg, Ag.*DEMg.cellsize^2 > A0);
    fluv.ksn = ksn(Sg, DEMg, Ag, mn, 100);
    fluv.xksn = Sg.x;
    fluv.yksn = Sg.y;
    fluv.chi = cell(numel(zb), 1);
    fluv.xchi = cell(numel(zb), 1);
    fluv.ychi = cell(numel(zb), 1);
    for i = 1:numel(zb)
        DEMi = DEMg;
        DEMi.Z(DEMi.Z < zb(i)) = nan;
        FDi = FLOWobj(DEMi, 'preprocess', 'c');
        Ai = flowacc(FDi);
        Si = STREAMobj(FDi, Ai.*DEMi.cellsize^2 > A0);
        fluv.chi{i} = chitransform(Si, Ai, 'a0', A0, 'mn', mn);
        fluv.xchi{i} = Si.x;
        fluv.ychi{i} = Si.y;
    end
    if verbose
        toc
    end
end

% --- Calculate Metric with Uncertainty ---
function [val, unc] = calcMetric(data, func, nbstrp)
    % Compute mean and uncertainty for a metric using bootstrapping.
    data = data(~isnan(data));
    if isempty(data) || numel(data) < 2
        val = nan;
        unc = nan;
    else
        x = bootstrp(nbstrp, func, data);
        val = mean(x);
        unc = max(x) - min(x);
    end
end

% --- Create DTable ---
function DTable = createDTable(DTable, metrics, info, RunCHt, zb, mn)
    % Create detailed metrics table with dynamic chi columns.
    [lat, lon] = projinv(info.CoordinateReferenceSystem, DTable.xd(:), DTable.yd(:));
    metrics.chi1 = transpose(metrics.chi1);
    metrics.chi2 = transpose(metrics.chi2);
    metrics.dchi = transpose(metrics.dchi);
    [~, ~, pair_num] = unique([DTable.id1(:), DTable.id2(:)], 'rows', 'stable');

    ksn_unit = sprintf('m^{%.2f}', round(2 * mn, 2));

    cols = {'Pair_number', pair_num; ...
            'distance_along_divide (m)', DTable.d(:); ...
            'Basin_ID1', DTable.id1(:); ...
            'Lat1 (°)', lat(:); ...
            'Lon1 (°)', lon(:); ...
            'G1 (m/m)', metrics.G1(:); ...
            'Z1 (m)', metrics.Z1(:); ...
            'H1 (m)', metrics.H1(:); ...
            ['ksn1 (' ksn_unit ')'], metrics.ksn1(:); ...
            ['ksn1_unc (' ksn_unit ')'], metrics.ksn1u(:); ...
            'Basin_ID2', DTable.id2(:); ...
            'Lat2 (°)', lat(:); ...
            'Lon2 (°)', lon(:); ...
            'G2 (m/m)', metrics.G2(:); ...
            'Z2 (m)', metrics.Z2(:); ...
            'H2 (m)', metrics.H2(:); ...
            ['ksn2 (' ksn_unit ')'], metrics.ksn2(:); ...
            ['ksn2_unc (' ksn_unit ')'], metrics.ksn2u(:); ...
            'dG (m/m)', metrics.dG(:); ...
            'dZ (m)', metrics.dZ(:); ...
            'dH (m)', metrics.dH(:); ...
            ['dksn (' ksn_unit ')'], metrics.dksn(:)};

    for i = 1:numel(zb)
        label = sprintf('%.0f', zb(i));
        cols = [cols; ...
                {['chi1-' label '(m)'], metrics.chi1(:,i)}; ...
                {['chi2-' label '(m)'], metrics.chi2(:,i)}; ...
                {['dchi-' label '(m)'], metrics.dchi(:,i)}];
    end

    if RunCHt
        cols = [cols; ...
                {'Ct1 (1/m)', metrics.Ct1(:); ...
                 'Ct1_unc (1/m)', metrics.Ct1u(:); ...
                 'E1 (m/yr)', metrics.E1(:); ...
                 'E1_unc (m/yr)', metrics.E1u(:); ...
                 'Ct2 (1/m)', metrics.Ct2(:); ...
                 'Ct2_unc (1/m)', metrics.Ct2u(:); ...
                 'E2 (m/yr)', metrics.E2(:); ...
                 'E2_unc (m/yr)', metrics.E2u(:); ...
                 'dE (m/yr)', metrics.dE(:); ...
                 'dCt (1/m)', metrics.Ct2(:) - metrics.Ct1(:); ...
                 'MIG (m/yr)', metrics.MIG(:); ...
                 'MIG_unc (m/yr)', metrics.MIGu(:)}];
    end

    DTable = table(cols{:,2}, 'VariableNames', cols(:,1));
end

% --- Create MTable ---
function MTable = createMTable(DTable, RunCHt)
    % Create summary table with mean metrics per unique pair.
    unique_pairs = unique(DTable.Pair_number(~isnan(DTable.Pair_number)));
    MTable = table(unique_pairs, 'VariableNames', {'Pair_number'});

    % Extract ksn unit from ksn1 column name
    ksn_col = DTable.Properties.VariableNames{find(contains(DTable.Properties.VariableNames, 'ksn1 ('), 1)};
    ksn_unit = regexp(ksn_col, 'm\^{\d+\.\d+}', 'match', 'once');
    if isempty(ksn_unit)
        ksn_unit = 'm^{0.90}'; % Fallback unit
    end

    % Define value-uncertainty pairs (column names)
    unc_names = {...
        ['ksn1_unc (' ksn_unit ')'], ...
        ['ksn2_unc (' ksn_unit ')']};
    val_unc_pairs = {...
        ['ksn1 (' ksn_unit ')'], ['ksn1_unc (' ksn_unit ')']; ...
        ['ksn2 (' ksn_unit ')'], ['ksn2_unc (' ksn_unit ')']};
    if RunCHt
        unc_names = [unc_names, ...
            {'Ct1_unc (1/m)', 'E1_unc (m/yr)', 'Ct2_unc (1/m)', 'E2_unc (m/yr)', 'MIG_unc (m/yr)'}];
        val_unc_pairs = [val_unc_pairs; ...
            {'Ct1 (1/m)', 'Ct1_unc (1/m)'}; ...
            {'E1 (m/yr)', 'E1_unc (m/yr)'}; ...
            {'Ct2 (1/m)', 'Ct2_unc (1/m)'}; ...
            {'E2 (m/yr)', 'E2_unc (m/yr)'}; ...
            {'MIG (m/yr)', 'MIG_unc (m/yr)'}];
    end

    % Map uncertainty names to indices
    unc_indices = zeros(1, length(unc_names));
    for k = 1:length(unc_names)
        idx = find(strcmp(DTable.Properties.VariableNames, unc_names{k}), 1);
        if ~isempty(idx)
            unc_indices(k) = idx;
        end
    end
    % Map value-uncertainty pairs to indices
    val_unc_indices = zeros(size(val_unc_pairs));
    for k = 1:size(val_unc_pairs, 1)
        val_idx = find(strcmp(DTable.Properties.VariableNames, val_unc_pairs{k,1}), 1);
        unc_idx = find(strcmp(DTable.Properties.VariableNames, val_unc_pairs{k,2}), 1);
        if ~isempty(val_idx) && ~isempty(unc_idx)
            val_unc_indices(k,:) = [val_idx, unc_idx];
        end
    end

    % Process columns (indices 2 to end)
    for i = 2:size(DTable, 2)
        col_name = DTable.Properties.VariableNames{i};
        MTable.(col_name) = nan(size(unique_pairs));

        % Basin_ID1 (index 3) or Basin_ID2 (index 11)
        if i == 3 || i == 11
            for j = 1:length(unique_pairs)
                p = unique_pairs(j);
                idx = find(DTable.Pair_number == p, 1, 'first');
                if ~isempty(idx)
                    MTable.(col_name)(j) = DTable{idx, i};
                end
            end
        % distance_along_divide (index 2)
        elseif i == 2
            for j = 1:length(unique_pairs)
                p = unique_pairs(j);
                data = DTable{DTable.Pair_number == p, i};
                if ~isempty(data) && ~all(isnan(data))
                    MTable.(col_name)(j) = mean(data, 'omitnan');
                end
            end
        % Uncertainty columns
        elseif ismember(i, unc_indices)
            val_idx = val_unc_indices(val_unc_indices(:,2) == i, 1);
            for j = 1:length(unique_pairs)
                p = unique_pairs(j);
                idx = DTable.Pair_number == p;
                val_data = DTable{idx, val_idx};
                unc_data = DTable{idx, i};
                if ~isempty(val_data) && ~isempty(unc_data) && ~all(isnan(val_data)) && ~all(isnan(unc_data))
                    MTable.(col_name)(j) = (max(val_data + unc_data) - min(val_data - unc_data)) / 2;
                end
            end
        % Default: Compute mean (includes chi columns and others)
        else
            for j = 1:length(unique_pairs)
                p = unique_pairs(j);
                data = DTable{DTable.Pair_number == p, i};
                if ~isempty(data) && ~all(isnan(data))
                    MTable.(col_name)(j) = mean(data, 'omitnan');
                end
            end
        end
    end
end