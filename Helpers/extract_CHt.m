function OUT = extract_CHt(params, useParallel, verbose)
    % Extract hilltop curvature and flow path metrics along a drainage divide.
    % This function processes DEMs to compute curvature, slope, relief, and flow path 
    % lengths at specified divide points, with surface roughness filtering and steepest 
    % descent tracing. Outputs are structured as a single structure for further analysis.
    % Hilltop nodes with no valid flow paths are skipped.
    %
    % Inputs:
    %   params       - Analysis parameters defined by the DefineParams function.
    %   useParallel  - Enable parallel processing (logical, true/false).
    %   verbose      - Display progress messages (logical, true/false).
    %
    % Outputs:
    %   OUT          - Structure with hilltop curvature and flow path metrics.

    % --- Initialize Parallel Pool ---
    if useParallel
        if isempty(gcp('nocreate'))
            parpool();
        end
    end

    % --- Surface Roughness Filtering ---
    id = filterRoughness(params.DEM1, params.xd, params.yd, params.wsR, params.Sth, verbose);
    xd = params.xd(id); yd = params.yd(id);

    % --- Steepest Descent Tracing ---
    [fp_coords, fq_coords, lSHt, lRHt, lLHt, lxHt, lyHt] = traceDescent(params.DEM2, params.FD2, params.S2, xd, yd, params.Lth, verbose);

    % --- Extract Hilltop Curvature ---
    OUT = computeCurvature(params.DEM1, params.ws, params.Gth, xd, yd, lSHt, lRHt, lLHt, lxHt, lyHt, fp_coords, fq_coords, verbose);
end

%% --- SUBFUNCTIONS ---

% --- Filter Roughness ---
function id = filterRoughness(DEM1, xd, yd, wsR, Sth, verbose)
    % Filter hilltop nodes based on surface roughness to identify valid hilltops.
    if verbose
        PB = ProgressBar(numel(xd)+1, 'taskname', 'Surface roughness filtering...', 'ui', 'cli');
        count(PB);
    end

    [xA, yA] = getcoordinates(DEM1, 'matrix');
    Z = DEM1.Z;
    id = [];

    for i = 1:numel(xd)
        if verbose
            count(PB);
        end
        [xi, yi] = coord2sub(DEM1, xd(i), yd(i));
        if isnan(xi) || isnan(yi) || isnan(Z(xi, yi))
            continue;
        end

        % Extract window for roughness analysis
        try
            Zwin = Z(xi-wsR:xi+wsR, yi-wsR:yi+wsR);
            xwin = xA(xi-wsR:xi+wsR, yi-wsR:yi+wsR);
            ywin = yA(xi-wsR:xi+wsR, yi-wsR:yi+wsR);
        catch
            continue;
        end

        if range(Zwin(:)) < 1e-6
            continue;
        end

        % Fit quadratic surface
        try
            xwin = double(xwin); ywin = double(ywin); Zwin = double(Zwin);
            xwin = (xwin - mean(xwin(:))) / std(xwin(:));
            ywin = (ywin - mean(ywin(:))) / std(ywin(:));
            Zwin = (Zwin - mean(Zwin(:))) / std(Zwin(:));
            fitr = fit([xwin(:), ywin(:)], Zwin(:), 'poly22');
            a = fitr.p20; b = fitr.p02; c = fitr.p11; d = fitr.p10; e = fitr.p01;
        catch
            continue;
        end

        % Compute normal vectors
        nx = -2*a*xwin - c*ywin - d;
        ny = -2*b*ywin - c*xwin - e;
        nz = ones(size(nx));
        norms = sqrt(nx.^2 + ny.^2 + nz.^2);
        nx = nx ./ norms; ny = ny ./ norms; nz = nz ./ norms;

        % Orientation matrix and eigenvalues
        T = zeros(3, 3);
        for k = 1:numel(nx)
            n = [nx(k), ny(k), nz(k)];
            T = T + n' * n;
        end
        ev = eig(T);
        Si = sort(ev, 'descend') / numel(nx);

        if min(Si) <= Sth
            id = [id, i];
        end
    end
end

% --- Trace Descent ---
function [fp_coords, fq_coords, lSHt, lRHt, lLHt, lxHt, lyHt] = traceDescent(DEM2, FD2, S2, xd, yd, Lth, verbose)
    % Trace steepest descent paths from divide points to streams, storing coordinates.
    % Skip hilltop nodes with no valid flow paths.
    if verbose
        PB = ProgressBar(numel(xd)+1, 'taskname', 'Tracing steepest descents...', 'ui', 'cli');
        count(PB);
    end

    % Initialize fast indexing
    if ~FD2.fastindexing
        FD2.fastindexing = true;
        FD2.ixcix = zeros(FD2.size, 'uint32');
        FD2.ixcix(FD2.ix) = uint32(1):uint32(numel(FD2.ix));
    end

    % Convert divide points to indices and filter invalid ones
    ixd = coord2ind(DEM2, xd, yd);
    valid = ~isnan(ixd) & ixd > 0 & ixd <= numel(DEM2.Z);
    ixd = ixd(valid);
    xd = xd(valid);
    yd = yd(valid);
    pix = ixd;
    qix = nan(1, numel(pix));

    if isempty(pix)
        fp_coords = cell(0, 1); fq_coords = cell(0, 1);
        lSHt = []; lRHt = []; lLHt = []; lxHt = []; lyHt = [];
        return;
    end

    % Build flow path map for unique pixels
    nbr_list = ixd;
    for i = 1:numel(xd)
        [r, c] = ind2sub(size(DEM2.Z), ixd(i));
        nbr = [r-1 c; r+1 c; r c-1; r c+1; r-1 c-1; r-1 c+1; r+1 c-1; r+1 c+1];
        valid_nbr = nbr(:,1) > 0 & nbr(:,1) <= size(DEM2.Z,1) & ...
                    nbr(:,2) > 0 & nbr(:,2) <= size(DEM2.Z,2);
        nbr = nbr(valid_nbr, :);
        nbr_ind = sub2ind(size(DEM2.Z), nbr(:,1), nbr(:,2));
        nbr_ind = nbr_ind(~isnan(nbr_ind) & nbr_ind > 0);
        nbr_list = [nbr_list; nbr_ind];
    end
    unq_nbr = unique(nbr_list);
    unq_nbr = unq_nbr(~isnan(unq_nbr) & unq_nbr > 0); % Remove NaN and invalid indices

    fp = cell(numel(unq_nbr), 1);
    valid_keys = false(numel(unq_nbr), 1);
    for i = 1:numel(unq_nbr)
        pth = uint32(unq_nbr(i));
        max_iterations = 10000; % Prevent infinite loops
        iter = 0;
        while FD2.ixcix(pth(end)) ~= 0 && iter < max_iterations
            next_pixel = FD2.ixc(FD2.ixcix(pth(end)));
            if isempty(next_pixel) || isnan(next_pixel) || next_pixel == pth(end)
                break;
            end
            pth(end+1) = next_pixel;
            if any(ismember(pth, S2.IXgrid))
                break;
            end
            iter = iter + 1;
        end
        if ~isempty(pth) && ~any(isnan(pth)) && numel(pth) > 1
            fp{i} = pth(:);
            valid_keys(i) = true;
        end
    end
    % Filter valid keys and flow paths
    unq_nbr = unq_nbr(valid_keys);
    fp = fp(valid_keys);
    if isempty(unq_nbr)
        fp_coords = cell(0, 1); fq_coords = cell(0, 1);
        lSHt = []; lRHt = []; lLHt = []; lxHt = []; lyHt = [];
        return;
    end
    fp_map = containers.Map(num2cell(unq_nbr), fp);

    % Find opposing flow paths
    for i = 1:numel(xd)
        if isKey(fp_map, ixd(i))
            ixc = fp_map(ixd(i));
            [r, c] = ind2sub(size(DEM2.Z), ixd(i));
            nbr = [r-1 c; r+1 c; r c-1; r c+1; r-1 c-1; r-1 c+1; r+1 c-1; r+1 c+1];
            valid_nbr = nbr(:,1) > 0 & nbr(:,1) <= size(DEM2.Z,1) & ...
                        nbr(:,2) > 0 & nbr(:,2) <= size(DEM2.Z,2);
            nbr = nbr(valid_nbr, :);
            max_d = 0;
            for j = 1:size(nbr, 1)
                ni = sub2ind(size(DEM2.Z), nbr(j,1), nbr(j,2));
                if isKey(fp_map, ni)
                    ixc2 = fp_map(ni);
                    diffi = numel(ixc2) - sum(ismember(ixc2, ixc));
                    if diffi > max_d
                        max_d = diffi;
                        qix(i) = ixc2(1);
                    end
                end
            end
        end
    end

    % Compute flow path metrics and store coordinates
    fp_coords = cell(numel(pix), 1); fq_coords = cell(numel(qix), 1);
    lSHt = nan(numel(pix), 1); lRHt = nan(numel(pix), 1); lLHt = nan(numel(pix), 1);
    lxHt = nan(numel(pix), 1); lyHt = nan(numel(pix), 1);

    for i = 1:numel(pix)
        if ~isnan(pix(i)) && isKey(fp_map, pix(i))
            if verbose
                count(PB);
            end
            % Primary flow path coordinates
            fp_indices = fp_map(pix(i));
            [fxp, fyp] = ind2coord(DEM2, fp_indices);
            fp_coords{i} = [fxp, fyp]; % Store [x, y] coordinates
            fp_coords{i}(isnan(fp_coords{i}(:,1)) | isnan(fp_coords{i}(:,2)), :) = [];

            % Opposing flow path coordinates
            if ~isnan(qix(i)) && isKey(fp_map, qix(i))
                fq_indices = fp_map(qix(i));
                [fxq, fyq] = ind2coord(DEM2, fq_indices);
                fq_coords{i} = [fxq, fyq]; % Store [x, y] coordinates
                fq_coords{i}(isnan(fq_coords{i}(:,1)) | isnan(fq_coords{i}(:,2)), :) = [];
            else
                fq_coords{i} = nan(0, 2);
            end
            
            % Metrics for primary path
            ix = coord2ind(DEM2, fp_coords{i}(:,1), fp_coords{i}(:,2));
            ix(isnan(ix)) = [];
            zi = DEM2.Z(ix);
            if isempty(zi) || any(isnan(zi))
                fp_coords{i} = nan(0, 2);
                fq_coords{i} = nan(0, 2);
                continue;
            end
            dx = diff(fp_coords{i}(:,1)); dy = diff(fp_coords{i}(:,2));
            d_fp = sqrt(dx.^2 + dy.^2);
            slopes_fp = diff(zi) ./ d_fp;
            mean_slope_fp = median(slopes_fp, 'omitnan');
            Ri_fp = max(zi) - min(zi);
            lhi_fp = sum(d_fp);

            % Metrics for opposing path
            if size(fq_coords{i}, 1) > 0
                ix = coord2ind(DEM2, fq_coords{i}(:,1), fq_coords{i}(:,2));
                ix(isnan(ix)) = [];
                zi = DEM2.Z(ix);
                if isempty(zi) || any(isnan(zi))
                    fp_coords{i} = nan(0, 2);
                    fq_coords{i} = nan(0, 2);
                    continue;
                end
                dx = diff(fq_coords{i}(:,1)); dy = diff(fq_coords{i}(:,2));
                d_fq = sqrt(dx.^2 + dy.^2);
                slopes_fq = diff(zi) ./ d_fq;
                mean_slope_fq = median(slopes_fq, 'omitnan');
                Ri_fq = max(zi) - min(zi);
                lhi_fq = sum(d_fq);
            else
                mean_slope_fq = nan;
                Ri_fq = nan;
                lhi_fq = nan;
            end

            % Average metrics if path length exceeds threshold
            cd = [lhi_fp; lhi_fq];
            if min(cd, [], 'omitnan') > Lth
                lxHt(i) = xd(i);
                lyHt(i) = yd(i);
                lLHt(i) = mean([lhi_fp, lhi_fq], 'omitnan');
                lRHt(i) = mean([Ri_fp, Ri_fq], 'omitnan');
                lSHt(i) = mean([mean_slope_fp, mean_slope_fq], 'omitnan');
            else
                fp_coords{i} = nan(0, 2);
                fq_coords{i} = nan(0, 2);
            end
        else
            fp_coords{i} = nan(0, 2);
            fq_coords{i} = nan(0, 2);
        end
    end
end

% --- Compute Curvature ---
function OUT = computeCurvature(DEM1, ws, Gth, xd, yd, lSHt, lRHt, lLHt, lxHt, lyHt, fp_coords, fq_coords, verbose)
    % Compute hilltop curvature for valid divide points across window sizes.
    if verbose
        PB = ProgressBar(numel(xd)+1, 'taskname', 'Extracting hilltop curvature...', 'ui', 'cli');
        count(PB);
    end

    cs = DEM1.cellsize;
    lfp = fp_coords(~isnan(lxHt)); lfq = fq_coords(~isnan(lxHt));
    xd = xd(~isnan(lxHt)); yd = yd(~isnan(lxHt));
    lSHt = lSHt(~isnan(lxHt)); lRHt = lRHt(~isnan(lxHt)); lLHt = lLHt(~isnan(lxHt));
    lxHt = lxHt(~isnan(lxHt)); lyHt = lyHt(~isnan(lxHt));

    % Preallocate results
    tCht = cell(numel(xd), 1); tSht = cell(numel(xd), 1); tRht = cell(numel(xd), 1);
    tLht = cell(numel(xd), 1); txHt = cell(numel(xd), 1); tyHt = cell(numel(xd), 1);

    for i = 1:numel(xd)
        lCHt_all = nan(1, numel(ws)); lSHt_all = nan(1, numel(ws));
        lRHt_all = nan(1, numel(ws)); lLHt_all = nan(1, numel(ws));
        lxHt_all = nan(1, numel(ws)); lyHt_all = nan(1, numel(ws));

        for i2 = 1:numel(ws)
            wsi = round(ws(i2)) / cs;
            [nmat, X2, Y2, X, Y, weight] = makenormeqn(wsi * 2 + 1, cs);

            xi = xd(i); yi = yd(i);
            [xi, yi] = coord2sub(DEM1, xi, yi);
            try
                z = double(DEM1.Z(xi-wsi:xi+wsi, yi-wsi:yi+wsi));
            catch
                continue;
            end

            % Compute normal equations
            zvec = [sum(sum(z .* weight .* X2)); sum(sum(z .* weight .* Y2)); ...
                    sum(sum(z .* weight .* X .* Y)); sum(sum(z .* weight .* X)); ...
                    sum(sum(z .* weight .* Y)); sum(sum(z .* weight))];
            A = nmat \ zvec;

            Gi = sqrt(A(4)^2 + A(5)^2);
            Ci = 2 * (A(1) + A(2));

            if Ci <= 0 && Gi <= Gth
                lCHt_all(i2) = Ci;
                lSHt_all(i2) = lSHt(i);
                lRHt_all(i2) = lRHt(i);
                lLHt_all(i2) = lLHt(i);
                lxHt_all(i2) = xd(i);
                lyHt_all(i2) = yd(i);
            end
        end

        tCht{i} = lCHt_all; tSht{i} = lSHt_all; tRht{i} = lRHt_all;
        tLht{i} = lLHt_all; txHt{i} = lxHt_all; tyHt{i} = lyHt_all;

        if verbose
            count(PB);
        end
    end

    % Assemble output
    OUT.Cht = vertcat(tCht{:}); OUT.Sht = vertcat(tSht{:}); OUT.Rht = vertcat(tRht{:});
    OUT.Lht = vertcat(tLht{:}); OUT.xht = vertcat(txHt{:}); OUT.yht = vertcat(tyHt{:});
    OUT.ws = ws;
    OUT.fpq = lfq; OUT.fpp = lfp;
end

% --- Make Normal Equations ---
function [nmat, X2, Y2, X, Y, weight] = makenormeqn(win, g)
    % Compute normal equations for fitting a 2D quadratic surface to DEM data.
    % Adapted from Struble and Roering (2021), Earth Surf. Dynam.
    if mod(win, 2) == 0 || win < 1
        error('Window size must be a positive odd integer');
    end
    if g <= 0
        error('Grid spacing must be positive');
    end

    edge = g * (win - 1) / 2;
    x = (-edge:g:edge)';
    [X, Y] = meshgrid(x, x);
    Y = flipud(Y);
    X2 = X.^2; Y2 = Y.^2;

    distance = sqrt(X.^2 + Y.^2);
    dcrit = sqrt(edge^2 + edge^2);
    weight = 1 ./ ((distance ./ dcrit) + 1).^(3/2);

    X3 = X.^3; Y3 = Y.^3;
    X4 = X.^4; Y4 = Y.^4;
    x4w = sum(sum(X4 .* weight)); y4w = sum(sum(Y4 .* weight));
    x2y2w = sum(sum(X2 .* Y2 .* weight)); x2w = sum(sum(X2 .* weight));
    y2w = sum(sum(Y2 .* weight)); ww = sum(sum(weight));
    x3yw = sum(sum(X3 .* Y .* weight)); x3w = sum(sum(X3 .* weight));
    x2yw = sum(sum(X2 .* Y .* weight)); xy3w = sum(sum(X .* Y3 .* weight));
    xy2w = sum(sum(X .* Y2 .* weight)); y3w = sum(sum(Y3 .* weight));
    xyw = sum(sum(X .* Y .* weight)); xw = sum(sum(X .* weight));
    yw = sum(sum(Y .* weight));

    nmat = [x4w x2y2w x3yw x3w x2yw x2w; ...
            x2y2w y4w xy3w xy2w y3w y2w; ...
            x3yw xy3w x2y2w x2yw xy2w xyw; ...
            x3w xy2w x2yw x2w xyw xw; ...
            x2yw y3w xy2w xyw y2w yw; ...
            x2w y2w xyw xw yw ww];
end