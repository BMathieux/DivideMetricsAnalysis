function HMetrics = CHt_metrics(BData, DEM1, params, varargin)
    % This function computes hilltop metrics for each basin defined in BData
    % using high-resolution and resampled DEMs. Curvature calculations are performed
    % only if runCHt is true. Parallelization is handled internally by the extract_CHt
    % function when useParallel is true.
    %
    % Inputs:
    %   BData    - Structure containing basin characteristics and divide data.
    %   DEM1     - High-resolution digital elevation model (GRIDobj).
    %   params   - Analysis parameters defined by the DefineParams function.
    %   varargin - Optional settings: 'Verbose' (true), 'UseParallel' (false).
    %
    % Outputs:
    %   HMetrics - Cell array of hillslope metrics (NaN if runCHt is false).

    % --- Parse Inputs ---
    p = inputParser;
    addParameter(p, 'Verbose', true, @islogical);
    addParameter(p, 'UseParallel', false, @islogical);
    parse(p, varargin{:});
    verbose = p.Results.Verbose;
    useParallel = p.Results.UseParallel;

    DEM2 = BData.DEM2;
    xd = BData.xdi;
    yd = BData.ydi;
    FD2 = FLOWobj(DEM2, 'preprocess', 'carve');

    HMetrics = cell(numel(BData.DEMi), 1);
    if verbose
        PB = ProgressBar(numel(BData.DEMi)+1, 'taskname', 'Extract hillslope metrics...', 'ui', 'cli');
        count(PB);
    end

    for i1 = 1:numel(BData.DEMi)
        D = dependencemap(FD2, BData.ixB(i1));
        mask = D.Z > 0;
        buf = imdilate(mask, strel('disk', 1));
        DEMi = DEM2;
        DEMi.Z(~buf) = NaN;
        DEMi = crop(DEMi, buf, NaN);
        [xbi, ybi] = getoutline(DEMi, true);
        valid = ~isnan(xbi) & ~isnan(ybi) & ~isinf(xbi) & ~isinf(ybi);
        xbi = xbi(valid); ybi = ybi(valid);
        in_box = xd >= min(xbi) & xd <= max(xbi) & yd >= min(ybi) & yd <= max(ybi);
        inP = false(size(xd));
        inP(in_box) = inpolygon(xd(in_box), yd(in_box), xbi, ybi);
        xdi = xd(inP);
        ydi = yd(inP);

        buf = imdilate(mask, strel('disk', max(params.ws)/DEM1.cellsize+1));
        DEMi = DEM2;
        DEMi.Z(~buf) = NaN;
        DEMi = crop(DEMi, buf, NaN);
        [xbi, ybi] = getoutline(DEMi, true);
        DEM1i = crop(DEM1, [min(xbi) max(xbi)], [min(ybi) max(ybi)]);

        DEM2i = resample(DEM1i, DEM2.cellsize);
        FD2i = FLOWobj(DEM2i, 'preprocess', 'carve');
        S2i = STREAMobj(FD2i, flowacc(FD2i).*DEM2.cellsize^2 > params.A0);

        % Prepare parameters for curvature calculation
        iparams = struct();
        iparams.DEM1 = DEM1i;
        iparams.DEM2 = DEM2i;
        iparams.FD2 = FD2i;
        iparams.S2 = S2i;
        iparams.xd = xdi;
        iparams.yd = ydi;
        iparams.wsR = params.wsR;
        iparams.Sth = params.Sth;
        iparams.Gth = params.Gth;
        iparams.Lth = params.Lth;
        iparams.ws = params.ws;

        HMetrics{i1} = extract_CHt(iparams, useParallel, false);

        if verbose
            count(PB);
        end
    end
end
