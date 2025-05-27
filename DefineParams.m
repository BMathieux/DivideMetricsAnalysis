function params = DefineParams(varargin)
    % Define parameters for divide metrics analysis.
    % This function sets default parameters for use in extractBasins, CHt_metrics, 
    % SelectWs, and CalcDivideTable functions, with optional overrides provided 
    % as name-value pairs. Parameters are grouped by the functions that call them.
    % Parameters Amin and wl can be chosen interactively: Amin in extractBasins 
    % and wl in SelectWs. These parameters are set to defaults if interactive options
    % are not chosen by the user.
    % The parameter nb_ite is used for interactive selection of Amin in extractBasins.
    %
    % Inputs:
    %   varargin - Optional name-value pairs to override default parameters.
    %
    % Outputs:
    %   params   - Structure containing analysis parameters.
    %
    % Parameter Descriptions:
    %
    % --- Parameters for extractBasins (basin delineation and divide extraction)
    %   A0       - Critical drainage area for stream networks [m²].
    %   res      - Resampled DEM resolution for basin extraction [m].
    %   Amin     - Default drainage area for basin extraction if not selected interactively [m², greater than A0].
    %   Dmin     - Minimum divide length (used in function DIVIDEobj/shrink) [m].
    %   nb_ite   - Number of candidate Amin thresholds for interactive selection of Amin [integer].
    %
    % --- Parameters for CHt_metrics (hilltop metric calculations)
    %   A0       - Critical drainage area for stream networks [m², shared with extractBasins and CalcDivideTable].
    %   ws       - Window size radius for curvature extraction [m].
    %   K        - Soil transport coefficient [m²/yr, shared with CalcDivideTable].
    %   dK       - Soil transport coefficient uncertainty [m²/yr].
    %   rp       - Soil-to-rock density ratio [dimensionless, shared with CalcDivideTable].
    %   nbstrp   - Number of bootstrap samples for uncertainty [integer, shared with CalcDivideTable].
    %   wsR      - Window size radius for roughness filtering [pixels].
    %   Sth      - Minimum eigenvalue threshold for roughness filtering [dimensionless].
    %   Gth      - Gradient threshold for curvature extraction [m/m].
    %   Lth      - Minimum flow path length for steepest descent traces [meters].
    %
    % --- Parameters for SelectWs (interactive window size selection)
    %   wl       - Default window size for curvature measurements if not selected interactively [m, within ws range].
    %
    % --- Parameters for CalcDivideTable (divide metrics for paired basins)
    %   A0       - Critical drainage area for stream networks [m², shared with extractBasins and CHt_metrics].
    %   mn       - m/n ratio for chi calculation [dimensionless, 0.3–0.6].
    %   zb       - Base levels for chi calculation [m].
    %   wH       - Local relief for Gilbert Metrics [m].
    %   nbstrp   - Number of bootstrap samples for uncertainty [integer, shared with CHt_metrics].
    %   K        - Soil transport coefficient [m²/yr, shared with CHt_metrics].
    %   rp       - Soil-to-rock density ratio [dimensionless, shared with CHt_metrics].

    %% --- Define Parameters ---

    % --- Parameters for extractBasins ---
    A0 = 9e4;       % Critical drainage area (m²)
    res = 30;       % Resampled DEM resolution (m)
    Amin = 1.75e6;   % Reference drainage area for basins (m²), default if not selected interactively
    Dmin = 250;     % Minimum divide length (m)
    nb_ite = 200;    % Number of candidate Amin thresholds for interactive selection of Amin in extractBasins

    % --- Parameters for CHt_metrics ---
    ws = [1:15, 16:4:40]; % Window sizes for curvature (m)
    K = 1.56e-2;          % Soil transport coefficient (m²/yr)
    dK = 1.01e-3;         % Soil transport coefficient uncertainty (m²/yr)
    rp = 0.5;             % Soil-to-rock density ratio
    nbstrp = 1e3;         % Number of bootstrap samples
    wsR = 2;              % Window size radius for roughness filtering (pixels)
    Sth = 0.015;          % Eigenvalue threshold for roughness filtering
    Gth = 0.4;            % Gradient threshold for curvature (m/m)
    Lth = 50;             % Minimum flow path length (m)

    % --- Parameters for SelectWs ---
    wl = 11;              % Default window length radius for hilltop curvature (m), default if not selected interactively

    % --- Parameters for CalcDivideTable ---
    mn = 0.45;            % m/n ratio for chi calculation
    zb = [300;450;600]; % Base levels for chi (m)
    wH = 500;             % Local relief for Gilbert Metrics (m)

    %% --- Parse Inputs ---
    p = inputParser;
    % extractBasins parameters
    addParameter(p, 'A0', A0, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'res', res, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'Amin', Amin, @(x) isnumeric(x) && x > A0);
    addParameter(p, 'Dmin', Dmin, @(x) isnumeric(x) && x > 100);
    addParameter(p, 'nb_ite', nb_ite, @(x) isnumeric(x) && x > 2);
    % CHt_metrics parameters
    addParameter(p, 'ws', ws, @(x) isnumeric(x) && all(x > 0));
    addParameter(p, 'K', K, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'dK', dK, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'rp', rp, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'nbstrp', nbstrp, @(x) isnumeric(x) && x > 0 && mod(x, 1) == 0);
    addParameter(p, 'wsR', wsR, @(x) isnumeric(x) && x > 0 && mod(x, 1) == 0);
    addParameter(p, 'Sth', Sth, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'Gth', Gth, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'Lth', Lth, @(x) isnumeric(x) && x > 0);
    % SelectWs parameters
    addParameter(p, 'wl', wl, @(x) isnumeric(x) && x > 0);
    % CalcDivideTable parameters
    addParameter(p, 'mn', mn, @(x) isnumeric(x) && x > 0 && x < 1);
    addParameter(p, 'zb', zb, @(x) isnumeric(x) && all(x > 0));
    addParameter(p, 'wH', wH, @(x) isnumeric(x) && x > 0);
    parse(p, varargin{:});
    params = p.Results;
end