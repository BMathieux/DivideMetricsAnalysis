# DivideMetricsAnalysis

Repository associated with "Drainage divide dynamics in a low-deformation setting: Insights from morphometry, hilltop curvature and cosmogenic nuclides" submitted to EPSL. 

## Overview 

The Divide Metrics Analysis toolkit processes DEMs to extract and analyze metrics for basins on either side of a drainage divide. It supports both automated and interactive workflows, allowing users to delineate basins, compute topographic and fluvial metrics, and select thresholds (e.g., reference drainage area for basin extraction) and divide via graphical interfaces. The toolkit is built around a modular set of MATLAB functions, leveraging the TopoToolbox library for robust geospatial processing of DEMs, flow routing, and divide network analysis.

## Structure

The main folder contains the entry point script `RunDivideAnalysis.m` to execute the entire workflow and `DefineParams.m` for parameters definitions. `Helpers` includes utility functions for progress display and user-assisted thresholds selection with graphical interfaces (`ProgressBar.m`,`SelectDivide.m`, `SelectWs.m`). `Core` houses the core analysis functions (`Cht_metrics.m`, `CalcDivideTable.m`, `extractBasins.m`) that perform the main computations.

## Installation

1. Clone the repository: `git clone <repository-url>` (replace `<repository-url>` with your repository URL) or download the ZIP file from GitHub and extract it.
2. In MATLAB, navigate to the `DivideMetricsAnalysis-main` folder and run:
   ```matlab
   addpath(genpath('DivideMetricsAnalysis-main'));
3. Ensure TopoToolbox v2 (https://github.com/TopoToolbox/topotoolbox) is installed and added to your MATLAB path. MATLAB R2020a or later is recommended (tested on R2022b). Statistics and Machine Learning Toolbox is required. The Parallel Computing Toolbox is optional for parallel processing.
   
Note: For now, the toolkit is not compatible with Topotoolbox v3.

## How to run the Code

1. Update parameters in DefineParams.m (see the H1 comment for details).
2. Follow the prompts for updating DEMs file paths, selecting interactively or not the reference drainage area for basin extraction and run or not hilltop curvature measurements with the extrapolation of associated metrics (i.e., denudation rates and divide migration rates) 
3. Navigate to the Main folder in MATLAB and run `RunDivideAnalysis.m`
4. The results will be saved as csv files.

## Contact

For questions, contact b.mathieux@unistra.fr or open an issue on GitHub. 

Bastien Mathieux - University of Strasbourg
