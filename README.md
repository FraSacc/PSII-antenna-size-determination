# PSII-antenna-size-determination
Determines plant Photosystem II (PSII) functional antenna cross-section using the Malkin method. The starting plant material must be treated with DCMU to block reaction centres.

Input: CSV file with the following columns:
- Time, s
- Time, ms
- raw (1 per replicate)
- offset (1 per replicate, with Fo set at the origin)

Output: 
- integral: dictionary for each sample containing the area above the curve for all replicates
- percentages: dictionary for each sample containing the average percentage (WT set to 100) and associated error
- Raw_full.svg: Graph with full fluorescence traces from PAM data
- Raw&averages.svg: Graph with normalised and offset traces. Panel1 raw, panel2 averages+- st err
- AreaAbove.svg: Example of the method
- PieCharts.svg: Pie charts with WT=100%
    
References:
- Malkin method - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8940271/
- Error propagation - https://www.geol.lsu.edu/jlorenzo/geophysics/uncertainties/Uncertaintiespart2.html
