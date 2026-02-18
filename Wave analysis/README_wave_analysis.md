# Bioluminescence-wave-analysis-trackmate
This script allows the analysis of traces obtained using the module Trackmate (ImageJ)<br/>
Python version  - 3.8.10<br/>
<br/>
The script uses .csv files derived from the ImageJ module Trackmate as input. The excel files (https://pubmed.ncbi.nlm.nih.gov/27713081/)<br/>
<br/>
Traces are selected based on length, filtered, smoothed and detrended to detect extrema. Then, period is calculated in non-detrended traces as the difference between consecutive half-maximum points. <br/>
