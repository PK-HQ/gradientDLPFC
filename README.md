# gradientDLPFC

MATLAB analysis code for Tan et al. (2023), *J. Neurosci.* — functional gradients 
and anterior-posterior parcellation of primate lateral prefrontal cortex (LPFC) 
during a visual working memory task.

## Usage
```
1. Select desired analyses in analyzeSpikeData.m
2. Run analyzeSpikeData in MATLAB
```

## Key scripts
- `analyzeSpikeData.m` — main entry point; calls gradient and decoding analyses
- `fitGradients.m` — fits functional gradients along the A-P axis
- `getFunctionalParc.m` — derives functional parcellation boundaries
- `calcFuncMeasures.m` — computes selectivity, latency, and other neural measures

## Publication
Tan PK, Tang C, Herikstad R, Pillay A, Libedinsky C (2023). Distinct lateral 
prefrontal regions are organized in an anterior-posterior functional gradient. 
*J. Neurosci.* [https://doi.org/10.1523/JNEUROSCI.XXXX]

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
