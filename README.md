# gradientDLPFC

I independently authored this MATLAB analysis code for our study of spatial working memory in macaque dorsolateral prefrontal cortex (dlPFC). The repository is an academic research archive. The published paper describes the experiments, analysis decisions and findings.

## Research question

Does function change smoothly along the anterior–posterior axis of dlPFC, or does the overall gradient arise from areas with distinct functional properties? We recorded single neurons while two macaques performed a memory-guided saccade task with an intervening distractor. We compared a smooth, single-segment gradient with piecewise models that allowed a boundary between areas.

Single-neuron functional measures favored an areal gradient. At the population level, target information during the delays was distributed and generally stronger posteriorly. The population analysis could not distinguish a smooth gradient from an areal gradient at that resolution. Decoders generalized better within a delay than between delays on opposite sides of the distractor, suggesting that the representation changed even where target location remained decodable. [Published results and full text](https://dspace.mit.edu/handle/1721.1/168246)

## MATLAB analyses

| Analysis | Implementation |
| --- | --- |
| Functional measures | Estimate response latency, target and distractor selectivity, distractor filtering, and delay-period measures. [calcFuncMeasures.m](calcFuncMeasures.m) |
| Gradient comparisons | The published comparison used principal components of single-neuron measures. The archived plotting workflow also fits smoothed individual measures with robust single-segment and segmented regressions, with shuffled or surrogate comparisons. Additional segment-count and HPC workflows are included. [plotGradientFits.m](plotGradientFits.m) · [fitGradients.m](fitGradients.m) · [fitGradientsHPC.m](fitGradientsHPC.m) |
| Functional parcels | Use estimated boundaries and electrode mappings to assign neurons to groups for population analysis. [getFunctionalParc.m](getFunctionalParc.m) |
| Population and cross-temporal decoding | Subsample to matched neuron counts across groups and decode location with linear discriminant analysis. The decoding code also includes PCA-based denoising. Evaluate on matching and different time bins to compare within-period decodability with cross-period generalization. [runCTDAnalysis.m](runCTDAnalysis.m) · [runCTD4.m](runCTD4.m) · [runCTD.m](runCTD.m) |

## Entry point and inputs

[analyzeSpikeData.m](analyzeSpikeData.m) is the original entry point. Add the repository’s MATLAB files to the MATLAB path, set `dataPath` to the original analysis root or its local equivalent, and select the analysis with the script’s `execute...` flags before running `analyzeSpikeData`. The checked-in settings enable gradient analysis and disable cross-temporal decoding. The decoding section defaults to preloading saved HPC results when enabled.

The driver changes into the analysis root and expects paths relative to it. The following are examples of the inputs required by the corresponding stages.

```text
<analysis-root>/
  data/
    raw/datasetsFull632_8R.mat
    neuronidx/                 # Functional measures, electrode maps and neuron indices
    em/                        # Gradient fits and intermediate results
    ctd/                       # Saved cross-temporal decoding results
  plots/                       # Output subdirectories specified in the driver
```

The driver contains an absolute macOS path. Some helpers also embed original workstation or HPC-era paths, so changing `dataPath` alone may leave unresolved file references. MATLAB’s Statistics and Machine Learning Toolbox is used for routines such as `fitlm` and `classify`. Different analysis switches can require different saved intermediates and supporting helpers.

The public checkout provides the MATLAB source for inspection. Executing the analyses requires the original recordings or the corresponding saved intermediates, which are absent from the checkout. The paper’s Data Availability statement directs data requests to the corresponding author.

## Publication

Tan PK, Tang C, Herikstad R, Pillay A, Libedinsky C. (2023). Distinct lateral prefrontal regions are organized in an anterior–posterior functional gradient. *Journal of Neuroscience*, 43(38), 6564–6572. [DOI 10.1523/JNEUROSCI.0007-23.2023](https://doi.org/10.1523/JNEUROSCI.0007-23.2023)

[Existing repository license notice · CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)
