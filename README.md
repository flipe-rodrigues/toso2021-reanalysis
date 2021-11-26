# toso2021-comment

Matlab code (tested on versions 2019b and 2020b) for reanalyzing the behavioral & neural [data](https://data.mendeley.com/datasets/wp9h39kbtv/2) from [Toso et al. 2021](https://doi.org/10.1016/j.neuron.2021.08.020).

### toso2021_wrapper.m  
- Loads the data;
- Sets *if* and *where* to save figures;
- Runs all other scripts in sequence (in the same order as they appear below);

### toso2021_preface.m
- Curates & parses the data;
- Sets aesthetic preferences for figures & axes;
- Sets all color schemes;

### toso2021_behavior.m
![Alt text](panels/sampling_scheme.svg?raw=true "psycurves")
![Alt text](panels/psychometric_curves_T2-T1.svg?raw=true "psycurves")
![Alt text](panels/contraction_bias.svg?raw=true "psycurves")
![Alt text](panels/psychometric_curves_T2.svg?raw=true "psycurves")
![Alt text](panels/choice_GLM.svg?raw=true "choiceGLM")

### toso2021_neuronSelection.m

### toso2021_trialTypeDistributions.m

### toso2021_overallModulation.m

### toso2021_rasters.m

### toso2021_PCA.m

### toso2021_neurometricCurves.m

### toso2021_naiveBayesDecoder.m
