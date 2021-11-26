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
- Plots stimulus pairs with the corresponding average performance;
![Alt text](panels/sampling_scheme.svg?raw=true)
- Same as before, plus a gradient with hypothesized continuous performance so as to allow for a visualization of *contraction bias* on T1.
![Alt text](panels/contraction_bias.svg?raw=true)
![Alt text](panels/psychometric_curves_i1.svg?raw=true)
![Alt text](panels/psychometric_curves_i2.svg?raw=true)
- Fits a generalized linear model (GLM) to choice data using T1, T2, I1 & I2 as predictors;
![Alt text](panels/choice_GLM.svg?raw=true)

### toso2021_neuronSelection.m

### toso2021_trialTypeDistributions.m

### toso2021_overallModulation.m

### toso2021_rasters.m

### toso2021_PCA.m

### toso2021_neurometricCurves.m

### toso2021_naiveBayesDecoder.m
