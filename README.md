# toso2021-comment

Matlab code (tested on versions 2019b and 2020b) implementing the reanalyzes of the behavioral & neural [data](https://data.mendeley.com/datasets/wp9h39kbtv/2) originally published in [Toso et al. 2021](https://doi.org/10.1016/j.neuron.2021.08.020).

### toso2021_main.m  
- Loads the data;
- Selects which task variant to analyze (delayed duration or intensity comparison);
- Sets _if_ and _where_ to save figures;
- Curates the data & prints _before_ & _after_ metrics;
- Parses the data;
- Sets neuron selection criteria;
- Sets aesthetic preferences for figures & axes;
- Sets all color schemes;
- Sets which experimental variable to use as contrast (e.g., I2);
- Runs the scripts corresponding to figures 1-7 (and S1-2) in sequence;
