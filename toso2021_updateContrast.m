contrasts = eval(contrast_str);
contrast_set = eval([contrast_str(1:end-1),'_set']);
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = eval([contrast_str,'_clrs']);
contrast_units = eval([contrast_str(1),'_units']);
contrast_lbl = [upper(contrast_str(1)),'_',contrast_str(2)];

