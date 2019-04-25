# Endodemog_2019<br/>
Population demography of grass fungal endophyte populations<br/>
Authors: Joshua C. Fowler and Tom E.X. Miller<br/>
ReadMe Update: 4/25/2019<br/>

Repo info:<br/>
This repo contains files for vital rate models and figure creation. Currently the vital rate models are in progress and these will be incorporated into an IPM. All stan models are written within the R scripts.

File info:<br/>
endodemog_surv_mixed_stan.R - survival vital rate model and start of ppc's<br/>
endodemog_grow_mixed_stan.R - growth vital rate model and start of ppc's<br/>
endodemog_flw_mixed_stan.R - flowering status vital rate model and start of ppc's<br/>
endodemog_fert_mixed_stan.R - fertility vital rate model. currently looking at prior on phi dispersion parameter which helped the model to converge. need to add ppc's<br/>
endodemog_seed_means_stan.R - calculate seed/spikelet means for all species can be used with fertility model to estimate seeds/plant
endodemog_seed_to_seedling_stan.R - recruitment vital rate model, currently run with Jenn's seed values and need to update to use new calculated seed values. Dealing with binomial model and plots where seeds<recruits.<br/>
endodemog_data_processing.R - data prep script to take in endo_demog_long.csv file, perform data processing for each vital rate model, and also to extract raw repro data from raw data excel files and merge this for the fertility models. Currently there is still data processing for the vital rate models in their respective scripts, but I have incorporated the surv and grw data processing into this script. currently working on flw_t1 records in endo_demog_long to generate flw_t for 2016 adn 2017.<br/>
endodemog_figure.R - combined figures script for surv/grw/flw vital rate models so far. Plan to transition this to an rMarkdown file and incorporate fertility models. (combines endodemog_surv_visualization.R, endodemog_grow_visualization.R, endodemog_flw_visualization.R)<br/>
Tom_figures.R - Tom file to redraw some figures<br/>
MeetingUpdates.tex - weekly meeting notes<br/>

Old/practice files:
endodemog_surv_visualization.R<br/>
endodemog_grow_visualization.R<br/>
endodemog_flw_visualization.R<br/>
endodemog_surv_full.stan<br/>
endodemog_surv_matrix.stan<br/>
endodemog_flw_full.stan<br/>
endodemog_surv.stan<br/>

