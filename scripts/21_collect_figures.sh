#!/bin/bash
# Copy the npj Parkinson's Disease figure set into figures/npjPD_final with final numbering.
cd "/path/to/project/figures" || exit 1
mkdir -p npjPD_final
rm -f npjPD_final/SuppFig*          # supplementary figures were renumbered by first citation (subgroups = 1, frequency = 2)
for ext in pdf png tif; do
  cp FigS1_flowchart.$ext               npjPD_final/Fig1_study_design.$ext
  cp Fig2_phenome_wide.$ext             npjPD_final/Fig2_phenome_wide.$ext
  cp Fig3_pd_transit_spectrum.$ext      npjPD_final/Fig3_dose_response_transit_spectrum.$ext
  cp Fig4_pd_robustness.$ext            npjPD_final/Fig4_robustness.$ext
  cp Fig5_pd_prediction.$ext            npjPD_final/Fig5_prediction_prodromal.$ext
  cp FigS2_pd_subgroups.$ext            npjPD_final/SuppFig1_subgroups.$ext
  cp Fig4_frequency_cardiometabolic.$ext npjPD_final/SuppFig2_frequency_other_outcomes.$ext
  cp FigS3_genetics_imaging.$ext        npjPD_final/SuppFig3_genetics_imaging.$ext
done
ls npjPD_final
