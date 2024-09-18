
# Reanalysis of Delavaux et al., "Mutualisms weaken the latitudinal diversity gradient among oceanic islands"

Changes in the original Delavaux scripts include:

- Myc_analyses.R: Save data that was used to produce the world map (Figure 1b in the original manuscript)
- Nfix_Analyses.R, Myc_Analyses.R, and Poll_Analyses.R: Add Random Forest estimator for expected mutualism ratios
- Nfix_Analyses.R, Myc_Analyses.R, and Poll_Analyses.R: Return observed mutualism ratio on islands
- Joint_Analyses.R: Calculate observed mutualism ratio on islands (new variable called biotic.ml_obs)
- Joint_Analyses.R: Calculate predicted/corrected mutualism ratio on islands (new variable called biotic.ml_rf)
- Joint_Analyses.R: Save data that was used to calculate the effects (e.g. Figure 2c in the original manuscript)

First, rerun all scripts from the original release to generate and save all necessary data.

Reanalysis can be found in `2_Analyses/Reanalysis.R` 