==============================================================
NeuroPower: Compute statistical power and sample size for fMRI
==============================================================

This power calculation tool was first introduced in this `OHBM poster <http://users.ugent.be/~jdurnez/ProsPeakPow_1406_OHBM.pdf>`_ .

This app applies functions to compute post-hoc power in a fMRI pilot study as described in `Durnez, Moerkerke & Nichols, 2014 <http://www.ncbi.nlm.nih.gov/pubmed/23927901>`_ and to compute the minimal sample size needed to obtain a preferred power rate in a subsequent fMRI study based on that pilot study. 

Use
---

In a first step, you can upload your statistical map (z-map or t-map) to extract peak information. Note that the data are NOT stored on any server once the peak information is extracted. This step usually takes less than 1 minute. The result of this step can be seen in the tab Peaks , where a table will be displayed with all local maxima above the preferred excursion threshold. The displayed x-,y- and z-coordinates are in voxel space.

Subsequently, the mixture model of inactive and active peaks is fit to the data. This model fitting includes an estimation of both the prevalence of active peaks, the mean (delta) and the standard deviation (sigma) of the distribution of active peaks. The fit of the model can be inpected in the tab Model fit . The left hand panel will show a histogram of the peak p-values with their estimated distribution, and the right hand panel will show a histogram of the peak heights (t or z) with their estimated null, alternative and combined distributions.

Finally, the power can be calculated. When entering the preferred thresholding procedure and level of significance, the post-hoc power will be estimated for the current data (result displayed in tab Post-hoc ), as well as the minimal sample size for a certain level of power in a new dataset (result displayed in tab Power ). 

Support and communication
-------------------------
If you have a problem or questions, please e-mail me at joke.durnez@gmail.com

