# Spectral Index Analysis

This repository contains the code and results of the spectral index analysis of the LMC region at 888 MHz, 1367 MHz, and 1419 MHz frequencies. The steps of the analysis are referred to Pennock et al (2021), as there are only three frequencies available for 30 Dor B in ASKAP public archives. The flux errors are estimated based on Prof. Sheng-Yuan's advice during the 2022 ISM class.

## Data Preparation

The 888 MHz, 1367 MHz, and 1419 MHz images were convolved to a common resolution of 18 arcsec by 18 arcsec to ensure consistent comparisons. The flux densities were extracted from these images, and a constraint was applied to select points with a signal-to-noise ratio (SNR) greater than or equal to 5 in both the 888 MHz and 1367 MHz images.

## Flux Calculation and Error Estimation

To estimate the background noise level, an empty region was selected, and the flux errors were calculated based on neighboring pixels within a distance of 2. This approach accounts for spatial variations in the noise and provides a more accurate estimation of the uncertainties in the flux measurements.

## Spectral Index Calculation

The spectral index was calculated as the logarithmic ratio of flux densities at different frequencies, following the equation: $$F_\nu \propto \nu^\alpha$$


## Histogram Analysis

A histogram was generated to visualize the distribution of the spectral index values. To identify statistically significant spectral indices, a rejection criterion based on a probability threshold of 0.05 is applied. Spectral indices failing to meet the criterion were masked, indicating their potential lack of significance. 

![image](https://github.com/boan-chen/30DorB-ASKAP/assets/108161781/a5b4b749-cbae-42d0-8fb3-7ab7688f79ea)

## Results

The analysis revealed two distinct components in the spectral index map. Region 1 represents the main structure and shows a steep negative spectral index, indicative of synchrotron emission. On the other hand, region 2 exhibits a filament-like structure with a positive spectral index, suggesting the presence of thermal radiated structures associated with free-free emission processes.

The provided figures illustrate the spectral index map for region 1 and region 2, highlighting their different characteristics. Additionally, the histogram of the spectral index distribution is included, showcasing the statistical significance of the spectral indices. The median value of the spectral index is -0.7, which is consistent with the result of Pennock et al (2021).

![image](https://github.com/boan-chen/30DorB-ASKAP/assets/108161781/5d486601-ca22-48a1-a071-2bf3739b88b4)

## Usage

The code and data files necessary to reproduce the analysis are provided in the repository. Detailed instructions on running the code are documented in the comments of pipline.py.

## References

Please acknowledge the following reference if you use this spectral index analysis in your research:

- Clara M Pennock and others, The ASKAP-EMU Early Science Project: 888 MHz radio continuum survey of the Large Magellanic Cloud, Monthly Notices of the Royal Astronomical Society, Volume 506, Issue 3, September 2021, Pages 3540–3559, https://doi.org/10.1093/mnras/stab1858

Feel free to explore the repository and utilize the results for further research and analysis.

## Note for me
Images are to be binned to 2*2 pixels, or 36*36 arcsec. Find Chung for convolving ATCA 5GHz (6cm telescope) data. After the steps, send Wang an email for advice.
