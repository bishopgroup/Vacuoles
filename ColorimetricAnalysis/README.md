1. RGB values of capillary images.ipynb
Extracts RGB values from microscopy images of capillaries containing solutions with unknown and known pH values. This serves as the raw input for the pH calibration process.
Input: Capillary images
Output: RGB values
2. pH calibration of capillary images.ipynb
Converts RGB values obtained from capillary images into the CIE 1931 xy chromaticity space. Uses the xy values of capillary images with known pH values to create a calibration curve that links color to pH. Uses calibration curve to estimate pH of capillary images with unknown pH values.
Input: RGB values of capillary images of known and known pH values.
Output: pH values of capillary mages with known pH values.
3. Drop detection, radius and angular variance.ipynb
Performs drop detection from microscopy images. For each drop, it calculates geometric features such as radius and angular variance, which can be used to identify vacuoles or deformations.
Input: Drop images
Output: Drop coordinates, radius, and angular variance, first point of vacuole formation
4. pH Mapping of drops by location and time.ipynb
Maps pH values to individual detected drops by matching each drop’s spatial position (x) and time label to the calibrated capillary data.
Input: Drop positions from (3) and capillary pH data from (2)
Output: Estimated pH value for each drop
