1. Drop detection, radius and angular variance.ipynb

Performs drop detection from microscopy images. For each drop, it calculates geometric features such as radius and angular variance, which can be used to identify vacuoles or deformations.

Input: Drop images

Output: Drop coordinates, radius, and angular variance, first point of vacuole formation

2. pH Mapping of drops by location and time.ipynb

Maps pH values to individual detected drops by matching each drop’s spatial position (x) and time label to the calibrated capillary data.

Input: Drop positions from (3) and capillary pH data from (2)

Output: Estimated pH value for each drop
