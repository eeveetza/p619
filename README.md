# MATLAB/Octave Implementation of Recommendation ITU-R P.619

<!--[![DOI](https://zenodo.org/badge/459560677.svg)](https://zenodo.org/badge/latestdoi/459560677) -->

This development code repository contains a MATLAB/Octave software implementation of [Recommendation ITU-R P.619-5](https://www.itu.int/rec/R-REC-P.619/en) "Propagation data required for the evaluation of interference between stations in space and those on the surface of the Earth".  

This is a very first implementation. It is not yet complete and is work under progress. 
It has not been tested and may contain errors and bugs.

The following table describes the structure of the folder `./matlab/` containing the MATLAB/Octave implementation of Recommendation ITU-R P.619.

<!--Methods using digital data maps need to be optimized. They run relatively slow on MATLAB/Octave on MacOS. They seem to be running OK on MATLAB on MS Windows OS.-->



| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`P619.m`                | MATLAB class with methods implementing Recommendation ITU-R P.619-5        |
|`initiate_digital_maps.m`| MATLAB script that processes the ITU-R maps and generates the necessary functions. It needs to be run prior to using this software implementation. For details, see [Integrating ITU Digital Products](#integrating-itu-digital-products). |
|`test_p619.m`          | MATLAB script (under development) that will be used to validate the implementation of this Recommendation  |


## Integrating ITU Digital Products

This software uses ITU digital products that are integral part of Recommendations. These products must not be reproduced or distributed without explicit written permission from the ITU.

### Setup Instructions

1. **Download and extract the required maps** to `./private/maps`:

   - From [ITU-R P.452-18](https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.452-18-202310-I!!ZIP-E.zip):
     - `DN50.TXT`

2. **Run the script** `initiate_digital_maps.m` to generate the necessary functions for retrieving and interpolating data from from the maps.

### Notes

- Ensure all files are placed in `./private/maps` before running the script.
- The script processes the maps, which are critical for the software’s functionality.
- The resulting `*.m` files are placed in the folder `./private`.

## Methods implemented in class P5P619
| Function          | Reference  | Description  |
|-------------------|------------|--------------|
|``| §   |  |



## Function Call


## Required input arguments

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|

## Optional input arguments
| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|


## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|



## Software Versions
The code was tested and runs on:
* MATLAB version 2022a 
* Octave version 6.1.0

## References

* [Recommendation ITU-R P.619](https://www.itu.int/rec/R-REC-P.619/en)
* 

