## QCAA - Query, Calculate, Age, Analyze
### Variability Analysis

Madyson G. Barber

This notebook computes the VarX values and VarX90 ages for a given list of stars using methods outlined in [Barber & Mann (2023)](paperlink). We calculate expected Gaia magnitude uncertainties using EDR3_Photometric_Uncertainties from [Riello et al. 2021](https://www.aanda.org/articles/aa/full_html/2021/05/aa39587-20/aa39587-20.html).

A seperate python file (analysisFunctions.py) contains all function definitions used below, as well as some additional helper functions. The sample cuts, as outlined in [Barber & Mann (2023)](paperlink), are included in the filters function. These cuts can be manipulated, but may yield differing results as these are the cuts the calibration is tuned to. 

Function uses are detailed below.

---
Function: query  
Description: Pull needed information from Gaia
- Inputs:
    - file path to csv including Gaia DR3 ids --> works with FriendFinder output
    - rewrite: whether or not to re-query the stars
        - True --> re-query stars
        - False --> grab reference to Gaia query --> defaults to False
- Outputs:
    - pandas dataframe of all stars in the sample with the right available data
    

Function: calculate  
Description: Basic var90 calculations
- Inputs:
    - pandas of gaia query results
    - Additional cut parameters 
        - voff --> abs(RV-Rv_pred) < Xkm/s
            - "on" --> X = 5
            - given value --> X = val
            - "off" or False or None --> do not cut
            - need to add -- need to get RV during Gaia query, need and input of expected RV
- Outputs:
    - adds varX values to dataframe
    - return varX90 values
    
    
Function: age  
Description: calculate age of group in given band or overall (combining all three bands)
- Inputs:
    - path to file with varX values
        - check to see if file exists -- call calculate if not
    - band -- G, BP, RP, overall, all
    - distance --> True/given value, False --> defaults to True
        - True --> use median distance of stars
        - given value --> override median dist
        - False --> age without distance
- Outputs:
    - return varX90 ages and overall age
 

Function: analyze  
Description: make plots
- Inputs:
    - path to file with varX values
    - Desired bandpass ('G', 'BP', 'RP')
- Outputs:
    - XYZ colored by varX
    - histogram of varX

---
