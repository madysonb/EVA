## QCAA - Query, Calculate, Age, Analyze
### Variability-Age Analysis

Madyson G. Barber

This notebook computes the VarX values and VarX90 ages for a given list of stars using methods outlined in Barber & Mann (in review). We calculate expected Gaia magnitude uncertainties using [EDR3_Photometric_Uncertainties](https://github.com/gaia-dpci/gaia-dr3-photometric-uncertainties) from [Riello et al. 2021](https://www.aanda.org/articles/aa/full_html/2021/05/aa39587-20/aa39587-20.html).

A seperate python file (analysisFunctions.py) contains all function definitions used below, as well as some additional helper functions. The sample cuts, as outlined in Barber & Mann (in review), are included in the filters function. These cuts can be manipulated, but may yield differing results as these are the cuts the calibration is tuned to.

Function uses are detailed below. An example use is in Var90Diagnostic.ipynb

---
Function: query  
Description: Pull needed information from Gaia
- Inputs:
    - file path to csv including Gaia DR3 ids --> works best with FriendFinder output
    - rewrite: whether or not to re-query the stars -- Default False
        - True --> re-query stars
        - False --> grab reference to Gaia query 
- Outputs:
    - pandas dataframe of all stars in the sample with the right available data
    

Function: calculate  
Description: Basic var90 calculations
- Inputs:
    - pandas dataframe of gaia query results
    - Additional cut parameters 
        - rv --> checks abs(RV-Rv_pred) < Xkm/s -- Default False
            - ONLY AVAILABLE WHEN STARTING WITH FRIENDFINDER OUTPUT FILE
            - True --> X = 5
            - given value --> X = val
            - False or None --> do not cut
- Outputs:
    - adds varX values to dataframe
    - return varX90 values
    
    
Function: age  
Description: calculate age of group in given band or overall (combining all three bands)
- Inputs:
    - pandas dataframe of gaia query results
    - band --> which age to report -- Default overall
        - G
        - BP
        - RP
        - overall --> weighted average of the three bandpass filters
        - all --> return order: overall, G, BP, RP
    - distance --> Include distance calibrations -- Default True
        - True --> use median distance of stars
        - given value --> override median dist
        - False --> age without distance
    - rv --> checks abs(RV-Rv_pred) < Xkm/s -- Default False
        - ONLY AVAILABLE WHEN STARTING WITH FRIENDFINDER OUTPUT FILE
        - True --> X = 5
        - given value --> X = val
        - False or None --> do not cut
- Outputs:
    - return desired varX90 age(s)
 

Function: analyze  
Description: make plots
- Inputs:
    - path to file with varX values
    - Desired bandpass ('G', 'BP', 'RP') -- Default G
- Outputs:
    - XYZ colored by varX
    - histogram of varX

---
### Dependencies  
- astropy  
- astroquery  
- galpy  
- matplotlib  
- numpy  
- pandas  
- scipy  
