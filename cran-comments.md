## Test environments



## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTES:

*  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Adam Loy <loyad01@gmail.com>'

    New submission
  
    Package was archived on CRAN
  
    Possibly misspelled words in DESCRIPTION:
      Goldstein (12:16)
      Leeden (11:16)
      REB (9:55)
      Rasbash (12:28)
      al (11:26)
      der (11:12)
      et (11:23)
      
        These are note misspelled words
  
    CRAN repository db overrides:
      X-CRAN-Comment: Archived on 2022-02-24 as requires archived package
        'catchr'.
        
        This submission removes the dependency on `catchr` to avoid future issues.
  
    Found the following (possibly) invalid DOIs:
      DOI: 10.1111/1467-9876.00415
        From: DESCRIPTION
        Status: Service Unavailable
        Message: 503
      
        This DOI is valid, and I did not receive an error navigating there in a web browser


## Downstream dependencies
I have also run R CMD check on downstream dependencies of httr 
(https://github.com/wch/checkresults/blob/master/httr/r-release). 
All packages that I could install passed except:
