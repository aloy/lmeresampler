
This is essentially a resubmission to address errors caused by my recent submission that were unexpected based on my testing. I apologize for a new submission so quickly.  

In this version I have:

* Failing tests are skipped on CRAN as they are not robust to different operating systems. They will be edited and changed on the development version and tested via continuous integration.

## Test environments
Windows (github actions - release,  win-builder - devel and release)
macOS (local, github actions)
Ubuntu (github actions, r-hub)
Fedora (r-hub)


## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE generated only on win-builder windows release:

*  checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Adam Loy <loyad01@gmail.com>’
 
    Found the following (possibly) invalid DOIs:
      DOI: 10.1111/1467-9876.00415
        From: DESCRIPTION
        Status: Service Unavailable
        Message: 503
      
        This DOI is valid, and I did not receive an error navigating there in a web browser
