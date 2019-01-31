## Test environments
* local OS X install, R 3.5.2
* ubuntu 12.04 (on travis-ci), R 3.5.2
* win-builder (devel and release)
* r_hub() Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results
There were no ERRORs, WARNINGs, or NOTES in R 3.5.2 or R-devel on the first three systems.

On r_hub() there was 1 ERROR

  Package required and available but unsuitable version: 'dplyr'

This update to lmeresampler is in response to the update to dplyr, so this should be expected, and I am depending on the improved version.
