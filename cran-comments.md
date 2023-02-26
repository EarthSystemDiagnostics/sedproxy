## Resubmission
This is a package update to address performance issues on r-devel. 

* issue involved comparing a timeseries object and numeric vector with %in%
* solved by wrapping stat::time() in as.numeric() to speed up %in% comparison

Additional changes

* R version requirement dropped to depend on R >= 3.5.0
* Remove tidyverse from suggests 
* To reduce dependencies, Shiny code moved to new package https://github.com/EarthSystemDiagnostics/shinysedproxy


## Checks

### R CMD check results

-- R CMD check results ----------- sedproxy 0.7.5 ----
Duration: 1m 54.9s

0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded


### Online checks:

Rhub has not been available for several days.

R CMD Check on Github completes with no errors on the following platforms. 

  - {os: macOS-latest,   r: 'release'}
  - {os: windows-latest, r: 'release'}
  - {os: ubuntu-latest,   r: 'devel', http-user-agent: 'release'}
  - {os: ubuntu-latest,   r: 'release'}
  - {os: ubuntu-latest,   r: 'oldrel-1'}



### Results from prior submission 

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Andrew Dolman <andrew.dolman@awi.de>'

Possibly misspelled words in DESCRIPTION:
  Alkenones (7:21)
  Dolman (10:39)
  Laepple (10:50)
  Sedproxy (9:39)

* These are correct spellings of a chemical, author names, and the package name.

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1029/2002GL016612
    From: man/N41.proxy.Rd
          inst/doc/introduction-to-sedproxy.html
    Status: 503
    Message: Service Unavailable
    
* I have checked this DOI is correct and working.

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
  
* This appears only on Windows platform checks.
  
  
## Downstream dependencies

There are no downstream dependencies

