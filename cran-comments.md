## Resubmission
This is a package update to address changes to testing for ggplot objects, see 
https://github.com/tidyverse/ggplot2/issues/6498

# sedproxy 0.7.5.2

* fixed test for ggplot2 objects to address https://github.com/tidyverse/ggplot2/issues/6498

Additional changes

* fixed additional warnings for using .data in tidyselect
* fixed warning for using @docType "package"
* fix bug in mixed layer modelling where top of core is assumed to be age == 0


## Checks

### R CMD check results

── R CMD check results ────── sedproxy 0.7.5.2 ────
Duration: 1m 18.2s

0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded


### Online checks:

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

