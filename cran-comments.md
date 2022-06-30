## Resubmission
This is a resubmission. In this version I have:

* Expanded the description in DESCRIPTION to a full paragraph.

* Checked all exported functions/methods for valid \value tags in their
.Rd files and added or expanded them where needed. 


## Checks

### R CMD check results

0 errors  | 0 warnings  | 0 notes 


### Online checks seem OK:

- check_rhub()

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Andrew Dolman <andrew.dolman@awi.de>'

New submission

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

