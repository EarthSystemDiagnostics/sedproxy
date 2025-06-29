# sedproxy 0.7.6
* fixed test for ggplot2 objects to address https://github.com/tidyverse/ggplot2/issues/6498
* fixed additional warnings for using .data in tidyselect
* fixed warning for using @docType "package"

# sedproxy 0.7.5.1
* fix bug in mixed layer modelling where top of core is assumed to be age == 0

# sedproxy 0.7.5
* depend on R >= 3.5.0
* remove tidyverse from suggests 
* move Shiny code to new package https://github.com/EarthSystemDiagnostics/shinysedproxy
* wrap stat::time() in as.numeric() to speed up %in% comparison

# sedproxy 0.7.4
* fix bug in ChunkMatrix where if a single timepoint is selected it would be off by one
* fix bug in BioturbationWeights to allow very small mixing depths with zero layer thickness
* fix bug in assigning points to the mixed layer when timepoint are not in increasing age order

# sedproxy 0.7.3
* refactor and rename growth rate function to be more intuitive

# sedproxy 0.7.2
* modifications to pass R CMD check
* depend on R >= 4.0.0

# sedproxy 0.7.1
* set growth to zero at temperatures below -2°C in growth_rate_l09_R (FAME)

# sedproxy 0.7.0
* feature - allow simulation of proxy in mixed layer

# sedproxy 0.6.6
* bugfix - allow duplicated requested timepoints
* additionally include n.samples in output

# sedproxy 0.6.5

* bugfix - fix case where bio.depth and layer width both == 0, single row of clim.signal should be used for each requested timepoint


# sedproxy 0.6.4

* implement growth rate threshold in growth_rate_l09_R (FAME)


# sedproxy 0.6.2

* fixed calibration parameters attribute

# sedproxy 0.6.1

* Release version for final copy edited version of CoTP paper
* Repository moved to Github


# sedproxy 0.6

* Fast version of ClimToProxyClim called internally when parameters are time-invariant
* Bugfix to calculation of bioturbation weights when layer.width == 0
* Other performance improvements
* Basic unit tests to detect changes to output for a standard simulation


# sedproxy 0.5.0

* Release version for revised CoTP Discussion paper
* Can now use dynamic habitat weights
* Now has explicit sensor stage with conversion from temperature to Mg/Ca or Uk'37 units
* Habitat weights can be calculated using temperature response function and parameterisation from FAME 1.0 module
* Calibration uncertainty modelled via sampling parameters from fitted calibration regression model for each replicate


# sedproxy 0.4.0

* Release version for CoTP Discussion paper


# sedproxy 0.3.1

* Added a `NEWS.md` file to track changes to the package.
* ClimToProxyClim now takes a ts object for the input climate signal.
* Bioturbation weights now take the thickness of the layer from which samples were picked/extracted into account when determining the time period over which a proxy integrates the climate signal. This is controlled by the new argument `layer.width`

