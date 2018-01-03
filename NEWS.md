# sedproxy 0.3.1

* Added a `NEWS.md` file to track changes to the package.
* ClimToProxyClim now takes a ts object for the input climate signal.
* Bioturbation weights now take the thickness of the layer from which samples were picked/extracted into account when determining the time period over which a proxy integrates the climate signal. This is controlled by the new argument `layer.width`


