# Changelog

## [v0.2.0]

### Enhancements
* Added compatibility with FINAM v1.
* Added support for 3D structured and unstructured grids.
* Added support for masked input and output data.
* Added support for spherical/lat-lon regridding and configurable regridding CRS via `RegridCRS`.
* Added `zero_region` argument to the regridder.
* Re-exported additional ESMPy constants: `LineType`, `NormType`, `PoleMethod`, and `Region`.

### Fixes
* Fixed CRS transformations to preserve XY axis order.
* Fixed spherical LocStream coordinate key names.
* Treated unstructured grids with only point data as LocStreams.
* Fixed structured grid generation order for transformed grids.

### Dependencies
* Removed the direct `xarray` dependency.
* Raised minimum dependency versions to FINAM 1.1.0 and ESMPy 8.7.

## [v0.1.0]

* Initial release of finam-regrid

[unpublished]: https://git.ufz.de/FINAM/finam-regrid/-/compare/v0.2.0...main
[v0.2.0]: https://git.ufz.de/FINAM/finam-regrid/-/compare/v0.1.0...v0.2.0
[v0.1.0]: https://git.ufz.de/FINAM/finam-regrid/-/commits/v0.1.0
