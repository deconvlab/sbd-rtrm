# Sparse blind deconvolution using the Riemannian Trust-Region Method (RTRM)
This package performs blind deconvolution under the data model
<p align="center">
    <img src="./docs/math/model.png"/ height="300">
</p>


## Setup
 1. Ensure the ManOpt package is installed for RTRM [(http://www.manopt.org)](http://www.manopt.org).
 2. Download the `sbd-rtrm` package, and in MATLAB run `init_sbd.m` each time you want to use the tools from this package.
 3. Run `examples/simple_SBD_example.m` or `full_SBD_and_STM_example.m`.

## Core functions
### `SBD.m`
### `core\Asolve_Manopt.m`
### `core\Xsolve_pdNCG.m`


## Future updates
 - Add `Xpos` option to `SBD.m`
 - Create a RTRM package for CDL