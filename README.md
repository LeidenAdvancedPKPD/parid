# parid: Calculation of parameter identifiability indicators for nonlinear mixed effects models

Methods for determining the identifiability of nonlinear mixed effects model parameters, namely:
* Sensitivity Matrix Method (SMM)
* Fisher Information Matrix Method (FIMM)
* Aliasing
The model is given as a set of differential equations, with specification of observed outputs, parameters and dosing.
The methods assess whether the parameters can be identified from the model output, and optionally indicate dependencies between parameters.

## Installation

The package can be installed using:

```R
devtools::install_github("LeidenAdvancedPKPD/parid").
```

## Documentation

A worked out [example](https://github.com/LeidenAdvancedPKPD/parid/blob/main/tests/example.r) is available in the package.
A showcase can be found [here](https://lapp.nl/lapp-software/parid.html).

## License

[GNU GPL 3.0 license](https://www.gnu.org/licenses/gpl-3.0.html)
