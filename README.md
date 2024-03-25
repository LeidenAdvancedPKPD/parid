# parid: Calculation of parameter identifiability indicators for nonlinear mixed effects models

Methods for determining the identifiability of nonlinear mixed effects model parameters, namely:

- Sensitivity Matrix Method (SMM)
- Fisher Information Matrix Method (FIMM)
- Aliasing

The model is given as a set of differential equations, with specification of observed outputs, parameters and dosing.
The methods assess whether the parameters can be identified from the model output, and optionally indicate dependencies between parameters.

## Unique features of the identifiability methods

The parameter identifiability methods in this package provide the following unique features.

- They take data limitations into account. That is, identifiability is assessed for a given set of subjects and their sample times. This is contrast to many other methods that assume that samples are available at all times. The parid methods let the user choose a discrete set of times, which is more in line with realistic applications.
- They provide categorical and continuous identifiability indicators. Categorical indicators provide a yes/no answer. Many other methods provide merely this, and while categorical indicators seem clear, it may be misleading, because model parameters may be formally identifiable, but hard to identify in practice. For that reason, continuous indicators are provided as well. They determine a level of identifiability. For the SMM and FIMM, this is a number between 0 and 1 where 0 means unidentifiable and non-zero means formally identifiable. However, small values indicate a low level of identifiability, which may lead to difficulties in estimating the parameters in practice. In this case, the influence of different parameters may be formally different but hard to distinguish on the chosen data set.
- They assess the identifiability of the parameter vector, not the individual parameters. That is, unidentifiability resulting from a combination of parameters can be detected. As a trivial example, a model with one data point and two parameters would not be identifiable. Some methods only consider individual parameters and might conclude that both parameters are identifiable (assuming tacitly that the other parameter is fixed).
- They do not need data values or a successful model fit. That is, they can be used a priori (before a model fit has been achieved) and even before the data has become available. The methods require merely the study design (sampling times and dosing information).

## Brief explanation of the three methods

- Sensitivity Matrix Method (SMM): Uses the sensitivity matrix dy/dp, where y is the model output and p are the parameters. The null space of this matrix characterizes the identifiability and determines the unidentifiable directions. Several continuous indicators have been developed based on this matrix. This method cannot take random-effects parameters into account.
- Fisher Information Matrix Method (FIMM): Uses the Fisher information matrix, which describes how the log-likelihood changes as function of the parameters. The eigensystem of this matrix characterizes the identifiability and determines the unidentifiable directions. Several continuous indicators have been developed based on this matrix. This method does take random-effects parameters into account.
- Aliasing: developed by [Augustin et al](https://github.com/mathworks-SimBiology/AliasingScoreApp/blob/master/AliasingScore_Poster.pdf) in Matlab, and reimplemented here in R.


For more details see the [Documentation section](#documentation).

## Installation

The package can be installed using:

```R
devtools::install_github("LeidenAdvancedPKPD/parid")
```

## Documentation

More information can be found on the [pkgdown](https://leidenadvancedpkpd.github.io/parid/) site.
A worked out [example](https://github.com/LeidenAdvancedPKPD/parid/blob/main/tests/example.r) is available in the package. This is the same example as on the pkgdown site.
A showcase can be found [here](https://lapp.nl/lapp-software/parid.html). This demonstrates a different example.

## License

[GNU GPL 3.0 license](https://www.gnu.org/licenses/gpl-3.0.html)
