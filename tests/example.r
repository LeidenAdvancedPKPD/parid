# parid package: Calculation of parameter identifiability
#                indicators for nonlinear mixed effects models.
#
# Copyright (C) 2023 LAP&P Consultants BV
# Contact: info@lapp.nl
#          www.lapp.nl
#          Archimedesweg 31
#          2333 CM Leiden, The Netherlands
# 
# This file is part of the parid package.
# The parid package is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# The parid package is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the parid package.
# If not, see <https://www.gnu.org/licenses/>. 


# ---

# This file provides an example application of the parid package.
# The example is a one compartment linear PK model with absorption and bio-availability.
# Two versions are tested, namely with bioavailability fixed or variable.
# The first should be identifiable, and the second should not.
library(parid)

# Model definition.
# Define functions describing the model (differential equations), parameters, initials and model output:
model      <- function(t, y, p) { c(-p[["Ka"]]*y[["x1"]], p[["Ka"]]*y[["x1"]] - p[["CL"]]/p[["V"]] * y[["x2"]]) }
p          <- function(theta, eta) { c(Ka = theta[["TVKa"]], CL = theta[["TVCL"]] * exp(eta[["iCL"]]),
                                       V = theta[["TVV"]], F = theta[["TVF"]], Dose = theta[["TVDose"]])  }
init       <- function(p) { c("x1" = p[["F"]] * p[["Dose"]], "x2" = 0) }
output     <- function(y, p, eps) { y[["x2"]] / p[["V"]] + eps[["addErr"]] }
# Define parameters. For thetas, define which ones should be considered variable, and which ones fixed:
theta      <- c("TVKa" = 1, "TVCL" = 0.2, "TVV" = 2, "TVF" = 0.7, "TVDose" = 10)
vartheta1  <- setdiff(names(theta), c("TVDose", "TVF"))     # Bioavailability F is fixed
vartheta2  <- setdiff(names(theta), "TVDose")               # Bioavailability F is variable
omega      <- diag(0.3, nrow = 1)
colnames(omega) <- row.names(omega) <- c("iCL")
sigma      <- diag(0.1, nrow = 1)
colnames(sigma) <- row.names(sigma) <- c("addErr")
# Define observation times. Initialization is always at time 0:
times      <- seq(0, 10, 1)

# Compute and plot variational equations:
vareq1 <- calcVariationsFim(model = model, p = p, init = init, output = output, times = times, symbolic = TRUE,
                           theta = theta, nmeta = colnames(omega), nmeps = colnames(sigma), vartheta = vartheta1)
plotVariations(vareq1)   # This function has parameters controlling which derivatives are plotted. See the documentation for details.
vareq2 <- calcVariationsFim(model = model, p = p, init = init, output = output, times = times, symbolic = TRUE,
                            theta = theta, nmeta = colnames(omega), nmeps = colnames(sigma), vartheta = vartheta2)
plotVariations(vareq2)

# Compute and plot Sensitivity Matrix Method (SMM) results:
sens1 <- calcSensitivityFromMatrix(outputs = list("N", "A", "R", "L", "M"), df = vareq1, vars = vartheta1)
simplifySensitivities(sens1)  # Identifiable if F is fixed (simplification sets near-zero values to zero)
plotSensitivities(sens1)      # Plot R, M, L indicators. See the documentation for details.

sens2 <- calcSensitivityFromMatrix(outputs = list("N", "A", "R", "L", "M"), df = vareq2, vars = vartheta2)
simplifySensitivities(sens2)  # Not identifiable if F is variable, with a relation involving CL, V, F.
plotSensitivities(sens2)

# Compute and plot Aliasing Method results:
alia1 <- calcAliasingScoresFromMatrix(df = vareq1, vars = vartheta1)
plotAliasing(alia1, elt = c("S", "T"))  # Identifiable if F is fixed (highest score is low at 24%). Plot normalized sensitivities and aliasing score matrix.

alia2 <- calcAliasingScoresFromMatrix(df = vareq2, vars = vartheta2)
plotAliasing(alia2, elt = c("S", "T"))  # If F is variable, highest score is 54%, suggesting still identifiable --> aliasing can not find relations between >2 parameters.

# Compute and plot Fisher Information Matrix Method (FIMM) results:
fim1 <- calcFimFromMatrix(df = vareq1, omega = omega, sigma = sigma, vartheta = vartheta1)   # Compute FIM
fimres1 <- fimIdent(fim = fim1, curvature = 1e-10, relChanges = TRUE, ci = 0.95, nsubj = 10) # Compute FIMM results
simplifyFimIdent(fimres1)  # Identifiable if F is fixed

fim2 <- calcFimFromMatrix(df = vareq2, omega = omega, sigma = sigma, vartheta = vartheta2)   # Compute FIM
fimres2 <- fimIdent(fim = fim2, curvature = 1e-10, relChanges = TRUE, ci = 0.95, nsubj = 10) # Compute FIMM results
simplifyFimIdent(fimres2)  # Low curvature --> not identifiable if F is variable
