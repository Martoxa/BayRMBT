# BayRMBT
BayRMBT is the R version of BayMBT (Tierney, Jessica E., 2021. BayMBT. https://github.com/jesstierney/BayMBT), a Bayesian calibration for the branched GDGT MBT5Me proxy in soils, peats, and lakes. In this package you'll find the calibration data, the code to construct the calibration, code to forward model MBT5Me from temperatures, and code to predict temperatures from MBT5Me. Please cite the following publication when using BayMBT for soils and peats:

Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M., & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model for branched glycerol dialkyl glycerol tetraethers in soils and peats. Geochimica et Cosmochimica Acta, 268, 142-159. https://doi.org/10.1016/j.gca.2019.09.043

And the following publication when using BayMBT for lakes:

Martínez-Sosa, P., Tierney, J. E., Stefanescu, I. C., Dearing Crampton-Flood, E., Shuman, B. N., Routson, C. A global Bayesian temperature calibration for lacustrine brGDGTs. Geochimica et Cosmochimica Acta https://doi.org/10.1016/j.gca.2021.04.038

#Description of contents

bayRmbt_model.R This contains the code to construct the calibration model. It loads either the soil or lake core top dataset and calculates model parameters using Bayesian linear regression. For everyday use you don't need to run this as the parameters are already calculated and provided here as .RDS files (params.RDS files). This is provided here for transparency. This function in R has three dependencies: the packages MASS, Rlab, and ks. 

bayRmbt_forward.R This calculates MBT5Me values from user-provided temperatures, using the parameters in the params.RDS files. There are no dependencies.

bayRmbt_predict.R This calculates temperatures from user-provided MBT5Me values, using the parameters in the params.RDS files, through Bayesian inference. It requires a prior mean and standard deviation on temperature. One can use either the "BayMBT" model (which is calibrated to mean annual temperatures) or the "BayMBT0" model (which is calibrated to mean annual temperatures above 0C). For both soils and lakes, we recommend using BayMBT0 for most applications, as it is a better fit to the core top data. The implicit assumption in BayMBT0 is that brGDGTs are not synthesized at temperatures below freezing.

To cite the code repository directly use:

Tierney, Jessica E., 2021. BayMBT. https://github.com/jesstierney/BayMBT.
