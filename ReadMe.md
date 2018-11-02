# Generating other versions of Figure 1 of report

In all cases the *R*<sup>2</sup> values are quite low (< 0.5).
Hence, I believe that if the observational datasets yield signinficantly higher *R*<sup>2</sup> values (e.g., > 0.8) for the same test, then the argument made in the report does not hold.

## Using uniform distributions

![image1a](https://github.com/briochemc/TellurideFigure/blob/master/dissolved_Fe_fraction_from_uniform_DFe_and_FeT.png)

Figure generated with uniform distributions of DFe and Fe<sub>T</sub> (then removing samples with DFe > Fe<sub>T</sub>).

## Using lognormal distributions

### lognormal DFe and Fe<sub>T</sub> that matches the mean and variance of the uniformly generated DFe and Fe<sub>T</sub> above

![image1b1](https://github.com/briochemc/TellurideFigure/blob/master/dissolved_Fe_fraction_from_lognormal_DFe_and_FeT.png)

Figure generated with lognormal distributions of DFe and Fe<sub>T</sub> (then removing samples with DFe > Fe<sub>T</sub>, but should not happen).

### lognormal DFe and Fe<sub>T</sub>, but with a higher mean and variance for DFe

![image1b2](https://github.com/briochemc/TellurideFigure/blob/master/dissolved_Fe_fraction_from_lognormal_DFe_and_FeT2.png)

Figure generated with lognormal distributions of DFe and Fe<sub>T</sub> (then removing samples with DFe > Fe<sub>T</sub>).

### lognormal DFe and PFe (more natural?)

![image1c](https://github.com/briochemc/TellurideFigure/blob/master/dissolved_Fe_fraction_from_lognormal_DFe_and_PFe.png)

Figure generated with uniform distributions of DFe and PFe (and then Fe<sub>T</sub> = DFe + PFe, no need to remove samples).
This last figure seems more natural to me because I believe generating DFe and Fe<sub>T</sub> and then imposing DFe < Fe<sub>T</sub> biases the sample distributions.
