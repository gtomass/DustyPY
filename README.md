# DustyPY

DustyPY is a handy package for running and adjusting SEDs using the Dusty radiative transfer modeling code [Dusty on GitHub](https://github.com/ivezic/dusty.git) version 4.0 and the [pymcmcstat](https://github.com/prmiles/pymcmcstat.git) package for MCMC fit.

## Installation

To install DustyPY, use the following command:

```bash
git clone https://github.com/Radron74/DustyPY.git
pip install .
```
## Import

```
from DustyPY import *
```

## Use

You can find an exemple of DustyPY utilisation in the notebook folder

## Note

As of today, there is a bug in dustyV4 that has not yet been corrected. This implies that the is not taken into account. Refer to the [raised issue](https://github.com/ivezic/dusty/issues/11) for solution.

## License

DustyPY is licensed under the MIT License. See the [LICENSE](https://github.com/Radron74/DustyPY/blob/main/LICENSE) file for more details.