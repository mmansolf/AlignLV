# AlignLV
Allows for multiple group factor analysis alignment a la Mplus to be applied to lists of single-group models estimated in lavaan or mirt

## Credit line
The alignment algorithm is described well on Mplus's website:
https://www.statmodel.com/Alignment.shtml

# Description
Coming soon!

# Example

# Installing from source
We recommend installing the latest development version of this package from Github to ensure use of the latest and greatest version. To do so:
1. Install the `devtools` package if you haven't. In R, this can be done using the following code:
`install.packages('devtools')`
2. Load the `devtools` package using the following code:
`library(devtools)`
3. Install the `AlignLV` package using `install_github()` using the following code:
```
devtools::install_github('AlignLV',build_vignettes=T)
```
The `,build_vignettes=T` is optional but recommended for viewing the vignette accompanying this package (COMING SOON).
The `AlignLV` package is not yet on CRAN, so trying to install it with `install.packages()` will not work. This is coming soon!

## Other installation options
While the above method should work for most users, there are alternatives:
* The `ghit` library has an analogous `ghit::install_github()` function to that in the `devtools` package
* Download the package as a .zip file, then run the following code and interactively select the .zip file to install:
`install.packages(file.choose(), repos = NULL, type = "source")`

# Additional license information
LICENSE: Creative Commons - Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
* For more information on this license type, see https://creativecommons.org/licenses/by-nc/4.0/
* For more information on licensing for this package in particular, see the included LICENSE.md file.

# Bugs and questions
Bug reports and feedback are always welcome. For reporting bugs and requesting minor improvements, please use the "Issues" functionality in Github. For larger requests for improvement or scholarly co-operation in expanding the true score imputation framework, please email me directly at maxwell.mansolf@northwestern.edu. When in doubt, either contact method is fine.
