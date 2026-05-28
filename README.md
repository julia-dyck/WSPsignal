# WSPsignal <img src="man/figures/Logo_WSPsignal.png" align="right" style="height:170px; width:auto;">

The family of Weibull Shape Parameter (WSP) tests was developed to detect signals of adverse events in electronic health records.

Most recently, the BPgWSP test was developed and presented in

Dyck, J., & Sauzet, O. (2025). **The BPgWSP test: a Bayesian Weibull Shape Parameter signal detection test for adverse drug reactions.** arXiv preprint arXiv:2412.05463.

available on [https://arxiv.org/abs/2412.05463](https://arxiv.org/abs/2412.05463)

triggering the effort to gather all WSP tests ready to apply in one package. A preprint describing the R package in more detail will be available soon.

## Installation of the WSPsignal package
You can install the WSPsignal package with

``` r
remotes::install_github(repo = "julia-dyck/WSPsignal")
```

## Tuning of the WSP test for your purpose
To tune the WSP test for your application, you can use the R scripts on [https://osf.io/h7edn/overview](https://osf.io/h7edn/overview) as a template, and adjust sample scenarios as well as test alternatives accordingly.

