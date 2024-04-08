# Kalman filtering to reduce measurement noise of sample entropy: An electroencephalographic study

This is the matlab code implementation of our paper "Kalman filtering to reduce measurement noise of  sample entropy: An electroencephalographic study"

## Datasets
The data used in this study includes the [sleep signals](https://physionet.org/content/sleep-edfx/1.0.0/) , [epilepsy signals](https://physionet.org/content/chbmit/1.0.0/), and synthetic signals (white noise, logistic map, rossel system). 

- Sleep signals: [Sleep-EDF Database Expanded](https://physionet.org/content/sleep-edfx/1.0.0/). DOI: https://doi.org/10.13026/C2X676
- Epilepsy signals: [CHB-MIT Scalp EEG Database](https://physionet.org/content/chbmit/1.0.0/). DOI: https://doi.org/10.13026/C2K01R

## Prerequisites
-  Matlab
-  [NNetEn calculator 1.0.0.4](https://www.researchgate.net/publication/366575849_NNetEn_calculator_1004_for_NNetEn_calculation_with_six_matrix_filling_methods)


## How to use this method?
- What if you want to use this method on sleep signals? Please execute this code.
```bash
run kalman_entropy_sleep.m
```

- What if you want to use this method on epilepsy signals? Please execute this code.
```bash
run kalman_entropy_epilepsy.m 
```


- What if you want to use this method on white noise signals? Please execute this code.
```bash
run kalman_entropy_whiteNoise.m
```

- What if you want to use this method on logistic map signals? Please execute this code.
```bash
run kalman_entropy_logisticMap.m
```

- What if you want to use this method on rossel system signals? Please execute this code.
```bash
run kalman_entropy_rossel_System.m 
```

