# SPAD counting model
This repository provides the software that goes with the publication [*Counting Statistics of Actively Quenched SPADs Under Continuous Illumination*](https://doi.org/10.1109/JLT.2020.2994654). It simulates the counting behavior of single-photon avalanche diodes (SPADs). Recovery time, afterpulsing and twilight pulsing are covered. Two commented programs are provided along with the necessary afterpulsing distributions.

Interactive CodeOcean capsule: https://doi.org/10.24433/CO.8487128.v1

If you have any questions or feedback, do not hesitate to contact me through the corresponding author's address or opening an issue here.

> Ivo Straka<br>
> ResearcherID V-2610-2017<br>
> ORCID 0000-0003-2675-6335

## mean_rate.py
Mean rate calculation implemented as a simple example written in Python. Numpy and SciPy are required.

## SPAD_simulation.c
A full simulation of counting behavior written in C. The output is counting statistics - the histogram of the number of detections per a 10-microsecond window. It assumes GCC with standard libraries; see headers. Makefile compilation is set to conventional compiling-machine-specific optimization.
