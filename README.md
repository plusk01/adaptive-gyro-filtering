Adaptive Gyro Notch Filtering
=============================

Signal conditioning is an important part of a control strategy. In the context of multirotor control, the on-board gyro provides angular rate measurements which are directly used to drive the actual angular rate to the desired angular rate. In this setting, gyro noise is fed into the actuators, creating visible attitude fluctuations (i.e., "wiggles") even when the vehicle is meant to be hovering. Additionally, noisy gyro signals are directly related to hot motors (which can lead to motor failure, but is also an indication of performance). The majority of gyro noise comes from vibrations, caused by the spinning of the motors.

An easy way to reduce noise is with a low-pass filter. A commonly employed strategy in control is an [RC-type low-pass filter](https://en.wikipedia.org/wiki/Low-pass_filter), which is equivalent to a first-order IIR filter (e.g., see [here](http://www.olliw.eu/2016/digital-filters/#chapter22)). The tradeoff in filtering is always attenuation vs sample delay. A low cutoff frequency will more strongly attenuate high frequency noise at the cost of increased sample delay. Therefore, instead of designing a low-pass filter with a low cutoff frequency to attenuate high-frequency motor noise, a more targeted approach should be used.

This repo provides an implementation of an adaptive notch filter for gyro noise suppresion. It is largely inspired by [Betaflight](https://github.com/betaflight/betaflight), where the algorithm is meant to run on a resource constrained microprocessor. While still efficient, this implementation enjoys the benefit of easy debugging and visualizations to help tune such an algorithm.

## Example Application

<p align="center"><img src=".github/dynnotch.gif" width="100%" /></p>

## Technical Description

### Resources

- [Betaflight's `gyroanalyse.c`](https://github.com/betaflight/betaflight/blob/master/src/main/flight/gyroanalyse.c)
- Biquads: [OlliW](http://www.olliw.eu/2016/digital-filters)
- Biquads: [Nigel Redmon](https://www.earlevel.com/main/2012/11/26/biquad-c-source-code/)
- Betaflight filtering: [rav's Python code](https://github.com/rav-rav/betaflightCalc/tree/master/src)
