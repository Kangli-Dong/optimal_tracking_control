# Introduction

This code is to performing the optimal tracking control of brain network dynamics.

Network control theory (NCT) has recently been used in neuroscience to understand brain stimulation effects. This paper considers stochastic brain dynamics and introduces optimal stochastic tracking control to synchronize brain dynamics to target dynamics rather than to a target state at a specific time point. A revised gradient descent optimization was used to estimate parameters for brain dynamic systems, and optimal stochastic tracking control was used to drive unhealthy dynamics by controlling a few nodes to synchronize with healthy dynamics. Results show that the tracking energy is negatively correlated with the average controllability of the brain network system, while the energy of the optimal state approaching control is significantly related to the target state value. For a 100-dimensional system, controlling five nodes with the lowest tracking energy can achieve acceptable control effects. These findings highlight the potential of stochastic tracking control as a novel approach for future brain stimulation interventions.


## Citing

```
@article{dong2025new,
  title={A new perspective on brain stimulation interventions: Optimal stochastic tracking control of brain network dynamics},
  author={Dong, Kangli and Chen, Siya and Dan, Ying and Zhang, Lu and Li, Xinyi and Liang, Wei and Zhao, Yue and Sun, Yu},
  journal={arXiv preprint arXiv:2501.08567},
  year={2025}
}
```
