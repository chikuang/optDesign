# Directional/partial derivative of each of the linear/non-linear regression models

---

## Models

1. Peleg model

Peleg is a statistical model used to investigate the relationship of water absorption for various kinds of food. This model is given by
```math
y_i=y_o+\frac{x_i}{\theta_1+\theta_2 x_i}+\epsilon_i,\quad i=1,2,\cdots,n,
```
where $`y_o`$ represents the initial moisture of the food, $`y_i`$ is the current level of moisture at current time $`x_i`$, $`\theta_1`$ is the Peleg's moisture rate constant, and $`\theta_2`$ is the asymptotic moisture as time increases. 

2. Polynomial model

A polynomial regression model of degree $q$ ($q \geq 1$) without intercept is given by
```math
y_i=\theta_1 x_i+\theta_2 x_i^2+\cdots+\theta_q x_i^q+\epsilon_i,~x_i\in S=[-1,+1],~i=1,2,\cdots,n.
```
Polynomial regression models are widely used when the response and regressors have a curvilinear relationship. Complex nonlinear relationships can be well approximated by polynomials over a small range of explanatory variables  (Montgomery et al., 2012, p. 223). 

3. EMax model

The model is given by
```math
y_i=\theta_1 e^{-\theta_2e^{-\theta_3 x_i}}+\epsilon_i,~i=1,2,\cdots,n,~\boldsymbol{\theta}
=(\theta_1,\theta_2,\theta_3)^\top,~x_i\in S,
```
where $`\theta_1`$ describes the maximum growing capacity, $`\theta_2`$ explains the initial status of the subject, $`\theta_3`$ determines the growth rate, $`y`$ is the overall growth at the current time point and $`x`$ is the time. Note that $`x`$, $`\theta_1`$, $`\theta_2`$ and $`\theta_3`$ are assumed to be positive in this context. We want to study how one subject's total growth is associated with time. The model has broad applications in biological science and cancer studies. 

4. Michaelis-Menton model
Michaelis-Menton model is one of the best-known models proposed by Leonor Michaelis and Maud Menten \citep{michaelis::menton1913kinetik} that studies the enzyme reactions between the enzyme and the substrate concentration. This model is given by $`y_i=\frac{\theta_1x_i}{\theta_2+x_i}+\epsilon_i.~ i=1,2,\cdots,n,~ x_i \in(0,                                                                k_o],~\theta_1,\theta_2\geq 0.`$ Here, $`\boldsymbol{\theta}=(\theta_1,\theta_2)^\top`$. Define $`z=\frac{1}{\theta_2+x}`$. Then 
$`\boldsymbol{f}(x,\boldsymbol{\theta})=(\frac{x}{\theta_2+x},\frac{-\theta_1x}{(\theta_2+x)^2})^\top= (1-z\theta_2, -z\theta_1 +\theta_1 \theta_2z^2)^\top.`$

## TODO
+ [ ] Generalised linear model (GLM)
