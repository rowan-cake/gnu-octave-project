# gnu-octave-project

This repo is a small Octave prototype for linking a seasonal epidemic model to a hospital surge-bed decision problem.

## Current Model

The epidemiology layer is a seasonal SIR model:

- `S(t)` = susceptible population
- `I(t)` = infected population
- `R(t)` = recovered / removed population

The state equations are:

```math
\frac{dS}{dt} = -\frac{\beta(t)}{N}IS
```

```math
\frac{dI}{dt} = \frac{\beta(t)}{N}IS - \gamma I
```

```math
\frac{dR}{dt} = \gamma I
```

with seasonal transmission rate

```math
\beta(t) = \beta_{mean} + \beta_{amp}\sin\left(\frac{2\pi t}{T}\right)
```

In the current code:

- `beta_mean = 0.15`
- `beta_amp = 0.10`
- `T = 50`
- `gamma = 0.05`

so `beta(t)` varies between `0.05` and `0.25`. These values were not informed by data but just toy ones to create a model.

## Running the Model

From the repo root:

```bash
octave run_sir.m
```

This regenerates:

- `sir_octave.out`
- `sir_plot.png`

To build and solve the surge-bed MDP:

```bash
octave run_mdp.m
```

This regenerates:

- `mdp_results.mat`
- `mdp_policy.png`

## Goal

The broader goal is not only to simulate disease spread, but to support the hospital decision:

`Do we open more surge beds right now?`

To answer that, the SIR model is used as a forecasting layer, and a Markov Decision Process (MDP) is used as the decision layer(like the research project would do).

## MDP Formulation

We model the hospital state as:

```math
s = (n, k)
```

where:

- `n` = current patient count requiring beds(informed through SIR model)
- `k` = current active surge-capacity level

This state is intentionally hospital facing. The MDP does not directly use `(S, I, R)` as its state. Instead, the epidemic model feeds the hospital model through predicted arrivals.

### State Space `S`

- `n in {0, 1, ..., Nmax}`
- `k in {0, 1, ..., Kmax}`

So each state tells us both current demand and currently activated extra capacity.

### Action Space `A`

At each decision time, the hospital manager chooses:

- `a = -1`: contract surge capacity by one level
- `a = 0`: maintain the current surge level
- `a = +1`: expand surge capacity by one level

The updated surge level is

```math
k' = \min(\max(k + a, 0), K_{max})
```

where:

- `k + a` applies the decision to the current surge level
- `max(k + a, 0)` prevents the model from dropping below zero active surge levels
- `min(..., Kmax)` prevents the model from opening more than the maximum allowed surge level

So this formula "clips" the next surge setting into the feasible range:

- if the model tries to contract below `0`, it stays at `0`
- if the model tries to expand above `Kmax`, it stays at `Kmax`
- otherwise, the surge level changes normally by one step

### Transition Model `P`

The transition model connects the epidemic layer to the hospital layer in a simple way:

- SIR gives infections
- `60%` of infected people generate hospital demand
- hospitalized patients leave at an average rate of `4` days
- the MDP decides whether to open or close surge beds

Let:

- `I(t)` be the infected fraction from the SIR model
- `N` be the total population size
- `n_t` be the current number of hospitalized patients

We model new arrivals during time step `t` as a Poisson random variable:

```math
A_t \sim \text{Poisson}(\lambda_t)
```

with expected arrival rate

```math
\lambda_t = 0.60 \, N \, I(t)
```

This means hospital demand rises and falls with the infected population.

We model departures from the hospital using a simple 4-day average length of stay. With a 1-day time step, the discharge rate is:

```math
\mu = \frac{1}{4}
```

so the expected number of patients leaving is:

```math
D_t = \mu n_t = \frac{1}{4} n_t
```

The patient count then updates as:

```math
n_{t+1} = \min(N_{max}, \max(0, n_t + A_t - D_t))
```

and the surge level updates as:

```math
k_{t+1} = \min(\max(k_t + a_t, 0), K_{max})
```

So each step works like this:

1. The SIR model gives the current infected level `I(t)`
2. That infected level determines expected hospital arrivals
3. Existing hospitalized patients leave at a 4-day average rate
4. The MDP updates the surge-bed level based on the chosen action
5. The next state becomes `(n_{t+1}, k_{t+1})`

### Reward / Cost Function `R`

The MDP should balance three competing pressures:

1. The cost of keeping surge beds open
2. The penalty for having more patients than available capacity
3. The cost of changing surge levels too often

We write the one step cost as:

```math
g(n, k, a) = c_1 k' + c_2 \max(0, n - C(k')) + c_3 |a|
```

where:

- `k'` is the post-action surge level
- `C(k')` is total available capacity under surge level `k'`
- `c1` = operating cost of active surge capacity
- `c2` = congestion penalty for unmet bed demand
- `c3` = switching cost for opening or closing surge beds

The most important modeling choice is that congestion should be penalized much more heavily than keeping surge beds open. This encourages the policy to expand before unsafe overload occurs.

## Interpretation

This creates a clean separation between the epidemic and decision pieces:

- SIR layer: forecasts disease pressure over time
- Hospital layer: converts that pressure into patient arrivals
- MDP layer: chooses whether to expand, maintain, or contract surge capacity

The final output of the MDP will be a policy that maps hospital state to action, for example:

- if patient count is high and surge level is low, expand
- if patient count is moderate and surge is already high, maintain
- if patient count is low and surge is high, contract

## Prototype Parameters

To build and visualize the first MDP, we need a small set of concrete parameter values. For now, these should be treated as arbitrary prototype values chosen to make the model easy to run and interpret. They are not calibrated to a real hospital.

Values:

- `N = 100`: total population represented by the SIR model
- `Nmax = 80`: maximum patient count tracked in the MDP state space
- `Kmax = 2`: maximum surge level
- `Cbase = 20`: baseline bed capacity with no surge activated
- `C(k) = Cbase + 10k`: each surge level adds `10` beds

This gives:

- `C(0) = 20`
- `C(1) = 30`
- `C(2) = 40`

For the one-step cost function

```math
g(n, k, a) = c_1 k' + c_2 \max(0, n - C(k')) + c_3 |a|
```

suggested prototype weights are:

- `c1 = 2`: operating cost for each active surge level
- `c2 = 25`: congestion penalty for patients above available capacity
- `c3 = 5`: switching cost for changing surge levels

These values are chosen to reflect the idea that:

- overload is much worse than running extra surge beds
- opening or closing surge beds has a real cost
- the model should avoid both unsafe crowding and unnecessary switching

This parameter set is mainly intended to support a first working MDP and a clear policy visualization.


## Reading the Plots

`sir_plot.png` shows the seasonal SIR simulation over time:

- `S` decreases as more people become infected
- `I` rises and falls in seasonal waves
- `R` increases as infected people recover and leave the infected group

`mdp_policy.png` shows the surge-bed policy chosen by the MDP:

- each panel corresponds to a current surge level `k = 0, 1, 2`
- the x-axis is day and the y-axis is the current patient count `n`
- the color shows the best action: contract, maintain, or expand

So the heatmap can be read as: given the current day in the epidemic, the current patient load, and the current surge setting, what should the hospital do next?


#### Sources

- SIR model reference https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology). 
- Used this to setup the ODE and solve (http://epirecip.es/epicookbook/chapters/sir/octave).
- For learing about the MDP (https://en.wikipedia.org/wiki/Markov_decision_process).