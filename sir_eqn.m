function xdot = sir_eqn(x, t)
  % SIR model with a seasonal transmission rate beta(t).
  %
  % The transmission rate is modeled as
  %   beta(t) = beta_mean + beta_amp * sin(2 * pi * t / T)
  % which creates repeating "surge" periods over time while keeping
  % transmission above a nonzero seasonal baseline.
  %
  % State equations:
  %   dS/dt = -(beta(t)/N) * I * S
  %   dI/dt =  (beta(t)/N) * I * S - gamma * I
  %   dR/dt =  gamma * I

  % Parameters
  beta_mean = 0.15;  % Midpoint of the seasonal transmission rate
  beta_amp = 0.10;   % Seasonal swing around the midpoint
  T = 50;            % Period of one seasonal cycle
  gamma = 0.05;

  % State variables
  S = x(1);
  I = x(2);
  R = x(3);

  % Total population
  N = S + I + R;

  % Seasonal transmission rate. This keeps beta(t) between 0.05 and 0.25:
  %   minimum = beta_mean - beta_amp = 0.05
  %   maximum = beta_mean + beta_amp = 0.25
  % When the sine term rises, infections spread more easily; when it
  % falls, transmission weakens but does not fully disappear.
  beta = beta_mean + beta_amp * sin(2 * pi * t / T);

  % ODE system
  dS = -(beta / N) * I * S;
  dI =  (beta / N) * I * S - gamma * I;
  dR =  gamma * I;

  % Return as a column vector for lsode
  xdot = [dS; dI; dR];
endfunction
