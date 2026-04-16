function xdot = sir_eqn(x, t)
  % SIR model with explicit population scaling:
  % dS/dt = -(beta/N) * I * S
  % dI/dt =  (beta/N) * I * S - gamma * I
  % dR/dt =  gamma * I

  % Parameters
  beta = 0.1;   % Will become beta(t) later
  gamma = 0.05;

  % State variables
  S = x(1);
  I = x(2);
  R = x(3);

  % Total population
  N = S + I + R;

  % ODE system
  dS = -(beta / N) * I * S;
  dI =  (beta / N) * I * S - gamma * I;
  dR =  gamma * I;

  % Return as a column vector for lsode
  xdot = [dS; dI; dR];
endfunction
