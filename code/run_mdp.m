% Build and solve a simple hospital surge-bed MDP driven by the seasonal
% SIR model. The state is s = (n, k), where:
%   n = number of hospitalized patients currently needing beds
%   k = active surge level
%
% The action is one of:
%   -1 = contract a surge level
%    0 = maintain the current surge level
%   +1 = expand a surge level
%
% This script uses a finite-horizon dynamic program because the arrival
% pressure changes over time as the seasonal SIR infection curve changes.

% Regenerate the SIR output if it does not exist yet.
if (exist("sir_octave.out", "file") != 2)
  run("run_sir.m");
endif

% Load the SIR output matrix with columns:
%   time, susceptible, infected, recovered
sir_data = load("-ascii", "sir_octave.out");

% Use one decision epoch per day. The SIR simulation is continuous in time,
% so we interpolate the infected fraction onto integer day values.
days = 1:200;
infected_daily = interp1(sir_data(:,1), sir_data(:,3), days, "linear");

% ------------------------
% Prototype MDP parameters
% ------------------------

% These are arbitrary prototype values for a first working model rather
% than calibrated hospital parameters.
N = 100;            % Population size represented by the SIR model
Nmax = 80;          % Maximum hospital census tracked by the MDP
Kmax = 2;           % Highest allowed surge level
Cbase = 20;         % Baseline staffed bed capacity
surge_step = 10;    % Additional beds opened by each surge level

% Hospital-flow assumptions discussed in the README.
p_hosp = 0.60;      % Fraction of infected people that need hospital care
mu = 1 / 4;         % Average 4-day stay -> 25% leave per day on average

% One-step cost weights.
c1 = 2;             % Cost of keeping surge beds active
c2 = 25;            % Penalty for patients above available capacity
c3 = 5;             % Cost of changing surge level
discount = 0.98;    % Mild discount on future cost

actions = [-1, 0, 1];
horizon = numel(days);

% Expected hospital arrivals on each day. This is the direct bridge from
% the SIR model to the hospital MDP.
lambda_daily = p_hosp * N * infected_daily;

% Value function and policy arrays:
%   V(t, n+1, k+1) stores the minimum expected cost-to-go from day t
%   policy(t, n+1, k+1) stores the best action at that state and day
V = zeros(horizon + 1, Nmax + 1, Kmax + 1);
policy = zeros(horizon, Nmax + 1, Kmax + 1);

% Backward dynamic programming over the finite time horizon.
for t = horizon:-1:1
  lambda_t = lambda_daily(t);

  % Truncate the Poisson arrival distribution to a manageable number of
  % buckets. The final bucket captures the entire right tail and maps to
  % the capped state n = Nmax.
  arrival_probs = poisson_tail_probs(lambda_t, Nmax);

  for n = 0:Nmax
    % Keep the discharge model deliberately simple: with a 4-day average
    % stay, the expected retained census after one day is 75% of n.
    retained_patients = round((1 - mu) * n);
    next_n_lookup = min(Nmax, max(0, retained_patients + (0:Nmax)));

    for k = 0:Kmax
      best_cost = Inf;
      best_action = 0;

      for action_index = 1:numel(actions)
        a = actions(action_index);

        % Apply the action and clip the result into the feasible surge
        % range, as described in the README.
        k_next = min(max(k + a, 0), Kmax);
        capacity = Cbase + surge_step * k_next;

        % The immediate cost penalizes active surge beds, unsafe overload,
        % and changing surge levels too often.
        overflow = max(0, n - capacity);
        immediate_cost = c1 * k_next + c2 * overflow + c3 * abs(a);

        % The only randomness here comes from arrivals, so the expected
        % future cost is just a weighted average over the arrival buckets.
        future_values = squeeze(V(t + 1, next_n_lookup + 1, k_next + 1))(:);
        expected_future_cost = arrival_probs' * future_values;

        total_cost = immediate_cost + discount * expected_future_cost;

        if (total_cost < best_cost)
          best_cost = total_cost;
          best_action = a;
        endif
      endfor

      V(t, n + 1, k + 1) = best_cost;
      policy(t, n + 1, k + 1) = best_action;
    endfor
  endfor
endfor

% Save the main MDP outputs for later inspection or reuse.
save("-mat", "mdp_results.mat", "days", "infected_daily", "lambda_daily", ...
     "policy", "V", "N", "Nmax", "Kmax", "Cbase", "surge_step", ...
     "p_hosp", "mu", "c1", "c2", "c3", "discount");

% --------------------------
% Create a policy heatmap PNG
% --------------------------

figure(2);
clf;
set(gcf, "position", [100, 100, 950, 900]);

for k = 0:Kmax
  subplot(Kmax + 1, 1, k + 1);

  % For each current surge level, show the best action over:
  %   x-axis = day in the seasonal epidemic
  %   y-axis = hospital census n
  % The color shows whether the MDP wants to contract, maintain, or expand.
  imagesc(days, 0:Nmax, squeeze(policy(:,:,k + 1))');
  axis xy;
  caxis([-1 1]);
  colormap(gca, [0.84 0.24 0.31; 0.92 0.92 0.92; 0.20 0.47 0.72]);
  colorbar("ytick", [-1, 0, 1], "yticklabel", {"contract", "maintain", "expand"});
  ylabel("Patients n");

  if (k < Kmax)
    set(gca, "xticklabel", []);
  else
    xlabel("Day");
  endif

  title(sprintf("Optimal Action Map for Current Surge Level k = %d", k));
endfor

print("mdp_policy.png", "-dpng");

% Print a compact summary in the terminal for quick inspection.
[peak_I, peak_index] = max(infected_daily);
printf("Saved MDP outputs to mdp_results.mat\n");
printf("Saved policy figure to mdp_policy.png\n");
printf("Peak infected fraction occurs on day %d with I(t) = %.4f\n", ...
       days(peak_index), peak_I);
printf("Expected hospital arrivals on that day: lambda = %.2f\n", ...
       lambda_daily(peak_index));
