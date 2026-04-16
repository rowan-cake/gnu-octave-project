% Simulate the SIR model over 200 time units with fine time resolution.
t = linspace(0, 200, 2001) + 0.1;

% Initial state: most of the population is susceptible, with a small
% infected seed and nobody recovered yet.
x0 = [0.99; 0.01; 0.00];

% Solve the SIR system using the seasonal beta(t) defined in sir_eqn.m.
x = lsode("sir_eqn", x0, t);

% Save the time series to disk as: time, susceptible, infected, recovered.
out = [t', x];
save -ascii sir_octave.out out;

% Plot the three compartments so the seasonal infection waves are visible.
figure(1);
plot(t, x(:,1), "-r", t, x(:,2), "-g", t, x(:,3), "-b");
xlim([0 200]);
xlabel("Time", "fontweight", "bold");
ylabel("Population", "fontweight", "bold");
legend("S", "I", "R");
grid on;

% Save the current figure so the plot is available in headless runs too.
print("sir_plot.png", "-dpng");
