t = linspace(0, 200, 2001) + 0.1;
x0 = [0.99; 0.01; 0.00];

% Solve the SIR system
x = lsode("sir_eqn", x0, t);

% Save results
out = [t', x];
save -ascii sir_octave.out out;

% Plot S, I, R
figure(1);
plot(t, x(:,1), "-r", t, x(:,2), "-g", t, x(:,3), "-b");
xlim([0 200]);
xlabel("Time", "fontweight", "bold");
ylabel("Population", "fontweight", "bold");
legend("S", "I", "R");
grid on;

% Save the current figure so the plot is available in headless runs too.
print("sir_plot.png", "-dpng");
