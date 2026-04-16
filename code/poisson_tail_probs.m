function probs = poisson_tail_probs(lambda, max_bucket)
  % Return probabilities for arrival counts:
  %   0, 1, 2, ..., max_bucket - 1, and a final tail bucket for >= max_bucket
  %
  % This keeps the dynamic program finite while still preserving total
  % probability mass. The tail bucket maps to the capped census n = Nmax.

  probs = zeros(max_bucket + 1, 1);

  % Start the Poisson recurrence at P(X = 0) = exp(-lambda).
  probs(1) = exp(-lambda);

  % Build exact probabilities for 1 through max_bucket - 1 using the
  % stable recurrence P(m) = P(m - 1) * lambda / m.
  for m = 1:(max_bucket - 1)
    probs(m + 1) = probs(m) * lambda / m;
  endfor

  % Capture the entire remaining tail in the final bucket.
  probs(end) = max(0, 1 - sum(probs(1:max_bucket)));
endfunction
