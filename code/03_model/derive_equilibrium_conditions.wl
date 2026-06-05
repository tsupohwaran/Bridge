(* Derive equilibrium existence and uniqueness checks for the current
   one-level logit monopsony model. Run with:

      wolframscript -file code/03_model/derive_equilibrium_conditions.wl

   If wolframscript cannot find the kernel on this Mac, run:

      /Applications/Wolfram.app/Contents/MacOS/WolframKernel \
        -script code/03_model/derive_equilibrium_conditions.wl

   The production curvature is fixed at alpha = 0.4 = 2/5.
*)

ClearAll["Global`*"];

alpha = 2/5;

Print["=== Model primitives ==="];
Print["alpha = ", alpha];
Print["Existence sufficient conditions:"];
Print["  theta > 0; eta finite (eta >= 0 for commuting disutility);"];
Print["  d_zj > 0; l_z > 0 with Sum_z l_z = 1; z_j > 0; finite a_j; J >= 2."];
Print[""];

Print["=== Softmax derivative check ==="];
Clear[th, x1, x2, x3];
den = Exp[th*x1] + Exp[th*x2] + Exp[th*x3];
p1 = Exp[th*x1]/den;
p2 = Exp[th*x2]/den;
p3 = Exp[th*x3]/den;

softmaxChecks = {
  FullSimplify[D[p1, x1] == th*p1*(1 - p1), Assumptions -> th > 0],
  FullSimplify[D[p1, x2] == -th*p1*p2, Assumptions -> th > 0],
  FullSimplify[D[p1, x3] == -th*p1*p3, Assumptions -> th > 0]
};
Print["D pi_j / D x_k = theta pi_j (1{j=k} - pi_k): ", softmaxChecks];
Print[""];

Print["=== General log-wage map condition ==="];
Print["Let x_j = log w_j and define"];
Print["  pi_zj = softmax_j(theta*(x_j + a_j - eta*log d_zj)),"];
Print["  L_j = Sum_z l_z pi_zj,"];
Print["  gamma_zj = l_z pi_zj / L_j,"];
Print["  S_jk = Sum_z gamma_zj pi_zk,"];
Print["  R_jk = Sum_z gamma_zj pi_zj pi_zk,"];
Print["  eps_j = theta*(1 - S_jj)."];
Print["Then Mathematica's softmax derivative implies"];
Print["  D log L_j / D x_k = theta*(1{j=k} - S_jk)."];
Print["It also implies"];
Print["  D S_jj / D x_k = theta*(1{j=k} S_jj - 2 R_jk + S_jj S_jk)."];
Print["For the unnormalized wage map"];
Print["  F_j(x) = log(alpha z_j) + (alpha-1) log L_j + log eps_j - log(1+eps_j),"];
Print["the Jacobian is"];
Print["  D F_j / D x_k = (alpha-1)*theta*(1{j=k} - S_jk)"];
Print["    - theta*(1{j=k} S_jj - 2 R_jk + S_jj S_jk)"];
Print["      / ((1 - S_jj)*(1 + theta*(1 - S_jj)))."];
Print["With alpha = 0.4, alpha - 1 = ", alpha - 1, "."];
Print["A sufficient global uniqueness condition is a contraction:"];
Print["  sup_x max_j Sum_k Abs[D normalized(F)_j / D x_k] < 1."];
Print["A practical local uniqueness/stability check at an equilibrium x* is:"];
Print["  SpectralRadius[D normalized(F)(x*)] < 1."];
Print["eta enters these conditions through pi_zj, hence through S_jk and eps_j."];
Print["Therefore the general model has no data-free cutoff using only theta and eta."];
Print["For a fixed eta and dataset, theta must be small enough that the above"];
Print["normalized Jacobian remains a contraction over the relevant wage domain."];
Print[""];

Print["=== Two-firm, one-origin benchmark ==="];
Clear[theta, u, p, c];
p = 1/(1 + Exp[-theta*u]);

(* In the two-firm one-origin case, pi_1=p, pi_2=1-p,
   eps_1=theta*(1-p), eps_2=theta*p. The log wage-ratio map is
   F(u) = c + log w_1 update - log w_2 update. *)
F = c + (alpha - 2)*Log[p/(1 - p)] +
    Log[(1 + theta*p)/(1 + theta*(1 - p))];

Fprime = FullSimplify[D[F, u], Assumptions -> theta > 0 && u \[Element] Reals];
Print["Derivative of two-firm log wage-ratio map:"];
Print["  F'(u) = ", Fprime];

Clear[q];
expr = theta*(alpha - 2) +
   theta^2*q*(1 - q)*(1/(1 + theta*q) + 1/(1 + theta*(1 - q)));

signCheck = FullSimplify[expr < 0, Assumptions -> theta > 0 && 0 < q < 1];
boundCheck = FullSimplify[
   -(2 - alpha)*theta <= expr <= -(1 - alpha)*theta,
   Assumptions -> theta > 0 && 0 < q < 1
];
contractionCondition = FullSimplify[(2 - alpha)*theta < 1, Assumptions -> theta > 0];

Print["For q in (0,1), F'(u)<0: ", signCheck];
Print["Derivative bound: ", boundCheck];
Print["Sufficient global contraction condition:"];
Print["  (2-alpha)*theta < 1  <=>  ", contractionCondition];
Print["With alpha=0.4 this is theta < ", N[1/(2 - alpha)], " = 5/8."];
Print["This benchmark has no direct eta bound; eta shifts q through commute costs."];
