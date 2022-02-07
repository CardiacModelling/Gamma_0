# Notes on main papers for Yann's Gamma-0 work

- [x] Attwell, Cohen, Eisner (1979) Membrane potential and ion concentration stability conditions for a cell with a restricted extracellular space
  - Examines steady state, assuming Vr is a fixed point.

- [x] Varghese, Sell (1997) A conservation principle and its effect on the formulation of Na–Ca exchanger current in cardiac cells
  - Examine SAN model Noble 89, Winslow 91.
  - "derivation of a conservation principle relating V and 5 concentrations"
  - "this conservation priniciple is not derived physical considerations; it is a consequence of the structure of the equations"
  - but "implications for physical interpretation .. presented in the Discussion"
  - divide model into voltage (`u`), gating variables (`v_i`), and concentrations (`w_i`)
  - Show that eq can be rewritten as V = linear sum of beta * concentration, minus a constant u0. And that -u0 is a lower bound for V (given that concentrations are positive).
  - Then give a few examples of V = ... + C0 for different models.
  - "In this paper we have shown the presence of a latent conservation principle in equations describing the electrochemical activity in cardiac cells"
  - "The conservation principle can be thought of as a reformulation of the familiar principle of Faraday: Q = CV, rewritten in terms of ionic concentrations rather than the charge."

- [x] Guan, Lu, Huang (1997) A discussion about the DiFrancesco-Noble model
  - Look at DN1985 model with channels, pumps, concentrations, V
  - Have read Varghese (in 1994 preprint form), but have independently arrived at something similar
  - Look at fixed point d/dx of all states = 0, show that there are five eqs but only 4 variables (C0 is missing)
  - Limit cycles are "not isolated": if a stable limit cycle has been obtained, a new stable limit cycle can be obtained by perturbing the value of Nai kr Ki from the original one (because C0 mis missing).
  - Say that Nai and Ki take ages to converge, and "therefore" have little effect on the other vars and can be treated as parameters
  - Claim this doesn't affect its "main dynamics"
  - In eq 1 miss out C0

- [x] Endresen, Hall, Hoye, Myrheim (2000) A theory for the membrane potential of living cells
  - Probably didn't know about Varghese until review stage, give it only a very brief mention in discussion.
  - "We give an explicit formula for Vm in terms of the ... concentrations"
  - We demonstrate that the work done by the pumps equals the change in potential energy of the cell, plus the energy lost in downhill ionic fluxes through the channels and exchangers."
  - Bold claims: "The model predicts the experimentally observed intracellular ionic concentration of potassium, calcium and sodium."
  - "We do not see any drift in the values of the concentrations in a long time simulation"
  - "and we obtain the same asymptotic values when starting from the full equilibrium situation with equal intracellular and extracellular ionic concentrations"
  - "we try to make the theory realistic by using equations that are compatible with, or can be derived from, basic physical principles."
  - "We observe that this differential equation can be integrated exactly, and argue that the integration constant is given by the requirement that the potential is zero when the ion concentrations on both sides of the membrane are equal, as the density of negative charge happens to be the same on both sides."
  - "a formula which is nothing but the one for an electric capacitance that follows from Gauss's law in electrostatics."
  - Using:
    - Boltzmann distribution: particle in thermal eq spends more time in lower energy states than in higher
    - Markov assumption: transition of stoch sys depends only on present state
    - Detailed balance: microscopic laws of physics are invariant w.r.t. direction of time
  - Inspired by
    - Ehrenstein and Lecar's model of channel gating (1977)
    - Nonner and Eisenberg's model for channel current (1998)
    - Mullins' model of the NaCa exchanger (1977), and
    - Chapman's model of the NaK pump (1978).
    - In particular the book by Mullins (1981) ``Ion transport in heart''
  - "By using an algebraic equation for the potential in place of the standard differential equation, as mentioned above, we obtain a model which is stable against a slow drift of the intracellular ion concentrations, sometimes seen in other models."
  - "Furthermore, by fixing the integration constant for the voltage we obtain from the model a prediction of the steady state ion concentrations in the cell."
  - "It is even possible to predict these steady state concentrations by starting with an initial state having equal concentrations inside and outside the cell, and
integrating the equations of motion over a long time interval."
  - Come up with integration constant v0
  - But then assume only K, Ca, and Na set membrane potential, so give expression for v0 in terms of external concentrations
  - Claim that "the constant v0 in Eq. (62) can be uncertain since the inside and the outside of the cell ... is different, but this can be neglected in [equation for v0] as anything like +/-100 mV has a tiny influence upon concentrations as the prefactor FV/C is very large. This we also verified by a simulation."
  
- [x] Hund, Kucera, Otani, Rudy (2001) Ionic Charge Conservation and Long-Term Steady State in the Luo-Rudy Dynamic Cell Model
  - "When ions carried by the stimulus current are taken into account, the algebraic and differential methods yield identical results and neither shows drift in computed parameters."
  - "The present study establishes the proper pacing protocol for simulation studies of cellular behavior during long periods of rapid pacing."

- [x] Genet, Costalat, Burger (2001) The influence of plasma membrane electrostatic properties on the stability of cell ionic composition
  - Looks at osmosic and Na+ relationships

- [x] Kneller, Ramirez et al., Nattel (2002) Time-dependent transients in an ionically based mathematical model of the canine atrial action potential
  - Add Cl concentration bookkeeping, after a previous model (Ramirez, Nattel, Courtemanche, 2000) of canine atrial cells included IClCa without accounting for Cl.

- [x] Jacquemet (2007) Steady-state solutions in mathematical models of atrial cell electrophysiology and their stability
 - Determines fixed points in Courtemanche and (modified) Nygren models
 - "Symbolic calculations were carried out as far as possible in order to prove the existence of these fixed points."
 - "In the Fenton–Karma model, a unique stable fixed point was found, namely the resting state."
 - "The Courtemanche model had an infinite number of fixed points. A bifurcation diagram was constructed by classifying these fixed points according to a conservation law."
 - "The Nygren model had ... an infinite number of fixed points, resulting in a bifurcation diagram similar to that of the Courtemanche model."
 - "...it is natural to start the analysis of such systems by identifying the fixed points... If these ... are stable, they correspond to resting states, whereas unstable fixed points may potentially induce pacemaker activity."
 - "the stability of the resting potential is often postulated in theoretical studies [4,5]" (Attwell, Genet)
 - To find the fixed state, they divide into V, concentrations c, gating variables y.
   - Then say fixed values for y are given by `y_inf(V_fixed, c_fixed)`
   - Then say V is a function of c and Q0, `V_fixed = f(c_fixed, Q0)`, so that there is one fixed V for any value of Q0 (Section 2.3)
   - Small perturbs: jacobian, large perturbs: AP simulation
 - Derive `C * V = F(concentrations) - Qns` where `Qns` is "non-specific charge", and is positive if Vr is negative.
   - `Qns = Vi F [NS]i`
   - Say Endresen assume Qns approx 0 at V=0,
   - while Hund et al., calculate it from chosen initial conditions
 - Explain Nygren model contained `phi_Na,en` parameter to account for influx from stimulus, but that this causes `dot(Qns) = -phi_Na,en` so that no fixed points exist (because dot(Qns) is never zero).
   - Which relates to Hund paper's point that stimulus needs to be accounted for in the concentrations.
  - Show that in Courtemanche there stable and unstable fixed points, and plot them as a function of `[NS]i`
    - Find a very different result in Nygren: two stable Vs (for realistic `[NS]i`), indicating a low and a high Vr are possible!
    - Says this occurs in Nygren but not in Courtemanche as a result of the different IK1 formulations!

- [x] Fraser, Huang (2007) Quantitative techniques for steady-state calculation and dynamic integrated modeling of membrane potential and intracellular ion concentrations
  - Argues for "charge-difference" modelling instead of "current-sum" models
  - Says stimulus addition will cause drift in models if not accounted for in concentrations

- [x] Livshitz, Rudy (2009) Uniqueness and Stability of Action Potential Models during Rest, Pacing, and Conduction Using Problem-Solving Environment
  - Rate-dependent long-term ion accumulation occurs in various disease processes, requires long-running simulations (several minutes)
  - Reaching steady-state is important before comparing normal and "diseased" model
  - People have "raised issues" including
    - "an apparent dependence of the solution on initial conditions"
    - "drift of the state variables (concentrations)"
    - discontinuities in kinetics (INa if statement)
  - Say stability issues are resolved by accepting it is a DAE; an ODE system with _contstraints_
    - say that "a consistent set of initial values for state variables and their derivatives must be supplied to insure stability and uniqueness of the solution"
    - (Not clear why only applies for DAE?)
    - constraint: conservation V and concs
    - constraint: conservation state occupancy in markov models
  - Say that "monphasic" stimulus carried by K+ was used before to not break rule 1
    - But "biphasic stimulation" occurs in fibers, by which I think they mean a negative current from the neighbours, followed by a positive current when this cell acts as a source
  - Five steps in this paper:
    - "uniqueness of solutions with a consistent set of initial conditions is tested for both quiescent and periodically paced modes"
    - "a conservative biphasic stimulation protocol for a cell is formulated"
    - "a convergence criterion for the steady-state periodic solution is defined, based on the incremental contribution of each ion species to the transmembrane potential during each beat."
    - "analytical algebraic expressions are derived to remove singularities in the state variables formulations.
    - "a MATLAB script is developed and implemented for single-cell and cardiac strand simulations"
  - Write Vm eq as `CmVm = -Q0 + Qstim + F sum_x sum_k {v_k z_x, [x]_{x,k}}`, where the subscript `k` is to take compartments into account (with volume `v_k`).
  - Says they use Matlab's ODE15s for DAEs. ODE15s does seem to allow ODEs, by specifying a "mass matrix" by which y' will be pre-multiplied. If you add zeros in this, you get `0 = some eq` in your RHS, which allows constraints to be added.
  - To set up simulations, they say to set all gates/markov models to a likely state, then set concentrations and V to what you think they should be, then let model relax - without pacing - into its "quiescent state". Then start from there.
    - They also call the "quiescent state" the "autonomic regime"
  - Show that pacing with a single pulse and adding it to K+ leads to different results than using a "biphasic beat" (quick pulse inwards, followed by long small pulse outward, net charge = 0), in terms of ion concentration bookkeeping. Unsure if they ascribe the biphasic one to Ki too.
    - Say that using this biphasic pulse is a good way to find initial conditions for multi-cell simulations (which act as sink, then as source, so also biphasic).
    - Say that otherwise, esp. at fast rates, you get bad estimates for what will happen to Ki and Nai ("bad" compared to biphasic stim).
  - Define stopping criteria for pacing that says "stop if contribution of change in ionic species X to change in Vm is less than some threshold".

- [x] Van Oosterom, Jacquemet (2009) Ensuring stability in models of atrial kinetics
  - Looks at quiescent steady state
  - Then introduces feedback mechanism to "anchor" Nai and Ki. Not great!



