# MTSA
Metabolic Thermodynamic Sensitivity Analysis (MTSA) is a novel computational framework designed to explore the thermodynamic and metabolic properties of complex biological systems. By integrating thermodynamic principles with constraint-based modeling, MTSA identifies key reactions and metabolites that drive metabolic fluxes under varying physiological conditions. The method begins with the estimation of physiologically relevant metabolite concentrations using Max/Min Driving Force (MDF) analysis, which assumes that all reactions operate at maximum thermodynamic driving force. This involves parsing reaction formulas, identifying reactants and products, and optimizing metabolite concentrations within physiologically relevant ranges (1 μM to 10 mM). The MDF analysis maximizes the minimum driving force across all reactions while ensuring thermodynamic feasibility based on Gibbs free energy values (ΔG°) and temperature-dependent terms. The resulting metabolite concentrations provide reference points for subsequent kinetic parameter estimation, allowing for a deeper understanding of metabolic network behavior.
Thermodynamic Analysis and Metabolite Concentration Estimation
The MTSA method begins with determining physiologically relevant metabolite concentrations using the Max/Min Driving
Force (MDF) analysis. We assume that each reaction operates at maximum driving force, as all reactions are in pseudosteady state and occur rapidly. So, for each reaction in our metabolic network, we first parse the reaction formula to identify
reactants and products, maintaining stoichiometric coefficients. The MDF analysis employs linear programming optimization
to maximize the thermodynamic driving force while satisfying concentration constraints. Metabolite concentrations are
constrained between 1 μM and 10 mM, reflecting physiologically relevant ranges. The optimization objective maximizes
the minimum driving force across all reactions while ensuring thermodynamic feasibility based on standard Gibbs free
energy values (ΔG°) and the calculated RT term for each temperature point. The MDF analysis yields optimal metabolite
concentrations that maximize the thermodynamic driving force for each reaction. These concentrations serve as reference
points for subsequent kinetic parameter estimation.
Kinetic Parameter Estimation
Following the MDF analysis, we estimate kinetic parameters using the Michaelis-Menten equation. The parameter estimation
process involves multiple steps to ensure robust and physiologically relevant results:
Initial Parameter Range Determination
For each reaction, we establish initial parameter ranges based on Flux Variability Analysis (FVA) results. The maximum flux
determined by FVA serves as a reference point for Vmax estimation. We consider a range of potential Vmax values spanning
from 50% to 150% of this maximum flux, generating five equally spaced estimates within this range. This approach ensures
that our parameter search space encompasses physiologically relevant values while accounting for potential variations in
enzyme activity.
Synthetic Data Generation
Using the optimal concentrations from MDF analysis, we generate a range of substrate concentrations spanning an order
of magnitude below and above the optimal concentration, with 20 logarithmically spaced points. For each substrate
concentration, we generate synthetic reaction rate data using the Michaelis-Menten equation, where Vmax is set to the FVA
upper bound, and Km is estimated as half the optimal substrate concentration. This estimation of Km as half the optimal
concentration is chosen because in the Michaelis-Menten kinetics, Km represents the substrate concentration at which the
reaction rate reaches half of its maximum value, making it a reasonable initial approximation for generating synthetic data
that reflects typical enzyme kinetic behavior. To simulate experimental variation, we add 30% Gaussian noise to the calculated
rates.
Parameter Fitting Process
The Michaelis-Menten equation describes the relationship
between reaction rate (v) and substrate concentration ([S]):
v = (Vmax[S])/(Km + [S])
where:
v is the reaction rate, Vmax is the maximum reaction
rate, [S] is the substrate concentration and Km is the
Michaelis constant, representing the substrate concentration
at which the reaction rate is half of Vmax (see Fig A in supplemetary file of the published paper).
The fitting process employs non-linear regression to
fit the Michaelis-Menten equation to the synthetic data. 
Multiple initial Vmax estimates are tested to avoid local minima in the optimization landscape. For each initial estimate, both
Vmax and Km values are determined while enforcing non-negative constraints. The quality of each fit is assessed using the
R-squared statistic, and the parameter set yielding the highest R-squared value is selected as the final estimate.
Rate-Limiting Substrate Identification
For reactions involving multiple substrates, we identify the rate-limiting substrate based on the kcat/Km ratio. The substrate
exhibiting the lowest kcat/Km ratio is designated as rate-limiting, and its kinetic parameters are used for subsequent
calculations. The kcat values are obtained using the DL-Kcat method, which provides enzyme-specific turnover number
estimates.
Temperature Dependence Integration
The temperature dependence of metabolic reactions is incorporated through several mechanisms:
1. Temperature effects on reaction thermodynamics are captured through the RT term in the MDF analysis
2. Metabolite concentrations and kinetic parameters are determined separately for each temperature point (36-40 ℃)
Quality Control and Validation
The MTSA method includes several quality control measures to ensure reliable results:
1. Validation of parameter bounds to ensure physically meaningful estimates (positive Vmax and Km values)
2. Assessment of fit quality through R-squared statistics
3. Identification and logging of failed parameter estimations
4. Consistency checks between FVA results and estimated kinetic parameters
This comprehensive approach enables systematic investigation of temperature effects on metabolic reaction rates while
maintaining consistency with thermodynamic constraints and network stoichiometry. The method accounts for both
the complexity of enzyme-catalyzed reactions and the limitations of available kinetic data, providing robust estimates of
temperature-dependent metabolic behavior
For further details on the MTSA methodology and its implementation, please refer to the supplementary files provided in this paper.
