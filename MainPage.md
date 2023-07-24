# OscProb  {#index}

OscProb is a small set of classes aimed at computing exact neutrino oscillation probabilities with a few different models.

OscProb contains a basic framework for computing neutrino oscillation probabilities. It is integrated into [ROOT](https://root.cern.ch/), so that each class can be used as you would any ROOT class.

Available classes are:
- **[PremModel](@ref OscProb::PremModel):** Used for determining neutrino paths through the earth
- **[PMNS_Fast](@ref OscProb::PMNS_Fast):** Standard 3-flavour oscillations
- **[PMNS_Iter](@ref OscProb::PMNS_Iter):** Standard 3-flavour oscillations (iterative)
- **[PMNS_Sterile](@ref OscProb::PMNS_Sterile):** Oscillations with any number of neutrinos
- **[PMNS_NSI](@ref OscProb::PMNS_NSI):** Oscillations with 3 flavours including vector Non-Standard Interactions
- **[PMNS_SNSI](@ref OscProb::PMNS_SNSI):** Oscillations with 3 flavours including scalar Non-Standard Interactions
- **[PMNS_Deco](@ref OscProb::PMNS_Deco):** Oscillations with 3 flavours including a simple decoherence model
- **[PMNS_LIV](@ref OscProb::PMNS_LIV):** Oscillations with 3 flavours including Lorentz Invariance Violations
- **[PMNS_Decay](@ref OscProb::PMNS_Decay):** Oscillations with 3 flavours including neutrino decays of the second and third neutrino mass states nu_2 and nu_3. [Requires external library Eigen3, see the instructions below.]
- **[Absorption](@ref OscProb::Absorption):** Computes absorption probabilities for high-energy neutrinos

A few example macros on how to use OscProb are available in a tutorial directory.
