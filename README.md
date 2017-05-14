# OscProb

OscProb is a small set of classes aimed at computing exact neutrino oscillation probabilities with a few different models.

OscProb contains a basic framework for computing neutrino oscillation probabilities. It is integrated into ROOT, so that each class can be used as you would any ROOT class.

Available classes are:
- PremModel: Used for determining neutrino paths through the earth
- PMNS_Fast: Standard 3-flavour oscillations
- PMNS_Sterile: Oscillations with any number of neutrinos
- PMNS_NSI: Oscillations with 3 flavours including Non-Standard Interactions

A few example macros on how to use OscProb are available in a tutorial directory.
