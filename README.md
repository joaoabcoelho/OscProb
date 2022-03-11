# OscProb

OscProb is a small set of classes aimed at computing exact neutrino oscillation probabilities with a few different models.

OscProb contains a basic framework for computing neutrino oscillation probabilities. It is integrated into [ROOT](https://root.cern.ch/), so that each class can be used as you would any ROOT class.

Available classes are:
- **PremModel:** Used for determining neutrino paths through the earth
- **PMNS_Fast:** Standard 3-flavour oscillations
- **PMNS_Iter:** Standard 3-flavour oscillations (iterative)
- **PMNS_Sterile:** Oscillations with any number of neutrinos
- **PMNS_NSI:** Oscillations with 3 flavours including Non-Standard Interactions
- **PMNS_Deco:** Oscillations with 3 flavours including a simple decoherence model
- **PMNS_LIV:** Oscillations with 3 flavours including Lorentz Invariance Violations
- **PMNS_Decay:** Oscillations with 3 flavours including neutrino decays of the second and third neutrino mass states \nu_2 and \nu_3. [Requires external library Eigen3, see the instructions below.]
- **Absorption:** Computes absorption probabilities for high-energy neutrinos

A few example macros on how to use OscProb are available in a tutorial directory.

# Installing OscProb

OscProb is very easy to install. The only requirements is to have ROOT installed with the GSL libraries.

**NEW: Thanks to Jacek Holeczek, OscProb now also builds with ROOT 6!!**

In order to compile the PMNS_Decay class, it is necessary to donwload the external Eigen library. This library is added as a submodule. There are two options:
- During cloning: `git clone --recurse-submodules https://github.com/joaoabcoelho/OscProb.git`
- After clonning: `git submodule update --init`

Once you have ROOT setup, simply do:
```sh
cd OscProb
make
```

A shared library will be produced: ```libOscProb.so```

This should take a few seconds and you are all set.

Just load the shared library in your ROOT macros with:
```cpp
gSystem->Load("/full/path/to/libOscProb.so");
```

Or use the ```LoadOscProb.C``` macro (see below).

# Tutorial

In the directory OscProb/tutorial you will find a few macros with examples using OscProb.

Two macros are particularly useful:
- ```simpleExamples.C``` : Contains some short pieces of code on how to perform different tasks.
- ```MakeOscillogram.C``` : Runs a full example of how to plot an oscillogram with the PREM model.

Additionally, these macros contain useful tools:
- ```LoadOscProb.C```: Searches for the OscProb library in your current directory, parent directory, or library path, and then loads it. It is called within the tutorial macros as a possible usage example.
- ```SetNiceStyle.C```: Provides simple tools to make your plots look nicer. Feel free to use it anytime you're making plots, even if you're not running OscProb. This is completely independent of OscProb.

To run macros in compiled mode you will need to preload the OscProb library, e.g:

```sh
root -l LoadOscProb.C MakeOscillogram.C+
```

# Documentation

More detailed documentation of the code can be found in a Doxygen page here:

[Doxygen](https://joaoabcoelho.github.io/OscProb/doxygen/html/annotated.html "OscProb Doxygen page")
