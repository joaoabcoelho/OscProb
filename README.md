# OscProb

OscProb is a small set of classes aimed at computing exact neutrino oscillation probabilities with a few different models.

OscProb contains a basic framework for computing neutrino oscillation probabilities. It is integrated into [ROOT](https://root.cern.ch/), so that each class can be used as you would any ROOT class.

Available classes are:
- **PremModel:** Used for determining neutrino paths through the earth
- **PMNS_Fast:** Standard 3-flavour oscillations
- **PMNS_Sterile:** Oscillations with any number of neutrinos
- **PMNS_NSI:** Oscillations with 3 flavours including Non-Standard Interactions
- **PMNS_Deco:** Oscillations with 3 flavours including a simple decoherence model

A few example macros on how to use OscProb are available in a tutorial directory.

# Installing OscProb

OscProb is very easy to install. The only requirements is to have ROOT installed with the GSL libraries.

**NEW: Thanks to Jacek Holeczek, OscProb now also builds with ROOT 6!!**

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

In some cases you may need to explicitly load GSL libraries. Just add something like this to your rootlogon, replacing the path to the GSL libraries if needed:
```cpp
gSystem->Load("/usr/lib/x86_64-linux-gnu/libgsl.so");
gSystem->Load("/usr/lib/x86_64-linux-gnu/libgslcblas.so");
```

# Tutorial

In the directory OscProb/tutorial you will find a few macros with examples using OscProb.

Two macros are particularly useful:
- ```simpleExamples.C``` : Contains some short pieces of code on how to perform different tasks.
- ```MakeOscillogram.C``` : Runs a full example of how to plot an oscillogram with the PREM model.

Additionally, the macro ```SetNiceStyle.C``` will provide simple tools to make your plots look nicer. Feel free to use it anytime you're making plots, even if you're not running OscProb. This is completely independent of OscProb.

# Documentation

More detailed documentation of the code can be found in a Doxygen page here:

[Doxygen](https://joaoabcoelho.github.io/OscProb/doxygen/html/annotated.html "OscProb Doxygen page")
