[![DOI](https://zenodo.org/badge/91241668.svg)](https://zenodo.org/badge/latestdoi/91241668)

# OscProb

OscProb is a small set of classes aimed at computing exact neutrino oscillation probabilities with a few different models.

OscProb contains a basic framework for computing neutrino oscillation probabilities.
It is integrated into [ROOT](https://root.cern.ch/), so that each class can be used as you would any ROOT class.

Available classes are:

- **PremModel:** Used for determining neutrino paths through the earth
- **PMNS_Fast:** Standard 3-flavour oscillations
- **PMNS_Iter:** Standard 3-flavour oscillations (iterative approximation)
- **PMNS_Sterile:** Oscillations with any number of neutrinos
- **PMNS_NSI:** Oscillations with 3 flavours including vector Non-Standard Interactions
- **PMNS_SNSI:** Oscillations with 3 flavours including scalar Non-Standard Interactions
- **PMNS_Deco:** Oscillations with 3 flavours including a simple decoherence model
- **PMNS_OQS:** Oscillations with 3 flavours including decoherence from Open Quantum Systems
- **PMNS_LIV:** Oscillations with 3 flavours including Lorentz Invariance Violations
- **PMNS_Decay:** Oscillations with 3 flavours including neutrino decays
- **PMNS_NUNM:** Oscillations with 3 flavours including non-unitary neutrino mixing
- **Absorption:** Computes absorption probabilities for high-energy neutrinos

A few example macros on how to use OscProb are available in a [tutorial](tutorial) directory.

# Installing OscProb

OscProb is very easy to install. The only requirements is to have ROOT installed.

**NEW: Thanks to Jacek Holeczek, OscProb now also builds with ROOT 6!!**

<details>
  <summary>Details: Eigen</summary>

  In order to compile the `PMNS_Decay` and `PMNS_Sterile` classes, it is necessary to
  download the external Eigen library. This library is added as a submodule and will be downloaded by make.

  Alternatively, you can trigger them manually with:

- During cloning: `git clone --recurse-submodules https://github.com/joaoabcoelho/OscProb.git`
- After cloning: `git submodule update --init`

</details>

Once you have ROOT setup, simply do:

```sh
cd OscProb
make
```

A shared library will be produced: ```lib/libOscProb.so```

This should take a few seconds and you are all set.

<details>
  <summary>Compile with cmake</summary>
  CMake is a modern tool to build package software. The main features are out-of-source compilation and the possibility of making the package easily available for dependencies.
  CMake will take care of installing the Eigen3 dependency, so there is not even need to include it as a submodule.
  To quickly install it use:

  ```shell
  build_dir=$PWD/build # Change with the build dir you want
  install_dir=$PWD/install # Change with the install dir you want
  cmake -B $build_dir --install-prefix $install_dir -DCMAKE_BUILD_TYPE=Release -DOSCPROB_ENABLE_TESTING=ON
  cmake --build $build_dir --target install
  ```

  To test the installation you can do:
  
  ```shell
  cd $install_dir/test
  root -b -q ../tutorial/LoadOscProb.C TestMethods.C
  root -b -q ../tutorial/LoadOscProb.C StressTest.C
  ```

</details>

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

- ```LoadOscProb.C```: Searches for the OscProb library in your current directory, parent directory, or library path, and then loads it.
- ```SetNiceStyle.C```: Provides simple tools to make your plots look nicer. Feel free to use it anytime you're making plots, even if you're not running OscProb. This is completely independent of OscProb.

To run macros you will need to preload the OscProb library, e.g:

## Interpreter mode

```sh
root -l LoadOscProb.C MakeOscillogram.C
```

## Compiler mode

```sh
root -l LoadOscProb.C MakeOscillogram.C+
```

# Python

You can use OscProb in Python using PyROOT by simply loading the shared library:

```py
import ROOT
ROOT.gSystem.Load('/path/to/OscProb/libOscProb.so')

p = ROOT.OscProb.PMNS_Fast()

print(p.Prob(1, 1, 5.0))
```

# Documentation

More detailed documentation of the code can be found in a Doxygen page here:

[Doxygen](https://joaoabcoelho.github.io/OscProb/ "OscProb Doxygen page")
