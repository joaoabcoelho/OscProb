#include "TF1.h"

#include "PremModel.h"
#include "PMNS_Fast.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"

#include "Oscillogram.C"



void TimeMeasurement()
{
    cout << "========== Fast ==========" << endl;
    Oscillogram(1,200,100,"fast");

    cout << "========== Taylor ==========" << endl;
    Oscillogram(1,200,100,"taylor");
}