#include "TF1.h"

#include "PremModel.h"
#include "PMNS_Fast.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"

#include "Oscillogram.C"



void TimeMeasurement()
{
    Oscillogram(1,2000,1000,"taylor");
}