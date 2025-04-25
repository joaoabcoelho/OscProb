#include "TF1.h"
#include "TGraph.h"

#include "PremModel.h"
#include "PMNS_Fast.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"

#include "Oscillogram.C"



void TimeMeasurement()
{

    auto g = new TGraph();

    for ( int i = 20 ; i<=200 ; i+= 20 ){

        cout << "========== Fast ==========  " << i << endl;
        vector<double> timeFast = Oscillogram(1,i,100,"fast");

        cout << "========== Taylor ==========" << endl;
        vector<double> timeTaylor = Oscillogram(1,i,100,"taylor");

        g->AddPoint( i , timeFast[0] / timeTaylor[0]);

    }

    g->Draw();

    cout<<"MMMM"<<endl;
}