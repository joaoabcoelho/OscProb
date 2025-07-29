#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"


#include "PremModel.h"
#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"

#include "Oscillogram.C"




void TimeMeasurement()
{

    auto r = new TGraph();


    OscProb::PMNS_Fast f;
    OscProb::PMNS_TaylorExp t;

    OscProb::PremModel prem;    

    double cosT = -0.9;
    prem.FillPath(cosT);
    f.SetPath(prem.GetNuPath());
    t.SetPath(prem.GetNuPath());

    TimeIt time;

    for (int i = 0 ; i<1000 ; i++)
    {
        time.count ++;
    }

    time.Print();
    double timeTotal = time.time();
    double timeIter = timeTotal / time.count;








































    /*
    
    //auto r = new TGraph();
    //auto f = new TGraph();
    //auto t = new TGraph();
    
    double aaa = 0;
    int ii = 0 ;

    cout << "========== Taylor ==========" << endl;
    vector<double> timeTaylor = Oscillogram(1,200,100,"taylor");



    //---------------------------------------------------------------------------
    for ( int i = 0 ; i<=500 ; i+= 1 ){

        cout << "========== Fast ==========  " << i << endl;
        //vector<double> timeFast = Oscillogram(1,i,100,"fast");

        cout << "========== Taylor ==========" << endl;
        vector<double> timeTaylor = Oscillogram(1,200,100,"taylor");

        //r->AddPoint( i , timeFast[0] / timeTaylor[0]);
        //f->AddPoint( i , timeFast[1] );
        t->AddPoint( i , timeTaylor[1] );

        cout<<endl;

        aaa += timeTaylor[0];
        ii++;
    }

    cout<<"temps moyen :" << aaa/ii<<endl;  
    
    TFile file1("graphTime.root", "RECREATE");

    r->SetName("ratio");
    r->SetTitle("ratio;nbrBinE;ratioTemps");
    r->Draw("AC*");
    r->Write();

    f->SetName("temps/iter fast");
    f->SetTitle("temps/iter fast;nbrBinE;temps/iter fast");
    f->Draw("AC*");
    f->Write();

    t->SetName("time/iter_taylor");
    t->SetTitle("time/iter_taylor;nbrBinE;time/iter_taylor");
    t->Draw("AC*");
    t->Write();

    file1.Close();
    
    */
    //---------------------------------------------------------------------------

    /*do{
        cout << "========== Fast ==========  " << i << endl;
        //vector<double> timeFast = Oscillogram(1,i,100,"fast");

        cout << "========== Taylor ==========" << endl;
        vector<double> timeTaylor = Oscillogram(1,200,100,"taylor");

    }while();*/
    //---------------------------------------------------------------------------

    // SEARCH WHEN FAST FASTER THAN TAYLOR
    /*vector<double> timeFast;
    vector<double> timeTaylor;
    int i = 20;
    do{
        timeFast.clear();
        timeTaylor.clear();

        cout << "========== Fast ==========  " << i << endl;
        timeFast = Oscillogram(1,i,100,"fast");

        cout << "========== Taylor ==========" << endl;
        timeTaylor = Oscillogram(1,i,100,"taylor");

        r->AddPoint( i , timeFast[0] / timeTaylor[0]);

        i += 10;

        cout<<endl;

    }while(timeFast[0] > timeTaylor[0]);*/
    //---------------------------------------------------------------------------
}