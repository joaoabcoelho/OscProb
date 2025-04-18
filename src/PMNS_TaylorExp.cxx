///////////////////////////////////////////////////////////////////////////////
//info
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "PMNS_TaylorExp.h"

using namespace OscProb;

using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_TaylorExp::PMNS_TaylorExp() : PMNS_Fast() 
{
    InitializeTaylorsVectors();

    SetwidthBin(0.1);
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_TaylorExp::~PMNS_TaylorExp() {}

//.............................................................................
///
///
///
void PMNS_TaylorExp::InitializeTaylorsVectors()
{
    fKE = matrixC(fNumNus, vectorC(fNumNus, 0));

    flambdaE = vectorD(fNumNus, 0);

    fVE = matrixC(fNumNus, vectorC(fNumNus, 0));

    fevolutionMatrixS = matrixC(fNumNus, vectorC(fNumNus, 0));
    for(int i= 0 ; i<fevolutionMatrixS.size(); i++){
        fevolutionMatrixS[i][i] = 1;
    }


    //REEL OU COMPLEX??????????
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::SetwidthBin(double widthBin)
{
    fwidthBin = widthBin;
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::rotateS(vectorC fPhases,matrixC& S)
{
    for(int j = 0 ; j<fNumNus ; j++) //PEUT ETRE AMELIORER
    {
        for(int i = 0 ; i<=j ; i++)
        {
            for(int k = 0 ; k<fNumNus ; k++)
            {
                S[i][j] += fEvec[i][k] * fPhases[k] * conj(fEvec[j][k]);
            }

            if(i != j){
                S[j][i] = -conj(S[i][j]);//??????????????????????
            }
        }
    }
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::rotateK(matrixC Kmass,matrixC& Kflavor)
{
    for(int j = 0 ; j<fNumNus ; j++)
    {
        for(int i = 0 ; i<=j ; i++)
        {
            for(int k = 0 ; k<fNumNus ; k++)
            {
                for(int l = 0 ; l<fNumNus ; l++)
                {
                    Kflavor[i][j] += fEvec[i][k] * Kmass[k][l] * conj(fEvec[j][l]);
                }
            }

            if(i != j){
                Kflavor[j][i] = conj(Kflavor[i][j]);
            }
        }
    }
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::BuildKE(double L , matrixC& K)
{
    for(int j = 0 ; j<fNumNus ; j++)
    {
        for(int i = 0 ; i<=j ; i++)
        {
            for(int k = 0 ; k<fNumNus ; k++)
            {
                for(int l = 0 ; l<fNumNus ; l++)
                {
                    K[i][j] += conj(fEvec[k][i]) * fHam[k][l] * fEvec[l][j];
                }
            }

            complexD C;
            if(i == j){
                C = complexD(1,0);
            }
            else{
                double argg = (fEval[i] - fEval[j]) * L;
                C = (complexD(cos(argg), sin(argg)) - complexD(1,0) ) / (complexD(0,argg)); //COMPLEX C???????
                cout<<C<<endl;
            }

            K[i][j] *= (-L / (2*fEnergy)) * K[i][j] * C; //UNIT2?????????

            if(i != j){
                K[j][i] = conj(K[i][j]);
            }

        }
    }
}

//.............................................................................
///
///
///
void PMNS_TaylorExp::MultiplicationRule(matrixC SLayer,matrixC KLayer)
{
    //K (ne pas modifier S avant de l'appliquer ici )

    //S
    for(int j = 0 ; j<fNumNus ; j++)
    {
        for(int i = 0 ; i<=j ; i++)
        {
            for(int k = 0 ; k<fNumNus ; k++)
            {
                if(j>k){
                    fevolutionMatrixS[i][j] -= SLayer[i][k] * conj(fevolutionMatrixS[j][k]);
                }
                else{
                    fevolutionMatrixS[i][j] += SLayer[i][k] * fevolutionMatrixS[k][j];
                }
            }
        }
    }
    //PEUT ETRE AMELIORER
    for(int j = 1 ; j<fNumNus ; j++)
    {
        for(int i = 0 ; i<j ; i++)
        {
            fevolutionMatrixS[j][i] = -conj(fevolutionMatrixS[i][j]);
        }
    }

}

//.............................................................................
///
///
///
double PMNS_TaylorExp::avrProbTaylor(double E , double widthBin)
{
    // reset K et S et Ve et lambdaE
    InitializeTaylorsVectors();

    SetEnergy(E);
    SetwidthBin(widthBin);

    //Propagate -> get S and K matrix (on the whole path)
    PropagateTaylor();

    //DiagolK -> get VE and lambdaE

    //return fct avr proba

    return 1;
}

//.............................................................................
///
/// Propagate neutrino state through full path
///
void PMNS_TaylorExp::PropagateTaylor()
{
  for (int i = 0; i < int(fNuPaths.size()); i++) { PropagatePathTaylor(fNuPaths[i]); }
}


//.............................................................................
///
/// Propagate the current neutrino state through a given path
/// @param p - A neutrino path segment
///
void PMNS_TaylorExp::PropagatePathTaylor(NuPath p)
{
    // Set the neutrino path
    SetCurPath(p);
    double L = p.length;

    // Solve for eigensystem
    SolveHam();

    // Get the evolution matrix in mass basis
    double LengthIneV = kKm2eV * p.length;
    for (int i = 0; i < fNumNus; i++) {
        double arg = fEval[i] * LengthIneV;
        fPhases[i] = complexD(cos(arg), -sin(arg));
    }

    // Rotate S in flavor basis
    matrixC Sflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
    rotateS(fPhases,Sflavor);

    // Build KE in mass basis
    matrixC Kmass = matrixC(fNumNus, vectorC(fNumNus, 0));
    BuildKE(L,Kmass);

    // Rotate KE in flavor basis
    matrixC Kflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
    rotateK(Kmass,Kflavor);

    //multiplication rule for K and S -> uptade K and S
    for(int j=0 ; j<3 ; j++)
    {
        for(int k=0 ; k<3 ; k++)
        {
            cout<<Kflavor[j][k]<<" "; //--> pas de probleme, cela vient donc de MultiplicationRule?? enfin si un prb car valeur trop faible summ!=1
        }
        cout<<endl;
    }


    MultiplicationRule(Sflavor,Kflavor);


    cout<<endl;


    for(int j=0 ; j<3 ; j++)
    {
        for(int k=0 ; k<3 ; k++)
        {
            cout<<fevolutionMatrixS[j][k]<<" ";
        }
        cout<<endl;
    }
 
}


