#ifndef SETTINGS_H
#define SETTINGS_H


scalar one(1.);  // stwórz obiekt klasy skalar o nazwie one zainicjalizowany 1
dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);   //  aby nie dzielić przez zero zdefiniowana wartość "small" z przechowywana informacja o wymiarze i warosci SMALL (zdefiniowana w OF)

//bool explicitSolver = runTime.controlDict().lookupOrDefault("explicit", true);
//bool mapFunLog      = runTime.controlDict().lookupOrDefault("mapFunLog", true);
bool limitFieldT    = runTime.controlDict().lookupOrDefault("limitFieldT", false);

double gamma = runTime.controlDict().lookupOrDefault("gamma", 1e-6);
double eps   = runTime.controlDict().lookupOrDefault("eps"  , 5e-16);
double pDim  = runTime.controlDict().lookupOrDefault("pDim" , 1.);
double par   = runTime.controlDict().lookupOrDefault("par"  , 2.);
double gradPsiLimit= runTime.controlDict().lookupOrDefault("gradPsiLimit", 18.); // gradPsinBound

int Ntau   = runTime.controlDict().lookupOrDefault("Ntau", 1);
double parEpsH   = runTime.controlDict().lookupOrDefault("parEpsH", 1.);
double pardTau   = runTime.controlDict().lookupOrDefault("pardTau", 1.);
double ilePkt = runTime.controlDict().lookupOrDefault("ilePkt", 64.);
par = parEpsH;

#endif // SETTINGS_H
