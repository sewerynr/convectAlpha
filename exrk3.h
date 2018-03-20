#ifndef EXRK4_H
#define EXRK4_H

#include <stdlib.h>
#include "localFunctions.h"

#include "fvCFD.H"
//#include "fvOptions.H"

void exRK3(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau,const volScalarField& PsiZero, volScalarField& Psi,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz, double gamma, const double ilePkt, const double gradPsiLimit);

#endif // EXRK4_H
