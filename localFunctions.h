#ifndef LOCALFUNCTIONS_H
#define LOCALFUNCTIONS_H

#include "fvCFD.H"
//#include "fvOptions.H"

void InitPsiXYZ(volScalarField& PsiZero, const fvMesh& mesh, scalar (*funIntP)(const double x,const double y, const double z) ) ;

void InitT(volScalarField& T, const volScalarField& PsiZero, const fvMesh& mesh, const dimensionedScalar epsH,
           scalar (*funInitT)(double psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z), double par
           );

void limitT(volScalarField& T,const fvMesh& mesh);

surfaceScalarField createPhiFieldEx(const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& PsiZero, const volScalarField& Psi );

surfaceScalarField createPhiFieldExC(const volScalarField& C, const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& Psi );

surfaceScalarField createPhiFieldImp(const volScalarField& C, const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& PsiZero, const volScalarField& Psi );

volScalarField createKField(string nazwa,const Time& runTime, const fvMesh& mesh);

void fileOpener(FILE** f, string outpath);     // wskaźnik wskazuje plik ale jest jak zmienna jezeli chce zmienic wartosc na ktory wskazuje to musze pobrac adres do wskaznika **

void AlphaToPsi(const volScalarField& T, volScalarField& Psi, const double& eps, const dimensionedScalar& epsH);

void AlphaToPsi2(const volScalarField& T, volScalarField& Psi, const double& eps, const dimensionedScalar& epsH, const fvMesh& mesh );

void AlphaToPsi3(const volScalarField& T, volScalarField& Psi, const double& eps, const dimensionedScalar& epsH);

void InitUXYZ( volVectorField& U, const fvMesh& mesh, vector (*funInitV) ( const double x, const double y, const double z) );

void PsiAlphaFace(const volScalarField& T,const volScalarField& Psi, const double& eps, const dimensionedScalar& epsH, const Time& runTime, const fvMesh& mesh,
                              surfaceScalarField& PsiF, surfaceScalarField& TF);
surfaceScalarField createSurfKField(string nazwa, const Time& runTime, const fvMesh& mesh);
// EXRK#
void updatemsnGradPsi(const fvMesh& mesh, const volScalarField & Psi, surfaceScalarField & msngradPsi);

void updatemGradPsi(const fvMesh& mesh, const volScalarField & Psi, volScalarField & mGradPsi);

void LimitGradPsi(const fvMesh& mesh, const volScalarField & Psi, volScalarField & gradPsi, double dx, const double gradPsiLimit, const volScalarField & PsiZero);

void LimitsnGradPsi(const fvMesh& mesh, const volScalarField & Psi, surfaceScalarField & msngradPsi, double dx, const double gradPsiLimit, const volScalarField & PsiZero, const Time& runTime);

#endif // LOCALFUNCTIONS_H
