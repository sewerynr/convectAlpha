#include <stdlib.h>
#include "localFunctions.h"

dimensionedScalar SMALL_NUMBER3("small", dimless, SMALL);

void InitPsiXYZ(volScalarField& PsiZero, const fvMesh& mesh, scalar (*funIntP)(const double x,const double y, const double z) )
{
    forAll(mesh.cellCentres(),cellI) // loop over cell centers
    {
        PsiZero[cellI] = funIntP( mesh.cellCentres()[cellI].x(),mesh.cellCentres()[cellI].y(),mesh.cellCentres()[cellI].z() );
    }

    forAll(mesh.boundary(), patchi)  // loop over all boundary patches
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& psiPatch = PsiZero.boundaryField()[patchi];
//        psiPatch == psiPatch.patchInternalField();  // wymusza wartosc na krawedziachz wbrzpochodnych
//        psiPatch.evaluate();
        forAll(patch, faceId)       //  loop over c.v. centers of the given patch
        {
//            psiPatch[faceId] = PsiZero[ patch.faceCells()[faceId] ];   // wymusza wartosc na krawedziachz wbrzpochodnych
            psiPatch[faceId] = funIntP( patch.Cf()[faceId].x(),  patch.Cf()[faceId].y(), patch.Cf()[faceId].z());
        }
//        psiPatch.evaluate();
    }
}

void InitT(volScalarField& T, const volScalarField& PsiZero, const fvMesh& mesh, const dimensionedScalar epsH,
           scalar (*funInitT)(double psi, const dimensionedScalar& epsH, double par,  double x, double y,  double z) , double par
           )
{
    forAll(mesh.cellCentres(),cellI)
    {
        double p = PsiZero[cellI];
        T[cellI] = funInitT(p, epsH, par, mesh.cellCentres()[cellI].x(),mesh.cellCentres()[cellI].y(),mesh.cellCentres()[cellI].z() );
    }

    forAll(mesh.boundary(), patchi)  // mesh.boundary() returns list of addresses to b.c.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& TPatch = T.boundaryField()[patchi];
        const fvPatchScalarField& psiPatch = PsiZero.boundaryField()[patchi];

        forAll(patch, faceId) // loop over centers of faces belonging to given patchi
        {
            double pBon = psiPatch[faceId];
//???       Patch[faceId] = funInitT( pBon, epsH, par, mesh.cellCentres()[faceId].x(),mesh.cellCentres()[faceId].y(),mesh.cellCentres()[faceId].z() );
            TPatch[faceId] = funInitT( pBon, epsH, par, patch.Cf()[faceId].x(),patch.Cf()[faceId].y(), patch.Cf()[faceId].z() );
        }
    }
}


void InitUXYZ( volVectorField& U, const fvMesh& mesh, vector (*funInitV) (const double x, const double y, const double z) )
{
    forAll (mesh.cellCentres(), I)
    {
        U[I] = funInitV( mesh.cellCentres()[I].x(), mesh.cellCentres()[I].y(), mesh.cellCentres()[I].z());
    }
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        fvPatchVectorField& UPach = U.boundaryField()[patchI];

        forAll(patch, faceId)
        {
            UPach[faceId] = funInitV( mesh.cellCentres()[faceId].x(), mesh.cellCentres()[faceId].y(), mesh.cellCentres()[faceId].z());
        }
    }
    U.correctBoundaryConditions();
}

void limitT(volScalarField& T, const fvMesh& mesh)
{
    T = Foam::min(Foam::max(T, scalar(0.) ), scalar(1.) );
//    forAll(mesh.boundary(), patchi)
//    {
//        const fvPatch& patch = mesh.boundary()[patchi];
//        fvPatchScalarField& TPatch = T.boundaryField()[patchi];
//        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//        {
//            if(TPatch[faceId] > 0.9999999999999)
//            {
//                TPatch[faceId] = scalar(1.0);
//            }
//        }
//    }
}


surfaceScalarField createPhiFieldEx(const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& PsiZero, const volScalarField& Psi )
{
    dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);
    surfaceScalarField phiR  // cv surface!
        (
          IOobject
          (
              "phiR",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          linearInterpolate( T*(scalar(1.) - T)*(mag(fvc::grad(Psi))- scalar(1.) ) /( mag(fvc::grad(Psi)) + SMALL_NUMBER ) ) *  fvc::snGrad(Psi) * mesh.magSf()
        );
    return phiR;
}


surfaceScalarField createPhiFieldExC(const volScalarField& C, const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& Psi )
{
    surfaceScalarField phiR  // cv surface!
        (
          IOobject
          (
              "phiR",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),  mesh,
              dimensionedScalar("phiR", dimless/dimTime, scalar(1000000000.))// ?????????????????? l^3/s
        );
    return phiR;
}


surfaceScalarField createPhiFieldImp(const volScalarField& C, const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& PsiZero, const volScalarField& Psi )
{
    dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);
    surfaceScalarField phiR  // cv surface!
        (
          IOobject
          (
              "phiR",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          linearInterpolate( C*( scalar(1.) - T )*( mag(fvc::grad(Psi))- scalar(1.) ) * (fvc::grad(Psi) /( mag(fvc::grad(Psi)) + SMALL_NUMBER ))  )  & mesh.Sf()
        );
    return phiR;
}

volScalarField createKField(string nazwa, const Time& runTime, const fvMesh& mesh)
{
    volScalarField k
    (
        IOobject
        (
            nazwa,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(nazwa, dimless, scalar(0))
    );
    return k;
}


void fileOpener(FILE** f, string outpath)
{
    *f = fopen(outpath.c_str(), "r");
    if(NULL == f)
       {
         fprintf(stderr, "!!! fopen() _nr_1 failed.\n");
       }
}

void AlphaToPsi(const volScalarField& T, volScalarField& Psi, const double& eps, const dimensionedScalar& epsH)
{
//    maxPsi = 1.e+12;
//    minPsi = -1.e+12;

//    if(abs(Psi) < 26*epsH)
     Psi == epsH*Foam::log((T+eps)/(1-T+eps));

//     forAll(Psi.mesh().boundary(), patchi)  // mesh.boundary() daje liste adresow do war. brzeg.
//     {
//         const fvPatch& patch = Psi.mesh().boundary()[patchi];
//         fvPatchScalarField& psiPatch = Psi.boundaryField()[patchi];
//         const fvPatchScalarField& TPatch = T.boundaryField()[patchi];

//         forAll(patch, faceId) // petla po centrach objetosci danego patcha
//         {
//             psiPatch[faceId] = epsH.value()*Foam::log((TPatch[faceId]+eps)/(1-TPatch[faceId]+eps));
//         }
//     }

//     maxPsi = max(Psi,maxPsi);
//     minPsi = min(Psi,minPsi);
//    else if(Psi > 0)
//     Psi = maxPsiC;
//    else if(Psi < 0)
//     Psi = minPsiC;
//
//    maxPsiC=maxPsi;
//    minPsiC=minPsi;*/

}

void PsiAlphaFace(const volScalarField& T,const volScalarField& Psi, const double& eps, const dimensionedScalar& epsH, const Time& runTime, const fvMesh& mesh,
                  surfaceScalarField& PsiF, surfaceScalarField& TF)
{
    PsiF == linearInterpolate(Psi);
    TF == scalar(1.) / ( scalar(1.) + Foam::exp(-PsiF/epsH));
}

surfaceScalarField createSurfKField(string nazwa, const Time& runTime, const fvMesh& mesh)
{
    surfaceScalarField ks
    (
        IOobject
        (
            nazwa,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(nazwa, dimless, scalar(0))
    );
    return ks;
}

// EXRK3
void updatemsnGradPsi(const fvMesh& mesh, const volScalarField & Psi, surfaceScalarField & msngradPsi)
{
    msngradPsi == mag(fvc::snGrad(Psi));
    forAll(msngradPsi, faceId)
    {
        if(msngradPsi[faceId] < SMALL_NUMBER3.value())
        {
            msngradPsi[faceId] = scalar(1.0) ;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvsPatchScalarField& gradPsiPatch = msngradPsi.boundaryField()[patchi];
        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
        {
            if(gradPsiPatch[faceId] < SMALL_NUMBER3.value())
            {
                gradPsiPatch[faceId] = scalar(1.0);
            }
        }
    }
}

void updatemGradPsi(const fvMesh& mesh, const volScalarField & Psi, volScalarField & mGradPsi)
{
    mGradPsi == mag(fvc::grad(Psi));
    forAll(mGradPsi, CellId)
    {
        if(mGradPsi[CellId] < SMALL_NUMBER3.value())
        {
            mGradPsi[CellId] = scalar(1.0);
        }
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& gradPsiPatch = mGradPsi.boundaryField()[patchi];
        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
        {
            if(gradPsiPatch[faceId] < SMALL_NUMBER3.value())
            {
                gradPsiPatch[faceId] = scalar(1.0);
            }
        }
    }
}

void LimitGradPsi(const fvMesh& mesh, const volScalarField & Psi, volScalarField & gradPsi, double dx, const double gradPsiLimit, const volScalarField & PsiZero)
{
    forAll( mesh.cellCentres(), cellId)
    {
        if( mag(PsiZero[cellId]) > gradPsiLimit*dx )
        {
            gradPsi[cellId] = 1.0;
        }
    }

    forAll(mesh.boundary(), patchi)  // mesh.boundary() returns addresses to b.c.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& gradPsiPatch = gradPsi.boundaryField()[patchi];
        const fvPatchScalarField& PsiZeroPatch = PsiZero.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if ( mag(PsiZeroPatch[faceId]) > gradPsiLimit*dx )
            {
                gradPsiPatch[faceId] = 1.0;
            }
        }
    }
}
void LimitsnGradPsi(const fvMesh& mesh, const volScalarField & Psi, surfaceScalarField & msngradPsi, double dx, const double gradPsiLimit, const volScalarField & PsiZero, const Time& runTime)
{
    surfaceScalarField sPsiZero
    (
        IOobject
        (
            "sPsiZero",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(PsiZero)
    );

    forAll(sPsiZero, faceId)
    {
        if( mag(sPsiZero[faceId]) > gradPsiLimit*dx )
        {
            msngradPsi[faceId] = scalar(1.0) ;
        }

    }

    forAll(mesh.boundary(), patchi)  // mesh.boundary() returns addresses to b.c.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvsPatchScalarField& msngradPsiPatch = msngradPsi.boundaryField()[patchi];
      //  const fvPatchScalarField& PsiZeroPatch = PsiZero.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if ( mag(patch.Cf()[faceId].x()) > gradPsiLimit*dx )
            {
                msngradPsiPatch[faceId] = 1.0;
            }
        }
    }
}