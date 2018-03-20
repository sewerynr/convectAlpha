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

void limitT(volScalarField& T)
{
    T = Foam::min(Foam::max(T, scalar(0.) ), scalar(1.) );
//    forAll(T.mesh().boundary(), patchi)
//    {
//        const fvPatch& patch = T.mesh().boundary()[patchi];
//        fvPatchScalarField& TPatch = T.boundaryField()[patchi];
//        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//        {
//            if(TPatch[faceId] > 1.0)
//            {
//                TPatch[faceId] = scalar(1.0);
//            }
//            else if(TPatch[faceId] < 0.0)
//            {
//                TPatch[faceId] = scalar(0.0);
//            }
//        }
//    }

//    forAll(T, celli)
//    {
//        if(T[celli] > 1.0)
//        {
//            T[celli] = scalar(1.0);
//        }
//        else if( T[celli] < 0.0)
//        {
//             T[celli] = scalar(0.0);
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
          linearInterpolate( T*(scalar(1.0) - T)*(mag(fvc::grad(Psi))- scalar(1.0) ) /( mag(fvc::grad(Psi)) + SMALL_NUMBER ) ) *  fvc::snGrad(Psi) * mesh.magSf()
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

void AlphaToPsi2(const volScalarField& T, volScalarField& Psi, const double& eps, const dimensionedScalar& epsH, const fvMesh& mesh )
{
    double small = 10.;
    double mult = 1.;
     forAll(T, CellId)
     {
         if(T[CellId] > SMALL_NUMBER3.value()*small && T[CellId] < (1.0 - SMALL_NUMBER3.value()*small  ) )
         {
            Psi[CellId] = epsH.value()*Foam::log((T[CellId])/(1.0 - T[CellId]));
         }
         else
         {
            Psi[CellId] = epsH.value()*Foam::log((T[CellId]+SMALL_NUMBER3.value()*mult)/(1.0 - T[CellId] + SMALL_NUMBER3.value()*mult));
         }
     }
     forAll(mesh.boundary(), patchi)
     {
         const fvPatch& patch = mesh.boundary()[patchi];
         fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];
         const fvPatchScalarField& TPatch = T.boundaryField()[patchi];
         forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
         {
             if(TPatch[faceId] > SMALL_NUMBER3.value()*small  && TPatch[faceId] < (1.0 - SMALL_NUMBER3.value()*small ) )
             {
                 PsiPatch[faceId] = epsH.value()*Foam::log((TPatch[faceId])/(1.0 - TPatch[faceId]));
             }
             else
             {
                 PsiPatch[faceId] = epsH.value()*Foam::log((TPatch[faceId]+SMALL_NUMBER3.value()*mult)/(1.0 - TPatch[faceId] + SMALL_NUMBER3.value()*mult));
             }
         }
     }
}

void AlphaToPsi3(const volScalarField& T, volScalarField& Psi, const double& eps, const dimensionedScalar& epsH)
{
    double small = 1.;
    double mult = 1.;
     forAll(T, CellId)
     {
         if(T[CellId] > SMALL_NUMBER3.value()*small && T[CellId] < (1.0 - SMALL_NUMBER3.value()*small  ) )
         {
            Psi[CellId] = epsH.value()*( Foam::log(T[CellId]) - Foam::log(1.0 - T[CellId])) ;
         }
         else
         {
            Psi[CellId] = epsH.value()*( Foam::log(T[CellId] +SMALL_NUMBER3.value()*mult) - Foam::log(1.0 - T[CellId] +SMALL_NUMBER3.value()*mult)) ;
         }
     }
     forAll(T.mesh().boundary(), patchi)
     {
         const fvPatch& patch = T.mesh().boundary()[patchi];
         fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];
         const fvPatchScalarField& TPatch = T.boundaryField()[patchi];
         forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
         {
             if(TPatch[faceId] > SMALL_NUMBER3.value()*small  && TPatch[faceId] < (1.0 - SMALL_NUMBER3.value()*small ) )
//                 if(TPatch[faceId] > 12  && TPatch[faceId] < (1.0 - 12 ) )
             {
                 PsiPatch[faceId] = epsH.value()*( Foam::log(TPatch[faceId]) - Foam::log(1.0 - TPatch[faceId]));
             }
             else
             {
                 PsiPatch[faceId] = epsH.value()*( Foam::log(TPatch[faceId]+SMALL_NUMBER3.value()*mult) - Foam::log(1.0 - TPatch[faceId]+SMALL_NUMBER3.value()*mult));
             }
         }
     }
}

void AlphaToPsi(const volScalarField& T, volScalarField& Psi, const double& eps, const dimensionedScalar& epsH)
{
//    if(abs(Psi) < 26*epsH)

    Psi == epsH*Foam::log((T + SMALL_NUMBER3/10000.)/(1.0 - T + SMALL_NUMBER3/10000.));

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
        dimensionedScalar(nazwa, dimless, scalar(0.))
    );
    return ks;
}

// EXRK3
void updatemsnGradPsi(const volScalarField & Psi, surfaceScalarField & msngradPsi)
{
    msngradPsi == mag(fvc::snGrad(Psi));
    forAll(msngradPsi, faceId)
    {
        if(msngradPsi[faceId] < SMALL_NUMBER3.value())
        {
            msngradPsi[faceId] = scalar(1.0) ;
        }
    }

//    forAll(Psi.mesh().boundary(), patchi)
//    {
//        const fvPatch& patch = Psi.mesh().boundary()[patchi];
//        fvsPatchScalarField& gradPsiPatch = msngradPsi.boundaryField()[patchi];
//        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//        {
//            if(gradPsiPatch[faceId] < SMALL_NUMBER3.value())
//            {
//                gradPsiPatch[faceId] = scalar(1.0);
//            }
//        }
//    }
// to pomaga !!!!!!!!!!!!! ?????????????????????? chyba
    forAll(msngradPsi, faceId)
    {
        if(msngradPsi[faceId] > scalar(1.0))
        {
            msngradPsi[faceId] = scalar(1.0) ;
        }
    }

//    forAll(Psi.mesh().boundary(), patchi)
//    {
//        const fvPatch& patch = Psi.mesh().boundary()[patchi];
//        fvsPatchScalarField& gradPsiPatch = msngradPsi.boundaryField()[patchi];
//        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//        {
//            if(gradPsiPatch[faceId] > scalar(1.0))
//            {
//                gradPsiPatch[faceId] = scalar(1.0);
//            }
//        }
//    }
}

void updatemGradPsi(const volScalarField & Psi, volScalarField & mGradPsi)
{
    mGradPsi == mag(fvc::grad(Psi));
    forAll(mGradPsi, CellId)
    {
        if(mGradPsi[CellId] < SMALL_NUMBER3.value())
        {
            mGradPsi[CellId] = scalar(1.0);
        }
    }
    forAll(Psi.mesh().boundary(), patchi)
    {
        const fvPatch& patch = Psi.mesh().boundary()[patchi];
        fvPatchScalarField& gradPsiPatch = mGradPsi.boundaryField()[patchi];
        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
        {
            if(gradPsiPatch[faceId] < SMALL_NUMBER3.value())
            {
                gradPsiPatch[faceId] = scalar(1.0);
            }
        }
    }

//    forAll(mGradPsi, CellId)
//    {
//        if(mGradPsi[CellId] > scalar(1.0))
//        {
//            mGradPsi[CellId] = scalar(1.0);
//        }
//    }
//    forAll(Psi.mesh().boundary(), patchi)
//    {
//        const fvPatch& patch = Psi.mesh().boundary()[patchi];
//        fvPatchScalarField& gradPsiPatch = mGradPsi.boundaryField()[patchi];
//        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//        {
//            if(gradPsiPatch[faceId] > scalar(1.0))
//            {
//                gradPsiPatch[faceId] = scalar(1.0);
//            }
//        }
//    }
}

void LimitGradPsi(const fvMesh& mesh, const volScalarField & Psi, volScalarField & gradPsi, double dx, const double gradPsiLimit, const volScalarField & PsiZero)
{
    forAll( mesh.cellCentres(), cellId)
    {
        if( mag(Psi[cellId]) > gradPsiLimit*dx )
        {
            gradPsi[cellId] = scalar(1.0);
        }
    }

    forAll(mesh.boundary(), patchi)  // mesh.boundary() returns addresses to b.c.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& gradPsiPatch = gradPsi.boundaryField()[patchi];
        const fvPatchScalarField& PsiZeroPatch = PsiZero.boundaryField()[patchi];
        const fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if ( mag(PsiPatch[faceId]) > gradPsiLimit*dx )
            {
                gradPsiPatch[faceId] = scalar(1.0);
            }
        }
    }
}
void LimitsnGradPsi(const fvMesh& mesh, const volScalarField & Psi, surfaceScalarField & msngradPsi, double dx, const double gradPsiLimit, const volScalarField & PsiZero, const Time& runTime)
{
//    surfaceScalarField sPsiZero
//    (
//        IOobject
//        (
//            "sPsiZero",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        linearInterpolate(PsiZero)
//    );

    surfaceScalarField sPsi
    (
        IOobject
        (
            "sPsi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(Psi)
    );

    forAll(sPsi, faceId)
    {
        if( mag(sPsi[faceId]) > gradPsiLimit*dx )
        {
            msngradPsi[faceId] = scalar(1.0) ;
        }
    }

    forAll(mesh.boundary(), patchi)  // mesh.boundary() returns addresses to b.c.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvsPatchScalarField& msngradPsiPatch = msngradPsi.boundaryField()[patchi];
        const fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if ( mag(PsiPatch[faceId]) > gradPsiLimit*dx )
            {
                msngradPsiPatch[faceId] = scalar(1.0) ;
            }
        }
    }
}
