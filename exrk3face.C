#include "exrk3face.h""
#include <iostream>
#include <fstream>
#include "basicexceptions.h"

dimensionedScalar SMALL_NUMBER2("small", dimless, SMALL);

void updateGradPsiF(const fvMesh& mesh, const volScalarField & Psi, surfaceScalarField & mGradPsiF)
{
    mGradPsiF == fvc::snGrad(Psi);

    forAll(mGradPsiF, FaceId)
    {
        if(mGradPsiF[FaceId] < SMALL_NUMBER2.value())
        {
            mGradPsiF[FaceId] = 1.;
        }
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvsPatchScalarField& gradPsiPatchF = mGradPsiF.boundaryField()[patchi];
        forAll(patch, faceId)                       // petla po centrach objetosci danego patcha
        {
            if(gradPsiPatchF[faceId] < SMALL_NUMBER2.value())
            {
                gradPsiPatchF[faceId] = 1.;
            }
        }
    }
}

void LimitGradPsiF(const fvMesh& mesh, const surfaceScalarField & PsiF, surfaceScalarField& mGradPsiF, const double dx, const double gradPsiLimit)
{
    forAll( mGradPsiF, FaceId )
    {
        if( (PsiF[FaceId] > gradPsiLimit*dx) || (PsiF[FaceId] < -gradPsiLimit*dx) )
        {
            mGradPsiF[FaceId] = 1.;
        }
    }
    forAll(mesh.boundary(), patchi)  // mesh.boundary() daje liste adresow do war. brzeg.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvsPatchScalarField& gradPsiPatchF = mGradPsiF.boundaryField()[patchi];
        const fvsPatchScalarField& PsiPatchF = PsiF.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if ( (PsiPatchF[faceId]  > gradPsiLimit*dx) || (PsiPatchF[faceId] < -gradPsiLimit*dx) )
            {
                gradPsiPatchF[faceId] = 1.;
            }
        }
    }
}

void exRK3Face(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau,const volScalarField& PsiZero, volScalarField& Psi,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz, double gamma, const bool& mapFunLog, const double ilePkt, const double gradPsiLimit)
    {
    std::fstream f;
    string nazwa, sciezka, metcalk;
    metcalk = "exRK3";
    std::ostringstream strs;
    strs << eps;
    std::string wielkosceps = strs.str();

    std::ostringstream strs2;
    strs2 << ilePkt;
    std::string ilePktStr = strs2.str();

    sciezka = "/home/sr/foam/sr-4.0/run/convectAlphaSTest2D_fe40";
    try
    {
    f.open("wynik.vtk", std::ios::out);
    if( f.good())
        Info << "plik otwarty" << endl;
    else
        throw FileException();
    }
    catch(BasicException& ex)
    { Info << ex.wyjatek << endl; }

    Info<< "Explicit RK3Face !!! "  << endl;
    scalar one(1);

    Psi == PsiZero;

    surfaceScalarField phiR2
    (
        IOobject
        (
            "phiR2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phiR2", (dimLength*dimLength*dimLength)/dimTime, scalar(1.))
    );

    volScalarField mGradPsi
    (
        IOobject
        (
            "gradPsi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(fvc::grad(Psi))
    );

    volScalarField gradTanalit
    (
        IOobject
        (
            "gradTanalit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag( ( one - Foam::tanh(PsiZero/(2*epsH))*Foam::tanh(PsiZero/(2*epsH)) )*(1/(4*epsH)) )
    );

    volScalarField magGradT
    (
        IOobject
        (
            "magGradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(fvc::grad(T))
    );

    volVectorField GradT
    (
        IOobject
        (
            "GradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T*(1-T)*(fvc::grad(Psi))/epsH
    );

    volScalarField sGradT
    (
        IOobject
        (
            "sGradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(fvc::grad(Psi))/epsH
    );

    volScalarField LaplaceT
    (
        IOobject
        (
            "LaplaceT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T*(1-T)*(fvc::laplacian(Psi))/epsH
    );

    volScalarField LaplaceTStraight
    (
        IOobject
        (
            "LaplaceTStraight",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::laplacian(T)
    );

    volScalarField TAnalit
    (
        IOobject
        (
            "TAnalit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*(1.+Foam::tanh(PsiZero/(2.*epsH)))
    );
    volScalarField TStartowe
    (
        IOobject
        (
            "TStartowe",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T
    );

    volScalarField LaplaceTW1D
    (
        IOobject
        (
            "LaplaceTW1D",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T*(1.-T)*( fvc::laplacian(Psi) + (1./epsH)*( fvc::grad(Psi)&fvc::grad(Psi) )*(1.-2.*T) )/epsH
     );

    volScalarField LaplaceAnalit
    (
        IOobject
        (
            "LaplaceAnalit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
            - ( Foam::tanh(PsiZero/(2.0*epsH)) - Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH)) )/(4*epsH*epsH)
     );

    surfaceScalarField surfgradT
    (
        IOobject
        (
            "surfgradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       linearInterpolate( T*(1-T)*(fvc::grad(Psi))/epsH ) & mesh.Sf()
    );

    volVectorField GradTW
    (
        IOobject
        (
            "GradTW",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T*(1.-T)*(fvc::grad(Psi))/epsH
    );

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NEW   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    surfaceScalarField PsiF
    (
        IOobject
        (
            "PsiF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("PsiF", dimLength, scalar(0))
    );

    surfaceScalarField TF
    (
        IOobject
        (
            "TF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("TF", dimless, scalar(0))
    );

    surfaceScalarField Cf
    (
        IOobject
        (
            "Cf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Cf", dimLength/dimTime, scalar(1.))
    );

    surfaceScalarField mGradPsiF
    (
        IOobject
        (
            "mGradPsiF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::snGrad(Psi)
    );

    for(int i = 0; i <= ilekrcz; ++i)
    {
        Told == T;

        AlphaToPsi(T, Psi, eps, epsH);
        PsiAlphaFace(T, Psi, eps, epsH, runTime, mesh, PsiF, TF);
        updateGradPsiF(mesh,Psi,mGradPsiF);
        LimitGradPsiF(mesh, PsiF, mGradPsiF, 1./ilePkt, gradPsiLimit);
        phiR2 = Cf*TF*( scalar(1.) - TF)*(mGradPsiF - scalar(1.))*( linearInterpolate(fvc::grad(Psi)) /( mGradPsiF )) & mesh.Sf();

        volScalarField k1 = createKField("k1", runTime, mesh);
        k1 ==  T + fvc::div(phiR2)*dtau;
        limitT(k1,mesh);
        AlphaToPsi(k1, Psi, eps, epsH);
        surfaceScalarField kF1 = createSurfKField("kF1",runTime, mesh);
        PsiAlphaFace(k1, Psi, eps, epsH, runTime, mesh, PsiF, kF1);
        updateGradPsiF(mesh,Psi,mGradPsiF);
        LimitGradPsiF(mesh, PsiF, mGradPsiF, 1./ilePkt, gradPsiLimit);

        phiR2 = Cf*kF1*( scalar(1.) - kF1)*(mGradPsiF - scalar(1.))*( linearInterpolate(fvc::grad(Psi)) /( mGradPsiF )) & mesh.Sf();

        volScalarField k2  = createKField("k2", runTime, mesh);
        k2 == 3./4* T + 1./4* k1 + 1./4* fvc::div(phiR2)*dtau;
        limitT(k2,mesh);
        AlphaToPsi(k2, Psi, eps, epsH);
        surfaceScalarField kF2 = createSurfKField("kF2",runTime, mesh);
        PsiAlphaFace(k2, Psi, eps, epsH, runTime, mesh, PsiF, kF2);
        updateGradPsiF(mesh,Psi,mGradPsiF);
        LimitGradPsiF(mesh, PsiF, mGradPsiF, 1./ilePkt, gradPsiLimit);

        phiR2 = Cf*kF2*( scalar(1.) - kF2)*(fvc::snGrad(Psi) - scalar(1.))*( linearInterpolate(fvc::grad(Psi)) /( mGradPsiF )) & mesh.Sf();

        T ==  1./3* T + 2./3* k2 + 2./3*fvc::div(phiR2)*dtau;
        limitT(T,mesh);

        double norm1c = Foam::sum(Foam::mag(T-Told)).value() / T.size();
        Info << "Norma 1 = " << norm1c << endl;

        forAll(mesh.cellCentres(), cellI )
        {
            sGradT[cellI] = GradT[cellI].x();
//            Info << GradT[cellI].x() << GradT[cellI].y() << endl;
        }
//        surfgradT = linearInterpolate( T*(1-T)*(fvc::grad(Psi))/epsH ) & mesh.Sf();
//        GradT = T*(1-T)*mag(fvc::grad(Psi))/epsH;

//        LaplaceT = fvc::div(surfgradT);
//        LaplaceTStraight = fvc::laplacian(T);

//        magGradT = mag(fvc::grad(T));
//        GradTW = T*(1.-T)*gradPsi/epsH;
//        LaplaceTStraight = fvc::laplacian(T);
//        LaplaceTW1D = T*(1.-T)*( fvc::laplacian(Psi) + (1./epsH)*( fvc::grad(Psi)&fvc::grad(Psi) )*(1.-2.*T) )/epsH;
//        gradTanalit =  mag( ( one - Foam::tanh(PsiZero/(2.*epsH))*Foam::tanh(PsiZero/(2.*epsH)) )*(1./(4.*epsH)) );

//        double norm1 = Foam::sum(Foam::mag(TAnalit-T)).value() / T.size();
//        double norm2 = Foam::pow(Foam::sum(Foam::pow(TAnalit-T, 2)).value(), 0.5 ) / T.size();
//        double norm3 = Foam::max(Foam::mag(TAnalit-T)).value();
//        Info << "Norma 1 analit = " << norm1 << endl;
//        double norm1gr = Foam::sum(Foam::mag(gradTanalit-GradTW)).value() / T.size();
//        double norm2gr = Foam::pow(Foam::sum(Foam::pow((gradTanalit-GradTW), 2)).value(), 0.5 ) / T.size();
//        double norm3gr = Foam::max(Foam::mag(gradTanalit-GradTW)).value();

//        double norm1lap = Foam::sum(Foam::mag(LaplaceAnalit-LaplaceTW1D)).value() / T.size();
//        double norm2lap = Foam::pow(Foam::sum(Foam::pow((LaplaceAnalit-LaplaceTW1D), 2)).value(), 0.5 ) / T.size();
//        double norm3lap = Foam::max(Foam::mag(LaplaceAnalit-LaplaceTW1D)).value();

//        f << i << " " << norm1 << " " << norm2 << " " << norm3 << " " << norm1gr << " " << norm2gr << " " << norm3gr << " " << norm1lap <<  " " << norm2lap << " " << norm3lap << " " << norm1c <<std::endl;

//        runTime.write();
//        ++runTime;
    }

    f.close();

//    try
//    {
//    f.open(sciezka + "gradAlpha_" + metcalk + "_eps_" + wielkosceps + "_" + ilePktStr + ".vtk" , std::ios::out);
//    if( f.good())
//        Info << "plik otwarty" << endl;
//    else
//        throw FileException();
//    }
//    catch(BasicException& ex)
//    { Info << ex.wyjatek << endl; }

//    double dx = 1./ilePkt;
//    double norm1_grad_pkt = 0.;
//    double norm1_grad_pkt_1TW = 0.;
//    forAll(mesh.cellCentres(), cellI )
//    {
//        if ( (mesh.cellCentres()[cellI].y() > 0) && (mesh.cellCentres()[cellI].y() < dx) )
//        {
//            norm1_grad_pkt_1TW = Foam::mag( gradTanalit[cellI] - GradTW[cellI] ) / ( Foam::mag( gradTanalit[cellI] ) + eps );
//            norm1_grad_pkt = Foam::mag( gradTanalit[cellI] - magGradT[cellI] ) / ( Foam::mag( gradTanalit[cellI] ) + eps );
//            f << mesh.cellCentres()[cellI].x() << " "<< T[cellI] << " " << TAnalit[cellI] << " "  << gradTanalit[cellI] << " " << magGradT[cellI] << " " <<GradTW[cellI] << " " << LaplaceAnalit[cellI] << " " << LaplaceTW1D[cellI] << " " << LaplaceTStraight[cellI] << " " << TStartowe[cellI] << " " << norm1_grad_pkt_1TW << " " << norm1_grad_pkt << std::endl;
//        }
//    }
//    f.close();

//   std::ofstream newFile(sciezka + "AnalizaZbierz.vtk", std::ios_base::app);

//       if(newFile.is_open())
//       {
//           double norm1 = Foam::sum(Foam::mag(TAnalit-T)).value() / T.size();
//           double norm2 = Foam::pow(Foam::sum(Foam::pow(TAnalit-T, 2)).value(), 0.5 ) / T.size();
//           double norm3 = Foam::max(Foam::mag(TAnalit-T)).value();

//           newFile << ilePktStr << " " << norm1 << " " << norm1 << " " << norm2 << " " << norm2 << " " << norm3 << " " << norm3 << std::endl;
//       }
//       else
//       {
//           Info << "wszyscy zginiemy!!! " << endl;
//       }
//       newFile.close();

}
