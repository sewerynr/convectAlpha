# include "implicit.h"
#include <iostream>
#include <fstream>
#include "basicexceptions.h"

void implicit(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau, const volScalarField& PsiZero, volScalarField& Psi,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilkrcz, double gamma, int ilePkt, const double gradPsiLimit)
{
    Info<< "Implicit !!! "  << endl;

    scalar one(1);
    dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);

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
     //   mesh,
      //  dimensionedScalar("TAnalit", dimless, scalar(-100000000.0))
        0.5*(1.+Foam::tanh(PsiZero/(2.*epsH)))
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
        mesh,
        dimensionedScalar("gradPsi", dimless, scalar(-100000000.0))
    );

//    volVectorField GradPsi
//    (
//        IOobject
//        (
//            "gradPsi",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh,
//        dimensionedVector("gradPsi", dimless, vector(-100000000.0, 0., 0.))
//    );

    string nazwa, sciezka, metcalk;
    metcalk = "exRK3";
    std::ostringstream strs;
    strs << eps;
    std::string wielkosceps = strs.str();

    std::ostringstream strs2;
    strs2 << ilePkt;
    std::string ilePktStr = strs2.str();

    sciezka = "/home/sr/foam/sr-4.0/run/convectAlphaSTest2D_fe40";

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OPEN FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::fstream f, fB, f0, fB0, fC, fAnZbie;
    try
    {
    fC.open("./convergence.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    f.open("./gradAlpha.vtk" , std::fstream::in | std::fstream::out | std::fstream::trunc);
    fB.open("./gradAlphaB.vtk" , std::fstream::in | std::fstream::out | std::fstream::trunc);
    f0.open("./gradAlpha0.vtk" , std::fstream::in | std::fstream::out | std::fstream::trunc);
    fB0.open("./gradAlphaB0.vtk" , std::fstream::in | std::fstream::out | std::fstream::trunc);
    fAnZbie.open("./AnalizaZbierz.vtk", std::fstream::in | std::fstream::out | std::fstream::app);

    if( f.good() && fB.good() && f0.good() && fB0.good() && fC.good() && fAnZbie.good())
        Info << "plik otwarty" << endl;
    else
        throw FileException();
    }
    catch(BasicException& ex)
    { Info << ex.wyjatek << endl; }

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SAVE VALUES t = t0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        forAll(mesh.cellCentres(), cellI )
        {
            f0 << mesh.cellCentres()[cellI].x()
               << " " << mesh.cellCentres()[cellI].y()
               << " " << std::setprecision(18) << T[cellI]
               << " " << std::setprecision(18) << TAnalit[cellI]
               << " " << std::setprecision(18) << PsiZero[cellI]
               << " " << std::setprecision(18) << Psi[cellI]
               << " " << std::setprecision(18) << mGradPsi[cellI] << std::endl;
        }
        forAll(mesh.boundary(), patchi )
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            if( patch.type() == "empty" )
            { continue; }
            fvPatchScalarField& TPatch = T.boundaryField()[patchi];
            fvPatchScalarField& TAPatch = TAnalit.boundaryField()[patchi];
            const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
            fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];
            fvPatchScalarField& PmGradPsi = mGradPsi.boundaryField()[patchi];

            forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
            {
                fB0 << patch.Cf()[faceId].x()
                    << " " << patch.Cf()[faceId].y()
                    << " " << std::setprecision(18) << TPatch[faceId]
                    << " " << std::setprecision(18) << TAPatch[faceId]
                    << " " << std::setprecision(18) << PZPatch[faceId]
                    << " " << std::setprecision(18) << PPatch[faceId]
                    << " " << std::setprecision(18) << PmGradPsi[faceId]  << std::endl;
            }
        }
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! reinitialization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Psi = epsH*Foam::log((T+eps)/(1-T+eps));
    surfaceScalarField phiR = createPhiFieldImp(C, runTime, mesh, T, PsiZero, Psi);

    fvScalarMatrix InitTEq
    (
        fvm::ddt(T)
        -fvm::div(phiR*dtau, T)
    );

    scalar res = Foam::max(Foam::mag(InitTEq.residual()));
    Info << "Start residual: " << res <<endl;
    double norm1 = 0.;
    double norm1c = 0.;

    for (int i = 1; i <= ilkrcz; ++i)
    {

        Told = T;
        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
            -fvm::div(phiR*dtau, T)
        );

        TEqn.solve();
        limitT(T);
        ++runTime;
     //   runTime.write();

        AlphaToPsi(T, Psi, eps, epsH);
//        AlphaToPsi2(T, Psi, eps, epsH, mesh);
//        AlphaToPsi3(T, Psi, eps, epsH);
        updatemGradPsi(Psi, mGradPsi);
        LimitGradPsi(mesh, Psi, mGradPsi, 1./ilePkt, gradPsiLimit, PsiZero);
        phiR = linearInterpolate( C*( scalar(1.) - T )*( mGradPsi- scalar(1.) ) * (fvc::grad(Psi) /( mGradPsi + SMALL_NUMBER ))  )  & mesh.Sf();

        norm1c = Foam::sum(Foam::mag(T-Told)).value() / T.size();
        norm1 = Foam::sum(Foam::mag(TAnalit-T)).value() / T.size();
        Info << "Norma conv = " << norm1c <<  "  Norma 1 analit = " << norm1 << "  krok=  " << i << endl;
//        double norm2 = Foam::pow(Foam::sum(Foam::pow(TAnalit-T, 2)).value(), 0.5 ) / T.size();
//        double norm3 = Foam::max(Foam::mag(TAnalit-T)).value();
    }

    double dx = 1.0/ilePkt;
    double norm1_grad_pkt = 0.0;
    double norm1_grad_pkt_1TW = 0.0;
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SAVE VALUES t = tk !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    forAll(mesh.cellCentres(), cellI)
    {
//        if ( (mesh.cellCentres()[cellI].y() > 0.0) && (mesh.cellCentres()[cellI].y() < dx) )
//        {
//            f << mesh.cellCentres()[cellI].x() << " "<< T[cellI] << " " << TAnalit[cellI] << " "  << gradTanalit[cellI] << " " << magGradT[cellI] << " " << GradTW[cellI] << " " << LaplaceAnalit[cellI] << " " << LaplaceTW1D[cellI] << " " << LaplaceTStraignt[cellI] << " " << norm1_grad_pkt_1TW << " " << norm1_grad_pkt << std::endl;
//            f << mesh.cellCentres()[cellI].x() << " " << std::setprecision(32) << T[cellI] << " " << TAnalit[cellI] << " " << PsiZero[cellI] << std::endl;
//        }
        f << mesh.cellCentres()[cellI].x()
          << " " << mesh.cellCentres()[cellI].y()
          << " " << std::setprecision(18) << T[cellI]
          << " " << std::setprecision(18) << TAnalit[cellI]
          << " " << std::setprecision(18) << PsiZero[cellI]
          << " " << std::setprecision(18) << Psi[cellI]
          << " " << std::setprecision(18) << mGradPsi[cellI]  << std::endl;
    }
    forAll(mesh.boundary(), patchi )
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        if( patch.type() == "empty" )
        {
            continue;
        }
        fvPatchScalarField& TPatch = T.boundaryField()[patchi];
        fvPatchScalarField& TAPatch = TAnalit.boundaryField()[patchi];
        const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
        fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];
        fvPatchScalarField& PmGradPsi = mGradPsi.boundaryField()[patchi];

        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
        {
            fB << patch.Cf()[faceId].x()
               << " " << patch.Cf()[faceId].y()
               << " " << std::setprecision(18) << TPatch[faceId]
               << " " << std::setprecision(18) << TAPatch[faceId]
               << " " << std::setprecision(18) << PZPatch[faceId]
               << " " << std::setprecision(18) << PPatch[faceId]
               << " " << std::setprecision(18) << PmGradPsi[faceId]  << std::endl;
        }
    }
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AnalizaZbierz save !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    norm1 = Foam::sum(Foam::mag(TAnalit-T)).value() / T.size();
    double norm2 = Foam::pow(Foam::sum(Foam::pow(TAnalit-T, 2)).value(), 0.5 ) / T.size();
    double norm3 = Foam::max(Foam::mag(TAnalit-T)).value();
    fAnZbie << ilePktStr << " " <<  std::setprecision(18) << norm1 << " " << norm1 << " " << norm2 << " " << norm2 << " " << norm3 << " " << norm3 << std::endl;

    fAnZbie.close();
    fC.close();
    f.close();
    fB.close();
    f0.close();
    fB0.close();
}
