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
              volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz, double gamma,
              const bool& mapFunLog, const double ilePkt, const double gradPsiLimit)
    {
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

    Info<< "Explicit RK3Face !!! "  << endl;
    scalar one(1);


    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CREATE FIELDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #include "createfieldsexrk3face.h"

    volScalarField k1 = createKField("k1", runTime, mesh);
    surfaceScalarField kF1 = createSurfKField("kF1",runTime, mesh);

    volScalarField k2  = createKField("k2", runTime, mesh);
    surfaceScalarField kF2 = createSurfKField("kF2",runTime, mesh);


    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OPEN FILES END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Psi == PsiZero;

    magGradT = mag(fvc::grad(T));
//    msnGradPsi == mag(fvc::snGrad(Psi));
//    mGradPsi = mag(fvc::grad(Psi));
//    surfgradT = linearInterpolate( T*(1.0-T)*(fvc::grad(Psi))/epsH ) & mesh.Sf();

//    GradTW = T*(1.0-T)*mag(fvc::grad(Psi))/epsH;
    gradTanalit = mag( ( one - Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH)) )*(1.0/(4.0*epsH)) );
//    LaplaceTStraignt = fvc::laplacian(T);
    TAnalit = 1.0/(1.0+Foam::exp(-PsiZero/epsH));
//    LaplaceTW1D = T*(1.0-T)*( fvc::laplacian(Psi) + (1.0/epsH)*( fvc::grad(Psi)&fvc::grad(Psi) )*(1.0-2.0*T) )/epsH;
//    LaplaceAnalit = - ( Foam::tanh(PsiZero/(2.0*epsH)) - Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH)) )/(4.0*epsH*epsH);


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
    for(int i = 0; i <= ilekrcz; ++i)
    {
        Told == T;

        limitT(T);
        AlphaToPsi(T, Psi, eps, epsH);
        PsiAlphaFace(T, Psi, eps, epsH, runTime, mesh, PsiF, TF);
        updateGradPsiF(mesh,Psi,mGradPsiF);
        LimitGradPsiF(mesh, PsiF, mGradPsiF, 1./ilePkt, gradPsiLimit);
        phiR2 = Cf*TF*( scalar(1.) - TF)*(mGradPsiF - scalar(1.))*( linearInterpolate(fvc::grad(Psi)) /( mGradPsiF )) & mesh.Sf();

        k1 ==  T + fvc::div(phiR2)*dtau;
        limitT(k1);
        AlphaToPsi3(k1, Psi, eps, epsH);
        PsiAlphaFace(k1, Psi, eps, epsH, runTime, mesh, PsiF, kF1);
        updateGradPsiF(mesh,Psi,mGradPsiF);
        LimitGradPsiF(mesh, PsiF, mGradPsiF, 1./ilePkt, gradPsiLimit);

        phiR2 = Cf*kF1*( scalar(1.) - kF1)*(mGradPsiF - scalar(1.))*( linearInterpolate(fvc::grad(Psi)) /( mGradPsiF )) & mesh.Sf();

        k2 == 3./4* T + 1./4* k1 + 1./4* fvc::div(phiR2)*dtau;
        limitT(k2);
        AlphaToPsi(k2, Psi, eps, epsH);
        PsiAlphaFace(k2, Psi, eps, epsH, runTime, mesh, PsiF, kF2);
        updateGradPsiF(mesh,Psi,mGradPsiF);
        LimitGradPsiF(mesh, PsiF, mGradPsiF, 1./ilePkt, gradPsiLimit);

        phiR2 = Cf*kF2*( scalar(1.) - kF2)*(fvc::snGrad(Psi) - scalar(1.))*( linearInterpolate(fvc::grad(Psi)) /( mGradPsiF )) & mesh.Sf();

        T ==  1./3* T + 2./3* k2 + 2./3*fvc::div(phiR2)*dtau;


        double norm1c = Foam::sum(Foam::mag(T-Told)).value() / T.size();
        Info << "Norma 1 = " << norm1c << endl;

        forAll(mesh.cellCentres(), cellI )
        {
            sGradT[cellI] = GradT[cellI].x();
//            Info << GradT[cellI].x() << GradT[cellI].y() << endl;
        }
    }
    double dx = 1.0/ilePkt;
    double norm1_grad_pkt = 0.0;
    double norm1_grad_pkt_1TW = 0.0;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SAVE VALUES t = tk !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    forAll(mesh.cellCentres(), cellI )
    {
        if ( (mesh.cellCentres()[cellI].y() > 0.0) && (mesh.cellCentres()[cellI].y() < dx) )
        {
//            f << mesh.cellCentres()[cellI].x() << " "<< T[cellI] << " " << TAnalit[cellI] << " "  << gradTanalit[cellI] << " " << magGradT[cellI] << " " << GradTW[cellI] << " " << LaplaceAnalit[cellI] << " " << LaplaceTW1D[cellI] << " " << LaplaceTStraignt[cellI] << " " << norm1_grad_pkt_1TW << " " << norm1_grad_pkt << std::endl;
//            f << mesh.cellCentres()[cellI].x() << " " << std::setprecision(32) << T[cellI] << " " << TAnalit[cellI] << " " << PsiZero[cellI] << std::endl;
        }
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
    double norm1 = Foam::sum(Foam::mag(TAnalit-T)).value() / T.size();
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
