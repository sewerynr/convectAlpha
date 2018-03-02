#include "exrk3.h"
#include <iostream>
#include <fstream>
#include "basicexceptions.h"

dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);

void exRK3(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau, const volScalarField& PsiZero, volScalarField& Psi,
                 volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz,
                 double gamma, const bool& mapFunLog, const double ilePkt, const double gradPsiLimit)
    {
    scalar one(1.0);
    Info<< "Explicit RK3 !!! "  << endl;
    Info<< "N_Tau = " << ilekrcz  << endl;
    std::ostringstream strs;
    strs << eps;
    std::string wielkosceps = strs.str();

    std::ostringstream strs2;
    strs2 << ilePkt;
    std::string ilePktStr = strs2.str();

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
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OPEN FILES END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Psi == PsiZero;

    #include "createfieldsexrk3.H"

    magGradT = mag(fvc::grad(T));
 //   msnGradPsi == mag(fvc::snGrad(Psi));
  //  mGradPsi = mag(fvc::grad(Psi));
//    surfgradT = linearInterpolate( T*(1.0-T)*(fvc::grad(Psi))/epsH ) & mesh.Sf();

//    GradTW = T*(1.0-T)*mag(fvc::grad(Psi))/epsH;
    gradTanalit = mag( ( one - Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH)) )*(1.0/(4.0*epsH)) );
//    LaplaceTStraignt = fvc::laplacian(T);
    TAnalit = 1.0/(1.0+Foam::exp(-PsiZero/epsH));
//    LaplaceTW1D = T*(1.0-T)*( fvc::laplacian(Psi) + (1.0/epsH)*( fvc::grad(Psi)&fvc::grad(Psi) )*(1.0-2.0*T) )/epsH;
//    LaplaceAnalit = - ( Foam::tanh(PsiZero/(2.0*epsH)) - Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH))*Foam::tanh(PsiZero/(2.0*epsH)) )/(4.0*epsH*epsH);

    updatemsnGradPsi(Psi, msnGradPsi);
    updatemGradPsi(Psi, mGradPsi);

    LimitGradPsi(mesh, Psi, mGradPsi, 1.0/ilePkt, gradPsiLimit, PsiZero);
    LimitsnGradPsi(mesh, Psi, msnGradPsi, 1.0/ilePkt, gradPsiLimit, PsiZero, runTime);

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SAVE VALUES t = t0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    forAll(mesh.cellCentres(), cellI )
    {
        f0 << mesh.cellCentres()[cellI].x()
           << " " << mesh.cellCentres()[cellI].y()
           << " " << std::setprecision(32) << T[cellI]
           << " " << std::setprecision(32) << TAnalit[cellI]
           << " " << std::setprecision(32) << PsiZero[cellI]
           << " " << std::setprecision(32) << Psi[cellI]
           << " " << std::setprecision(32) << mGradPsi[cellI] << std::endl;
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
                << " " << std::setprecision(32) << TPatch[faceId]
                << " " << std::setprecision(32) << TAPatch[faceId]
                << " " << std::setprecision(32) << PZPatch[faceId]
                << " " << std::setprecision(32) << PPatch[faceId]
                << " " << std::setprecision(32) << PmGradPsi[faceId]  << std::endl;
        }
    }
    volScalarField k1 = createKField("k1", runTime, mesh);
    volScalarField k2 = createKField("k2", runTime, mesh);
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! reinitialization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for(int i = 0; i < ilekrcz; ++i)
    {
        Info<<"reinitialization iter "<< i << endl;

        Told == T; // == means copy with b.c.

        phiR = linearInterpolate( C*T*( scalar(1.0) - T ) * (fvc::grad(Psi) / ( mGradPsi )) ) & mesh.Sf();
        phiR == phiR * ( msnGradPsi - scalar(1.0) );

        k1 == T + fvc::div(phiR)*dtau;

        limitT(k1);
//        AlphaToPsi(k1, Psi, eps, epsH);
//        AlphaToPsi2(k1, Psi, eps, epsH, mesh);
        AlphaToPsi3(k1, Psi, eps, epsH);

        updatemGradPsi(Psi, mGradPsi);
        updatemsnGradPsi(Psi, msnGradPsi);
        LimitGradPsi(mesh, Psi, mGradPsi, 1./ilePkt, gradPsiLimit, PsiZero);
        LimitsnGradPsi(mesh, Psi, msnGradPsi, 1.0/ilePkt, gradPsiLimit, PsiZero, runTime);

        phiR = linearInterpolate( C*k1*( scalar(1.0) - k1 ) * (fvc::grad(Psi) / ( mGradPsi )) ) & mesh.Sf();
        phiR == phiR * ( msnGradPsi - scalar(1.0) );

        k2 == 3.0/4.0* T + 1.0/4.0* k1 + 1.0/4.0* fvc::div(phiR)*dtau;

        limitT(k2);
//        AlphaToPsi(k1, Psi, eps, epsH);
//        AlphaToPsi2(k1, Psi, eps, epsH, mesh);
        AlphaToPsi3(k2, Psi, eps, epsH);

        updatemsnGradPsi(Psi, msnGradPsi);
        updatemGradPsi(Psi, mGradPsi);
        LimitGradPsi(mesh, Psi, mGradPsi, 1./ilePkt, gradPsiLimit, PsiZero);
        LimitsnGradPsi(mesh, Psi, msnGradPsi, 1.0/ilePkt, gradPsiLimit, PsiZero, runTime);


        phiR = linearInterpolate( C*k2*( scalar(1.0) - k2 ) * (fvc::grad(Psi) / ( mGradPsi )) ) & mesh.Sf();
        phiR == phiR * ( msnGradPsi - scalar(1.0) );

        T ==  1.0/3.0* Told + 2.0/3.0* k2 + 2.0/3.0*fvc::div(phiR)*dtau;

        limitT(T);
//        AlphaToPsi(k1, Psi, eps, epsH);
//        AlphaToPsi2(k1, Psi, eps, epsH, mesh);
        AlphaToPsi3(T, Psi, eps, epsH);

        updatemsnGradPsi(Psi, msnGradPsi);
        updatemGradPsi(Psi, mGradPsi);
        LimitGradPsi(mesh, Psi, mGradPsi, 1./ilePkt, gradPsiLimit, PsiZero);
        LimitsnGradPsi(mesh, Psi, msnGradPsi, 1.0/ilePkt, gradPsiLimit, PsiZero, runTime);

        double norm1c = Foam::sum(Foam::mag(T-Told)).value() / T.size();
        Info << "Norma conv = " << norm1c << endl;
//        surfgradT = linearInterpolate( T*(1.0-T)*(fvc::grad(Psi))/epsH ) & mesh.Sf();

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NORMY ANALITYCZNE w KAZDYM KROKU dTau !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//        magGradT == mag(fvc::grad(T));
//        GradTW = T*(1.0-T)*mGradPsi/epsH;
//        LaplaceTStraignt = fvc::laplacian(T);
//        LaplaceTW1D = T*(1.0-T)*( fvc::laplacian(Psi) + (1.0/epsH)*( fvc::grad(Psi)&fvc::grad(Psi) )*(1.0-2.0*T) )/epsH;

        double norm1 = Foam::sum(Foam::mag(TAnalit-T)).value() / T.size();
//        double norm2 = Foam::pow(Foam::sum(Foam::pow(TAnalit-T, 2)).value(), 0.5 ) / T.size();
//        double norm3 = Foam::max(Foam::mag(TAnalit-T)).value();
        Info << "Norma 1 analit = " << norm1 << endl;
//        double norm1gr = Foam::sum(Foam::mag(gradTanalit-GradTW)).value() / T.size();
//        double norm2gr = Foam::pow(Foam::sum(Foam::pow((gradTanalit-GradTW), 2)).value(), 0.5 ) / T.size();
//        double norm3gr = Foam::max(Foam::mag(gradTanalit-GradTW)).value();

//        double norm1lap = Foam::sum(Foam::mag(LaplaceAnalit-LaplaceTW1D)).value() / T.size();
//        double norm2lap = Foam::pow(Foam::sum(Foam::pow((LaplaceAnalit-LaplaceTW1D), 2)).value(), 0.5 ) / T.size();
//        double norm3lap = Foam::max(Foam::mag(LaplaceAnalit-LaplaceTW1D)).value();

//        fC << i << " " << norm1 << " " << norm2 << " " << norm3 << " " << norm1gr << " " << norm2gr << " " << norm3gr << " " << norm1lap <<  " " << norm2lap << " " << norm3lap << " " << norm1c <<std::endl;
        runTime.write();
//        ++runTime;
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
