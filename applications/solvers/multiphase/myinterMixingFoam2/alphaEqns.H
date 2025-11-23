{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // 界面压缩通量：沿用 interface.nHatf()（已在类中有 deltaN 护栏）
    surfaceScalarField phir
    (
        IOobject("phir", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        phic*interface.nHatf()
    );

    // 相变源（已在 mixture->vDotAlphal() 内做了 Δt 限幅；单位 [1/s]）
    Pair<tmp<volScalarField>> vDotAlphal = mixture->vDotAlphal();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();   // vaporisation / pyrolysis
    // condensation 端为零：const volScalarField& vDotcAlphal = vDotAlphal[0]();

    // MULES 预解/修正阶段不引入源项
    volScalarField SpZero
    (
        IOobject("SpZero", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar("zero", dimless/dimTime, 0.0)
    );

    // 预解：仅对 α1 做 ddt+conv，不含源项
    tmp<surfaceScalarField> talphaPhi;
    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (mesh, phi, upwind<scalar>(mesh, phi)).fvmDiv(phi, alpha1)
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << nl;

        talphaPhi = alpha1Eqn.flux();
    }

    volScalarField alpha10("alpha10", alpha1);

    // α1 的 MULES 修正（仍不带源）
    for (int aCorr = 0; aCorr < nAlphaCorr; ++aCorr)
    {
        tmp<surfaceScalarField> talphaPhiCorr
        (
            fvc::flux(phi, alpha1, alphaScheme)
          + fvc::flux(-fvc::flux(-phir, alpha2, alpharScheme), alpha1, alpharScheme)
          + fvc::flux(-fvc::flux(-phir, alpha3, alpharScheme), alpha1, alpharScheme)
        );

        if (MULESCorr)
        {
            talphaPhiCorr.ref() -= talphaPhi();

            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),     // "psi" = 1
                alpha1,
                talphaPhi.ref(),
                talphaPhiCorr.ref(),
                SpZero,                                 // Sp = 0
                (divU*(alpha10 - alpha100))(),          // Su
                1,                                      // nLimiterIter
                0
            );

            if (aCorr == 0)
            {
                talphaPhi.ref() += talphaPhiCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;      // 稳定化
                talphaPhi.ref() += 0.5*talphaPhiCorr();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                talphaPhiCorr.ref(),
                SpZero,
                (divU*alpha1)(),
                1,
                0
            );

            talphaPhi = talphaPhiCorr;
        }
    }

    // 最终修正：同时推进 α1（带源）与 α2（可含 α2-α3 扩散）
    for (int gCorr = 0; gCorr < nAlphaCorr; ++gCorr)
    {
        // ---------- α1：构造非有界通量并生成独立 limiter ----------
        surfaceScalarField alphaPhi1
        (
            fvc::flux(phi, alpha1, alphaScheme)
          + fvc::flux(-fvc::flux(-phir, alpha2, alpharScheme), alpha1, alpharScheme)
          + fvc::flux(-fvc::flux(-phir, alpha3, alpharScheme), alpha1, alpharScheme)
        );
        surfaceScalarField alphaPhi1BD( upwind<scalar>(mesh, phi).flux(alpha1) );
        alphaPhi1 -= alphaPhi1BD;

        scalarField lam1F(mesh.nFaces(), 1.0);
        slicedSurfaceScalarField lam1
        (
            IOobject("lambda1", mesh.time().timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh, dimless, lam1F, false
        );

        if (LTS)
        {
            const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);
            MULES::limiter(lam1F, rDeltaT, geometricOneField(),
                           alpha1, alphaPhi1BD, alphaPhi1, zeroField(), zeroField(), 1, 0);
        }
        else
        {
            MULES::limiter(lam1F, 1.0/runTime.deltaT().value(), geometricOneField(),
                           alpha1, alphaPhi1BD, alphaPhi1, zeroField(), zeroField(), 1, 0);
        }

        // ---------- α2：独立 limiter（不再与 α1 复用/覆盖） ----------
        surfaceScalarField alphaPhi2
        (
            fvc::flux(phi, alpha2, alphaScheme)
          + fvc::flux(-fvc::flux(phir, alpha1, alpharScheme), alpha2, alpharScheme)
        );
        surfaceScalarField alphaPhi2BD( upwind<scalar>(mesh, phi).flux(alpha2) );
        alphaPhi2 -= alphaPhi2BD;

        scalarField lam2F(mesh.nFaces(), 1.0);
        slicedSurfaceScalarField lam2
        (
            IOobject("lambda2", mesh.time().timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh, dimless, lam2F, false
        );

        if (LTS)
        {
            const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);
            MULES::limiter(lam2F, rDeltaT, geometricOneField(),
                           alpha2, alphaPhi2BD, alphaPhi2, zeroField(), zeroField(), 1, 0);
        }
        else
        {
            MULES::limiter(lam2F, 1.0/runTime.deltaT().value(), geometricOneField(),
                           alpha2, alphaPhi2BD, alphaPhi2, zeroField(), zeroField(), 1, 0);
        }

        // ---------- 组装有限通量 ----------
        alphaPhi1 = alphaPhi1BD + lam1*alphaPhi1;
        alphaPhi2 = alphaPhi2BD + lam2*alphaPhi2;

        // ---------- α1 方程（此处一次性加入体积分率源） ----------
        const volScalarField aGate = pos(alpha1 - scalar(1e-6));
        fvScalarMatrix alpha1Eqn( fvm::ddt(alpha1) + fvc::div(alphaPhi1) - aGate*vDotvAlphal );
        alpha1Eqn.solve();

        // ---------- α2 方程（α2↔α3 扩散；D23 可为 0） ----------
        volScalarField Dc23(D23*max(alpha3, scalar(0))*pos(alpha2));
        volScalarField Dc32(D23*max(alpha2, scalar(0))*pos(alpha3));

        fvScalarMatrix alpha2Eqn
        (
            fvm::ddt(alpha2)
          + fvc::div(alphaPhi2)
          - fvm::laplacian(Dc23 + Dc32, alpha2)
        );
        alpha2Eqn.solve();

        // ---------- 有界与闭合 ----------
        alpha1 = min(max(alpha1, scalar(0)), scalar(1));
        alpha2 = min(max(alpha2, scalar(0)), scalar(1) - alpha1);
        alpha3 = scalar(1) - alpha1 - alpha2;

        alpha1.correctBoundaryConditions();
        alpha2.correctBoundaryConditions();
        alpha3.correctBoundaryConditions();

        // ---------- 更新质量通量（只用对流通量；不掺矩阵通量） ----------
        rhoPhi =
              alphaPhi1*(rho1 - rho3)
            + alphaPhi2*(rho2 - rho3)
            + phi*rho3;
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << nl;

    Info<< "N2 phase volume fraction = "
        << alpha2.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
        << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
        << nl;
}

