{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    // ★ 防 1/0：所有除以 rAUf 的地方，统一用带护栏的 rAUfPos
    const surfaceScalarField rAUfPos
    (
        "rAUfPos",
        max(rAUf, dimensionedScalar("minrAUf", rAUf.dimensions(), SMALL))
    );

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p_rgh);

    // 表面张力+重力产生的“源通量”
    surfaceScalarField phig
    (
        (
            interface.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUfPos*mesh.magSf()
    );

    
    // 源通量合力一致：phig
    phiHbyA += phig;
    

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - MRF.relative(mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUfPos.boundaryField())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUfPos, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            // 最终速度修正：压力修正 + 体力Fdrag（向量，不丢切向）
            U = HbyA
              + rAU*(
                     fvc::reconstruct((phi - phiHbyA)/rAUfPos)
                  + Fdrag
                );

            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    p = p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}

