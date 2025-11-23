// ---- TEqn.H (recommended form) ----

// conductivity on faces: k [W/(m·K)]
const surfaceScalarField kf = mixture->kappaf();

// volumetric heat capacity: rhoCp [J/(m3·K)]
volScalarField rhoCpV = max(
    alpha1*rho1*cp1 + alpha2*rho2*cp2 + alpha3*rho3*cp3,
    dimensionedScalar("rhoCpFloor", rho1.dimensions()*cp1.dimensions(), 1e-6)
);

// effective thermal diffusivity: D_T = k / (rhoCp) [m2/s]
const surfaceScalarField Deff = kf / fvc::interpolate(rhoCpV);

// pyrolysis rate (only in liquid), units W/m3 after multiplying Hr
const scalar alphaCut = 1e-6;
const Pair<tmp<volScalarField> > vDot = mixture->vDotAlphal();   // [1/s] after limiter
volScalarField Rr =
    rho1 * max(vDot[1](), dimensionedScalar("zero", dimless/dimTime, 0.0))
         * pos(alpha1 - alphaCut);  // kg/m3/s

// if vDot[1] was given in 1/min, uncomment the next line:
// Rr *= (1.0/60.0);

dimensionedScalar Hr("Hr", dimEnergy/dimMass, 5.414e5); // 541.4 kJ/kg → J/kg
/*volScalarField Qreact
(
    IOobject("Qreact", mesh.time().timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    -Hr * Rr   // W/m3
);                        */

//volScalarField Qreact = -Hr * Rr;   // 现在是临时量
Qreact = -Hr * Rr;

fvScalarMatrix TEqn
(
    fvm::ddt(T)
  + fvm::div(phi, T)
  - fvm::laplacian(Deff, T)
 == Qreact / rhoCpV
);

/*fvScalarMatrix TEqn
(
    fvm::ddt(rhoCpV, T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kf, T)
 ==
    Qreact
);*/

TEqn.relax();
fvOptions.constrain(TEqn);   // 确保 fvOptions 里没有 T 的体热源，否则会重复
TEqn.solve();
fvOptions.correct(T);
T.correctBoundaryConditions();

// ---- end TEqn.H ----


// （可选，仅在起步阶段短暂开启以防炸解，稳定后注释掉）
// T = min(max(T, dimensionedScalar("TminPhys", dimTemperature, 299.0)),
//         dimensionedScalar("TmaxPhys", dimTemperature, 1200.0));

// ===== end TEqn.H =====

