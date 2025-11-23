Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading volume\n" << endl;
volScalarField volume
(
    IOobject
    (
        "volume",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroV", dimVolume, 0.0)   // m^3
);

volume.primitiveFieldRef() = mesh.V();
volume.correctBoundaryConditions();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
autoPtr<phaseChangeTwoPhaseMixture> mixture =
    phaseChangeTwoPhaseMixture::New(U, phi);

volScalarField& alpha1(mixture->alpha1());
volScalarField& alpha2(mixture->alpha2());
volScalarField& alpha3(mixture->alpha3());

const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();
const dimensionedScalar& rho3 = mixture->rho3();

const dimensionedScalar& cp1 = mixture->cp1();
const dimensionedScalar& cp2 = mixture->cp2();
const dimensionedScalar& cp3 = mixture->cp3();

dimensionedScalar D23(mixture->lookup("D23"));

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3,
    alpha1.boundaryField().types()
);
rho.oldTime();

Info<< "Reading / calculating rho*cp\n" << endl;
// Need to store rhoCp for ddt(rhoCp, T)
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1*cp1 + alpha2*rho2*cp2 + alpha3*rho3*cp3,
    alpha1.boundaryField().types()
);
rhoCp.oldTime();

Info<< "Reading / calculating rho*phi*cp\n" << endl;
surfaceScalarField rhoCpPhi
(
    IOobject
    (
        "rhoCpPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);

Info<< "Creating reaction heat source field Qreact\n" << endl;
volScalarField Qreact
(
    IOobject
    (
        "Qreact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE   // 需要写盘就留着；不需要写盘也可用 NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroQ", dimEnergy/dimTime/dimVolume, 0.0)  // W/m^3
);

// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho1*phi
);

// Construct interface from alpha1/2/3
threePhaseInterfaceProperties interface(alpha1, alpha2, alpha3, U, mixture());

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "nearWallDist.H"

// 只创建一次 nearWallDist，UEqn.H 中直接复用 yDist.y()
nearWallDist yDist(mesh);

// 新增：cell 中心的壁距（专给 δ 用，避免在 PIMPLE 里临时构建）
#include "wallDist.H"
// 找到侧壁补丁
const polyBoundaryMesh& pbm = mesh.boundaryMesh();
label heatedID = pbm.findPatchID("heatedWall");
if (heatedID == -1)
{
    FatalErrorInFunction
        << "Patch 'heatedWall' not found in constant/polyMesh/boundary"
        << exit(FatalError);
}

// 只包含 heatedWall 的集合
labelHashSet sideOnly(1);
sideOnly.insert(heatedID);

// 第三个参数给 "wall"，wallDist 会去 fvSchemes 里找 "wallDist { ... }"
wallDist yToHeated(mesh, sideOnly, "wall");

// 发布一个可写的体场 y（到 heatedWall 的几何距离，单位 m）
volScalarField y
(
    IOobject("y", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    yToHeated.y()
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

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

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());
mesh.setFluxRequired(alpha2.name());
mesh.setFluxRequired(alpha3.name());   // 建议也声明

