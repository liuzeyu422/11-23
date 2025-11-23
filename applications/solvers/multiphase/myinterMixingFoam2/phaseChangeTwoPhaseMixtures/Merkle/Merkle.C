/*---------------------------------------------------------------------------*\
| Merkle.C
\*---------------------------------------------------------------------------*/

#include "Merkle.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Merkle, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, Merkle, components);
}
}

// ===== ctor =====
Foam::phaseChangeTwoPhaseMixtures::Merkle::Merkle
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),
    preExp_( phaseChangeTwoPhaseMixtureCoeffs_.lookup("preExp") ),   // [1/s]
    EoverR_( phaseChangeTwoPhaseMixtureCoeffs_.lookup("EoverR") ),    // [K]
    // defaults (overridden below if found in dict)
    TminArr_("TminArr", dimTemperature, 753.0),
    TmaxArr_("TmaxArr", dimTemperature, 1200.0),
    TexpMin_("TexpMin", dimTemperature, 300.0),
    lnExpMax_(50.0)   // exp(50) ~ 3e21
{
    // optional guards as FOUR-PART (dimensioned) entries are supported
    if (phaseChangeTwoPhaseMixtureCoeffs_.found("TminArr"))
    {
        TminArr_ = dimensionedScalar
        (
            phaseChangeTwoPhaseMixtureCoeffs_.lookup("TminArr")
        );
    }
    if (phaseChangeTwoPhaseMixtureCoeffs_.found("TmaxArr"))
    {
        TmaxArr_ = dimensionedScalar
        (
            phaseChangeTwoPhaseMixtureCoeffs_.lookup("TmaxArr")
        );
    }
    if (phaseChangeTwoPhaseMixtureCoeffs_.found("TexpMin"))
    {
        TexpMin_ = dimensionedScalar
        (
            phaseChangeTwoPhaseMixtureCoeffs_.lookup("TexpMin")
        );
    }
    if (phaseChangeTwoPhaseMixtureCoeffs_.found("lnExpMax"))
    {
        // allow lnExpMax as FOUR-PART (dimless) entry
        dimensionedScalar lnE
        (
            phaseChangeTwoPhaseMixtureCoeffs_.lookup("lnExpMax")
        );
        lnExpMax_ = lnE.value();  // store as plain scalar
    }

    correct();
}

// ===== mDotAlphal: return (condensation, vaporisation) in kg/m3/s =====
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Merkle::mDotAlphal() const
{
    const volScalarField& T =
        alpha1_.db().lookupObject<volScalarField>("T");

    // condensation = 0
    tmp<volScalarField> mDotCond
    (
        new volScalarField
        (
            IOobject("mDotCond", T.time().timeName(), T.mesh(),
                     IOobject::NO_READ, IOobject::NO_WRITE),
            T.mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    // guards for Arrhenius
    volScalarField Tden = max(T, TexpMin_);     // keep denominator > 0
    volScalarField Tcap = min(Tden, TmaxArr_);  // upper cap
    tmp<volScalarField> gateTmp = pos(T - TminArr_);  // below Tmin -> 0
    const volScalarField& gate = gateTmp();

    // 指数幅值上限：argPos = min(E/R/T, lnExpMax)
    const volScalarField argPos =
        min(EoverR_/Tcap, dimensionedScalar("lnExpMax", dimless, lnExpMax_));

    // mDot_vap [kg/m3/s] = ρ_l * A * exp(-argPos) * gate * 1[α1 > aCut]
    const dimensionedScalar aCut("aCut", dimless, 1e-6);

    tmp<volScalarField> mDotVap
    (
        new volScalarField
        (
            IOobject("mDotVap", T.time().timeName(), T.mesh(),
                     IOobject::NO_READ, IOobject::NO_WRITE),
            rho1()*preExp_*exp(-argPos)*gate*pos(alpha1_ - aCut)
        )
    );

    return Pair<tmp<volScalarField> >(mDotCond, mDotVap);
}

// ===== pressure-driven part not used =====
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Merkle::mDotP() const
{
    const fvMesh& mesh = alpha1_.mesh();

    tmp<volScalarField> zero1
    (
        new volScalarField
        (
            IOobject("mDotP1", mesh.time().timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    tmp<volScalarField> zero2
    (
        new volScalarField
        (
            IOobject("mDotP2", mesh.time().timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    return Pair<tmp<volScalarField> >(zero1, zero2);
}

void Foam::phaseChangeTwoPhaseMixtures::Merkle::correct(){}

bool Foam::phaseChangeTwoPhaseMixtures::Merkle::read()
{
    if (!phaseChangeTwoPhaseMixture::read()) return false;

    phaseChangeTwoPhaseMixtureCoeffs_.lookup("preExp")  >> preExp_;
    phaseChangeTwoPhaseMixtureCoeffs_.lookup("EoverR")  >> EoverR_;

    if (phaseChangeTwoPhaseMixtureCoeffs_.found("TminArr"))
    {
        TminArr_ = dimensionedScalar
        (
            phaseChangeTwoPhaseMixtureCoeffs_.lookup("TminArr")
        );
    }
    if (phaseChangeTwoPhaseMixtureCoeffs_.found("TmaxArr"))
    {
        TmaxArr_ = dimensionedScalar
        (
            phaseChangeTwoPhaseMixtureCoeffs_.lookup("TmaxArr")
        );
    }
    if (phaseChangeTwoPhaseMixtureCoeffs_.found("TexpMin"))
    {
        TexpMin_ = dimensionedScalar
        (
            phaseChangeTwoPhaseMixtureCoeffs_.lookup("TexpMin")
        );
    }
    if (phaseChangeTwoPhaseMixtureCoeffs_.found("lnExpMax"))
    {
        dimensionedScalar lnE
        (
            phaseChangeTwoPhaseMixtureCoeffs_.lookup("lnExpMax")
        );
        lnExpMax_ = lnE.value();
    }
    return true;
}
// ************************************************************************* //

