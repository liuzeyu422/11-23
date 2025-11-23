/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(phaseChangeTwoPhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixture::phaseChangeTwoPhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    incompressibleThreePhaseMixture(U, phi),
    phaseChangeTwoPhaseMixtureCoeffs_(subDict(type + "Coeffs")),
    pSat_("pSat", dimPressure, lookup("pSat")),
    // 置为液相密度，使 vDot = (alpha1/rho0)*mDotAlphal 成为 [1/s]
    rho0_(rho1())
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotAlphal() const
{
    const fvMesh& mesh = alpha1_.mesh();

    // 每步最多转化的体积分数比例 Ck（默认 0.10；保持你当前取值）
    const scalar Ck =
        phaseChangeTwoPhaseMixtureCoeffs_.lookupOrDefault<scalar>
        ("vDotLimiterCoeff", 0.05);

    // kMax = Ck/Δt  [1/s]
    const dimensionedScalar kMax
    (
        "kMax", dimless/dimTime, Ck/mesh.time().deltaTValue()
    );

    // ① 取“质量源”mDot（子类提供），单位 [kg/m3/s]
    Pair<tmp<volScalarField> > mDot = this->mDotAlphal();

    // ② α1 裁剪到 [0,1]
    tmp<volScalarField> a1clipTmp = min(max(alpha1_, scalar(0)), scalar(1));
    const volScalarField& a1clip = a1clipTmp();

    // ③ 体积分率蒸发源（带 Δt 限幅，且随 α1 缩放）：[1/s]
    //    vDot_vap = min( (α1/rho0)*mDot_vap , α1*kMax )
    tmp<volScalarField> vDotVapLimited
    (
        min
        (
            (a1clip/rho0_) * mDot[1],
            a1clip * kMax
        )
    );

    // ④ 凝结端置零场（[1/s]），避免对 tmp<const-ref> 赋值
    tmp<volScalarField> vDotCond
    (
        new volScalarField
        (
            IOobject
            (
                "vDotCond",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        )
    );

    return Pair<tmp<volScalarField> >(vDotCond, vDotVapLimited);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotP() const
{
    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
    Pair<tmp<volScalarField> > mDotP = this->mDotP();

    return Pair<tmp<volScalarField> >(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}

bool Foam::phaseChangeTwoPhaseMixture::read()
{
    if (incompressibleThreePhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        lookup("pSat") >> pSat_;
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //

