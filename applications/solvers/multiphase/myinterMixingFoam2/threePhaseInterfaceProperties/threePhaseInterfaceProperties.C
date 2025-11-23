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

#include "threePhaseInterfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::threePhaseInterfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::threePhaseInterfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& alpha1 = alpha1_.boundaryField();
    const volScalarField::Boundary& alpha2 = alpha2_.boundaryField();
    const volScalarField::Boundary& alpha3 = alpha3_.boundaryField();
    const volVectorField::Boundary& U = U_.boundaryField();

    const fvMesh& mesh = U_.mesh();
    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha1[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& a2cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha2[patchi]);

            const alphaContactAngleFvPatchScalarField& a3cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha3[patchi]);

            scalarField twoPhaseAlpha2(max(a2cap, scalar(0)));
            scalarField twoPhaseAlpha3(max(a3cap, scalar(0)));

            scalarField sumTwoPhaseAlpha
            (
                twoPhaseAlpha2 + twoPhaseAlpha3 + SMALL
            );

            twoPhaseAlpha2 /= sumTwoPhaseAlpha;
            twoPhaseAlpha3 /= sumTwoPhaseAlpha;

            fvsPatchVectorField& nHatp = nHatb[patchi];

            scalarField theta
            (
                convertToRad
              * (
                   twoPhaseAlpha2*(180 - a2cap.theta(U[patchi], nHatp))
                 + twoPhaseAlpha3*(180 - a3cap.theta(U[patchi], nHatp))
                )
            );

            vectorField nf(boundary[patchi].nf());

            // Reset nHatPatch to correspond to the contact angle

            scalarField a12(nHatp & nf);

            scalarField b1(cos(theta));

            scalarField b2(nHatp.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            
            nHatp /= (mag(nHatp) + SMALL);// 用无量纲极小值归一，避免量纲混合

            //nHatp /= (mag(nHatp) + deltaN_.value());
        }
    }
}


void Foam::threePhaseInterfaceProperties::calculateK()
{
    const volScalarField& alpha1 = alpha1_;

    const fvMesh& mesh = alpha1.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    volVectorField gradAlpha(fvc::grad(alpha1));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    correctContactAngle(nHatfv.boundaryFieldRef());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    // volVectorField nHat = gradAlpha/(mag(gradAlpha) + deltaN_);
    // nHat.boundaryField() = nHatfv.boundaryField();
    // K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threePhaseInterfaceProperties::threePhaseInterfaceProperties
(
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const volScalarField& alpha3,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        readScalar
        (
            U.mesh().solverDict
            (
                alpha1.name()
            ).lookup("cAlpha")
        )
    ),
    sigma12_("sigma12", dimensionSet(1, 0, -2, 0, 0), dict),
    sigma13_("sigma13", dimensionSet(1, 0, -2, 0, 0), dict),

    deltaN_
    (
        "deltaN",
        dimless/dimLength,                               // [1/Length]
        1e-8
        //scalar(1e-8)/pow(average(U.mesh().V()).value(), 1.0/3.0)   // 纯 double
    ),

    alpha1_(alpha1),
    alpha2_(alpha2),
    alpha3_(alpha3),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::threePhaseInterfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::threePhaseInterfaceProperties::nearInterface() const
{
    return max
    (
        pos(alpha1_ - 0.01)*pos(0.99 - alpha1_),
        pos(alpha2_ - 0.01)*pos(0.99 - alpha2_)
    );
}


// ************************************************************************* //
