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

#include "incompressibleThreePhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void Foam::incompressibleThreePhaseMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();
    nuModel3_->correct();

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(alpha1_*rho1_ + alpha2_*rho2_ + alpha3_*rho3_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleThreePhaseMixture::incompressibleThreePhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phase1Name_(wordList(lookup("phases"))[0]),
    phase2Name_(wordList(lookup("phases"))[1]),
    phase3Name_(wordList(lookup("phases"))[2]),

    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha3_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase3Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U.time().timeName(),
            U.db()
        ),
        U.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    nuModel3_
    (
        viscosityModel::New
        (
            "nu3",
            subDict(phase3Name_),
            U,
            phi
        )
    ),

    rho1_("rho", dimDensity, nuModel1_->viscosityProperties()),
    rho2_("rho", dimDensity, nuModel2_->viscosityProperties()),
    rho3_("rho", dimDensity, nuModel3_->viscosityProperties()),
    cp1_("cp", dimensionSet(0, 2, -2, -1, 0, 0, 0), nuModel1_->viscosityProperties().lookup("cp")),
    cp2_("cp", dimensionSet(0, 2, -2, -1, 0, 0, 0), nuModel2_->viscosityProperties().lookup("cp")),
    cp3_("cp", dimensionSet(0, 2, -2, -1, 0, 0, 0), nuModel3_->viscosityProperties().lookup("cp")),
    Pr1_("Pr", dimensionSet(0, 0, 0, 0, 0, 0, 0), nuModel1_->viscosityProperties().lookup("Pr")),
    Pr2_("Pr", dimensionSet(0, 0, 0, 0, 0, 0, 0), nuModel2_->viscosityProperties().lookup("Pr")),
    Pr3_("Pr", dimensionSet(0, 0, 0, 0, 0, 0, 0), nuModel3_->viscosityProperties().lookup("Pr"))
{
    alpha3_ = 1.0 - alpha1_ - alpha2_;
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::mu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            alpha1_*rho1_*nuModel1_->nu()
          + alpha2_*rho2_*nuModel2_->nu()
          + alpha3_*rho3_*nuModel3_->nu()
        )
    );
}

Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::muf() const
{
    const surfaceScalarField a1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );
    const surfaceScalarField a2f
    (
        min(max(fvc::interpolate(alpha2_), scalar(0)), scalar(1))
    );

    const surfaceScalarField sum12 = a1f + a2f;
    const surfaceScalarField scale =
        min(scalar(1), scalar(1)/max(sum12, scalar(SMALL)));

    const surfaceScalarField a1fn = a1f*scale;
    const surfaceScalarField a2fn = a2f*scale;
    const surfaceScalarField a3fn = scalar(1) - a1fn - a2fn;

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "mu",
            a1fn*rho1_*fvc::interpolate(nuModel1_->nu())
          + a2fn*rho2_*fvc::interpolate(nuModel2_->nu())
          + a3fn*rho3_*fvc::interpolate(nuModel3_->nu())
        )
    );
}

Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::nuf() const
{
    const surfaceScalarField a1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );
    const surfaceScalarField a2f
    (
        min(max(fvc::interpolate(alpha2_), scalar(0)), scalar(1))
    );
    const surfaceScalarField sum12 = a1f + a2f;
    const surfaceScalarField scale =
        min(scalar(1), scalar(1)/max(sum12, scalar(SMALL)));
    const surfaceScalarField a1fn = a1f*scale;
    const surfaceScalarField a2fn = a2f*scale;
    const surfaceScalarField a3fn = scalar(1) - a1fn - a2fn;

    const surfaceScalarField muface =
        a1fn*rho1_*fvc::interpolate(nuModel1_->nu())
      + a2fn*rho2_*fvc::interpolate(nuModel2_->nu())
      + a3fn*rho3_*fvc::interpolate(nuModel3_->nu());

    const surfaceScalarField rhoface =
        a1fn*rho1_ + a2fn*rho2_ + a3fn*rho3_;

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nu",
            muface/max
            (
                rhoface,
                dimensionedScalar("tinyRho", rho1_.dimensions(), SMALL)
            )
        )
    );
}


//原始版本1.0
/*Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::muf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "mu",
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
          + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::nuf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nu",
            (
                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
              + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
              + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
            )/(alpha1f*rho1_ + alpha2f*rho2_ + alpha3f*rho3_)
        )
    );
}*/

Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::kappaf() const
{
	const surfaceScalarField alpha1f
	(
		min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
	);

        const surfaceScalarField alpha2f
	(
		min(max(fvc::interpolate(alpha2_), scalar(0)), scalar(1))
	);
	
	const surfaceScalarField alpha3f
        (
            // 用“剩余量”并钳为非负：a3f = max(0, 1 - a1f - a2f)
            max(scalar(1) - alpha1f - alpha2f, scalar(0))
        );

	return tmp<surfaceScalarField>
	(
		new surfaceScalarField
		(
			"kappaf",
			(
				alpha1f*rho1_*cp1_*(1/Pr1_)
				*fvc::interpolate(nuModel1_->nu())
				+ alpha2f*rho2_*cp2_*(1/Pr2_)
                                *fvc::interpolate(nuModel2_->nu())
                                + alpha3f*rho3_*cp3_*(1/Pr3_)
                                *fvc::interpolate(nuModel3_->nu())
                                /*+ (scalar(1) - alpha1f-alpha2f)*rho3_*cp3_
				*(1/Pr3_)*fvc::interpolate(nuModel3_->nu())*/
			)
		)
	);
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::nu1() const
{
    // 直接转发相1(聚合物)模型给出的 kinematic ν
    return nuModel1_->nu();
}

Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::nu1f() const
{
    return fvc::interpolate(nuModel1_->nu());
}


bool Foam::incompressibleThreePhaseMixture::read()
{
    if (transportModel::read())
    {
        if
        (
            nuModel1_().read(*this)
         && nuModel2_().read(*this)
         && nuModel3_().read(*this)
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;
            nuModel3_->viscosityProperties().lookup("rho") >> rho3_;
            nuModel1_->viscosityProperties().lookup("cp") >> cp1_;
	    nuModel2_->viscosityProperties().lookup("cp") >> cp2_;
	    nuModel3_->viscosityProperties().lookup("cp") >> cp3_;
            nuModel1_->viscosityProperties().lookup("Pr") >> Pr1_;
	    nuModel2_->viscosityProperties().lookup("Pr") >> Pr2_;
	    nuModel3_->viscosityProperties().lookup("Pr") >> Pr3_;            
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
