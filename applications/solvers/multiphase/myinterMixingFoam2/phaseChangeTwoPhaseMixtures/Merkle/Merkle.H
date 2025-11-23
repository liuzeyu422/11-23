/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::phaseChangeTwoPhaseMixtures::Merkle

Description
    Merkle-type pyrolysis/cavitation model rewritten with robust Arrhenius
    protection (temperature clipping and exponent limiting).

SourceFiles
    Merkle.C
\*---------------------------------------------------------------------------*/

#ifndef Merkle_H
#define Merkle_H

#include "phaseChangeTwoPhaseMixture/phaseChangeTwoPhaseMixture.H"

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

class Merkle : public phaseChangeTwoPhaseMixture
{
    // Kinetics parameters
    dimensionedScalar preExp_;     // A  [s^-1]
    dimensionedScalar EoverR_;     // E/R [K]

    // Numerical/physical guards for Arrhenius term
    dimensionedScalar TminArr_;    // below this, reaction off [K]
    dimensionedScalar TmaxArr_;    // cap Arrhenius evaluation [K]
    dimensionedScalar TexpMin_;    // min T used in exponent  [K]
    scalar            lnExpMax_;   // max exponent for exp(), dimensionless

public:
    TypeName("Merkle");

    Merkle
    (
        const volVectorField& U,
        const surfaceScalarField& phi
    );

    virtual ~Merkle() {}

    // Return pair of "coefficient-type" mass sources (kg/m3/s), not multiplied by alpha1
    virtual Pair<tmp<volScalarField> > mDotAlphal() const;

    // No pressure-driven part in this model
    virtual Pair<tmp<volScalarField> > mDotP() const;

    virtual void correct();
    virtual bool read();
};

} // namespace phaseChangeTwoPhaseMixtures
} // namespace Foam

#endif

// ************************************************************************* //

