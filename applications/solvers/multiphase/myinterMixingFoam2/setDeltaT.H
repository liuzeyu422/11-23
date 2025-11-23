/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

{
    if (adjustTimeStep)
    {
        const dictionary& pDict = pimple.dict();

        const scalar maxCo      = pDict.lookupOrDefault<scalar>("maxCo", 0.9);
        const scalar maxAlphaCo = pDict.lookupOrDefault<scalar>("maxAlphaCo", 0.9);
        const scalar maxDeltaT  = pDict.lookupOrDefault<scalar>("maxDeltaT", GREAT);

        // CourantNo.H / alphaCourantNo.H 已经给出了 CoNum、alphaCoNum（标量）
        scalar deltaTFact = min
        (
            maxCo/(CoNum + SMALL),
            maxAlphaCo/(alphaCoNum + SMALL)
        );

        // 限制每步增大幅度（可选，避免放大过快）
        deltaTFact = max(deltaTFact, scalar(0.1));   // 允许减小到 0.1 倍
        deltaTFact = min(deltaTFact, scalar(1.2));   // 最多放大到 1.2 倍

        runTime.setDeltaT( min(deltaTFact*runTime.deltaTValue(), maxDeltaT) );
    }

    Info<< "deltaT = " << runTime.deltaTValue() << endl;
}


// ************************************************************************* //
