/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Author
    Federico Municchi, Purdue University (2018). All rights reserved.
\*---------------------------------------------------------------------------*/

#include "kineticShearRate.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvc.H"
#include "twoPhaseSystem.H"
#include "slipFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace nonLocalShearRateModels
{
    defineTypeNameAndDebug(kineticShearRate, 0);

    addToRunTimeSelectionTable
    (
        nonLocalShearRateModel,
        kineticShearRate,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonLocalShearRateModels::kineticShearRate::
kineticShearRate
(
    const dictionary& dict
)
:
    nonLocalShearRateModel(dict),
    coeffDict_(dict.subDict(typeName + "Coeffs"))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonLocalShearRateModels::kineticShearRate::
~kineticShearRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::nonLocalShearRateModels::kineticShearRate::
gammaDotEff
(
   const phaseModel& phase,
   const volScalarField& Theta
) const
{
    tmp<volScalarField> tSR
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("tSR", phase.name()),
                Theta.time().timeName(),
                Theta.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Theta.mesh(),
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        )
    );

    volScalarField& SR = tSR.ref();

    //- First eveluate local shear rate
    volVectorField U = phase.fluid().U();
    volSymmTensorField shear = symm(dev(fvc::grad(U)));

    //- Evaluate effective shear rate
    SR =
    (
        sqrt(scalar(2) * (shear&&shear))
      + (sqrt(scalar(2)*Theta)/phase.d())
    );

    return tSR;

}



bool Foam::kineticTheoryModels::nonLocalShearRateModels::kineticShearRate::
read()
{
    coeffDict_ <<= dict_.subDict(typeName + "Coeffs");


    return true;
}


// ************************************************************************* //
