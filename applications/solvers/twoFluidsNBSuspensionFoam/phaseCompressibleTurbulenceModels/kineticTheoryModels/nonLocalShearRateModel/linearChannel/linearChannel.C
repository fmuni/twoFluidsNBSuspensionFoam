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

#include "linearChannel.H"
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
    defineTypeNameAndDebug(linearChannel, 0);

    addToRunTimeSelectionTable
    (
        nonLocalShearRateModel,
        linearChannel,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonLocalShearRateModels::linearChannel::
linearChannel
(
    const dictionary& dict
)
:
    nonLocalShearRateModel(dict),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    H_
    (
        "channelHeight",
        dimLength,
        readScalar(coeffDict_.lookup("channelHeight"))
    ),
    suspensionBased_(
        coeffDict_.lookupOrDefault<Switch>("suspensionShear",false)
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonLocalShearRateModels::linearChannel::
~linearChannel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::nonLocalShearRateModels::linearChannel::
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
    volVectorField U = phase.U();

    if(suspensionBased_)
    {
        U = phase.fluid().U();
    }

    volSymmTensorField shear = symm(dev(fvc::grad(U)));

    SR = sqrt(scalar(2) * (shear&&shear));

    //- Get delta velocity over channel
    dimensionedScalar deltaU = max(mag(U)) - min(mag(U));

    //- Finally add the uniform constant value
    SR += phase.d()/H_ * (scalar(2)*deltaU/H_) ;

    return tSR;

}



bool Foam::kineticTheoryModels::nonLocalShearRateModels::linearChannel::
read()
{
    coeffDict_ <<= dict_.subDict(typeName + "Coeffs");

    H_.read(coeffDict_);

    return true;
}


// ************************************************************************* //
