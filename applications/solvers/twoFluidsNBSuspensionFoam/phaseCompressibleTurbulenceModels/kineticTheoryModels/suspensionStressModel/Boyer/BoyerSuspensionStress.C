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

#include "BoyerSuspensionStress.H"
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
namespace suspensionStressModels
{
    defineTypeNameAndDebug(Boyer, 0);

    addToRunTimeSelectionTable
    (
        suspensionStressModel,
        Boyer,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::suspensionStressModels::Boyer::
Boyer
(
    const dictionary& dict
)
:
    suspensionStressModel(dict),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    mu1_("mu1", dimless, coeffDict_),
    mu2_("mu2", dimless, coeffDict_),
    I0_("I0", dimless, coeffDict_),
    lambda_("lambdaAnisotropic", dimless, coeffDict_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::suspensionStressModels::Boyer::
~Boyer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::suspensionStressModels::Boyer::computeQ
(
    const phaseModel& phase
) const
{

    return Q(phase.U(),lambda_.value());
}

Foam::tmp<Foam::volVectorField>
Foam::kineticTheoryModels::suspensionStressModels::Boyer::shearForce
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMax,
    const volScalarField& gammaDotEff,
    const volSymmTensorField& fabT
) const
{
    //- Collect and calculate auxiliary quantities

    const volScalarField& alpha = phase;
    const volScalarField& muf   =
        phase.fluid().otherPhase(phase).mu();

    volScalarField invA
    (
        scalar(1)
        /
        (
            alphaMax
            /
            max(alpha,phase.residualAlpha())
            - scalar(1)
        )
    );

    return
    (

        -invA*invA*muf/(phase.rho())
      * fvc::div(gammaDotEff*fabT)
    );
}

Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::suspensionStressModels::Boyer::SigmaS
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMax,
    const volScalarField& gammaDotEff,
    const volSymmTensorField& fabT
) const
{
    //- Collect and calculate auxiliary quantities

    const volScalarField& alpha = phase;
    const volScalarField& muf   =
        phase.fluid().otherPhase(phase).mu();

    volScalarField invA
    (
        scalar(1)
        /
        (
            alphaMax
            /
            max(alpha,phase.residualAlpha())
            - scalar(1)
        )
    );

    return
    (

        -invA*invA*muf/(phase.rho())
      * gammaDotEff*fabT
    );
}


Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::suspensionStressModels::Boyer::DSigmaDalpha
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMax,
    const volScalarField& gammaDotEff,
    const volSymmTensorField& fabT
) const
{
    //- Collect and calculate auxiliary quantities
    const volScalarField& alpha = phase;
    const volScalarField& muf   =
        phase.fluid().otherPhase(phase).mu();


    return(

        scalar(2)*muf*(alphaMax*alpha)
      / (
            pow
            (
                max(alphaMax-alpha, phase.residualAlpha()),
                3
            )
        )
     *  (gammaDotEff*fabT)

    );

}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::suspensionStressModels::Boyer::nu
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMax
) const
{
    //- Auxiliary quantities
    const volScalarField& alpha = phase;
    const volScalarField& muf   =
        phase.fluid().otherPhase(phase).mu();

    volScalarField invA
    (
        scalar(1)
        /
        (
            alphaMax
            /
            max(alpha,phase.residualAlpha())
            - scalar(1)
        )
    );

    return
    (
         (muf/phase.rho())
         *
         (
            invA* scalar(5)/scalar(2)*alphaMax
            +
            invA*invA*(mu1_ + (mu2_ - mu1_)/(scalar(1)+I0_*invA*invA ))
         )

    );

}

bool Foam::kineticTheoryModels::suspensionStressModels::Boyer::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    mu1_.read(coeffDict_);
    mu2_.read(coeffDict_);
    I0_.read(coeffDict_);

    lambda_.read(coeffDict_);

    return true;
}


// ************************************************************************* //
