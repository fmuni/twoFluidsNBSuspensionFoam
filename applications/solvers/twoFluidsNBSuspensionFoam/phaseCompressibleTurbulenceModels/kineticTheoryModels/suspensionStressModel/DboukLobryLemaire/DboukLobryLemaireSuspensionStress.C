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

#include "DboukLobryLemaireSuspensionStress.H"
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
    defineTypeNameAndDebug(DLL, 0);

    addToRunTimeSelectionTable
    (
        suspensionStressModel,
        DLL,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::suspensionStressModels::DLL::
DLL
(
    const dictionary& dict
)
:
    suspensionStressModel(dict),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    a_("a", dimless, coeffDict_),
    b_("b", dimless, coeffDict_),
    c_("c", dimless, coeffDict_),
    KN_("KN", dimless, coeffDict_),
    fixedLambdas_(coeffDict_.lookup("fixedLambdas"))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::suspensionStressModels::DLL::
~DLL()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::suspensionStressModels::DLL::computeQ
(
    const phaseModel& phase
) const
{
    const volScalarField& alpha = phase;
    scalar alphaMax(phase.alphaMax());
    volScalarField lambda1 = alpha*0 + scalar(1);
    volScalarField lambda2 = scalar(0.81)*alpha/alphaMax + scalar(0.66);
    volScalarField lambda3 = -scalar(0.0088)*alpha/alphaMax + scalar(0.54);

    if(fixedLambdas_)
    {
        lambda2 = scalar(0.8);
        lambda3 = scalar(0.5);
    }

    return Q(phase.U(),lambda1,lambda2,lambda3);
}

Foam::tmp<Foam::volVectorField>
Foam::kineticTheoryModels::suspensionStressModels::DLL::shearForce
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

        -KN_*invA*invA*muf/(phase.rho())
      * fvc::div(gammaDotEff*fabT)
    );
}

Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::suspensionStressModels::DLL::SigmaS
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

        -KN_*invA*invA*muf/(phase.rho())
      * gammaDotEff*fabT
    );
}

Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::suspensionStressModels::DLL::DSigmaDalpha
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

        KN_*scalar(2)*muf*(alphaMax*alpha)
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
Foam::kineticTheoryModels::suspensionStressModels::DLL::nu
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

    volScalarField alphaS
    (
        alphaMax
        /
        max(alpha,phase.residualAlpha())
    );

    return
    (
         (muf/phase.rho())
         *
         (
            a_
            +
            invA* b_*alphaMax
            +
            c_*alphaS*alphaS*invA*invA
         )

    );

}

bool Foam::kineticTheoryModels::suspensionStressModels::DLL::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    a_.read(coeffDict_);
    b_.read(coeffDict_);
    c_.read(coeffDict_);
    KN_.read(coeffDict_);


    return true;
}


// ************************************************************************* //
