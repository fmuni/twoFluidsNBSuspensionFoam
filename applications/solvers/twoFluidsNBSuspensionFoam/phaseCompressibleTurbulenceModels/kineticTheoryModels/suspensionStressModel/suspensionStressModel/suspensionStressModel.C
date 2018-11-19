/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "suspensionStressModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(suspensionStressModel, 0);

    defineRunTimeSelectionTable(suspensionStressModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::suspensionStressModel::suspensionStressModel
(
    const dictionary& dict
)
:
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::suspensionStressModel::~suspensionStressModel()
{}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::suspensionStressModel::Q
(
    const volVectorField& U,
    const vector& lambda
) const
{
    
    volVectorField dirU
    (
        U
        /
        (
            max(
                    mag(U),
                    dimensionedScalar
                    (
                        "V",
                        dimVelocity,
                        1e-18
                    )
               )
        )
    );

    volVectorField dirCurlU
    (
        fvc::curl(U)
        /
        (
            max(
                    mag(fvc::curl(U)),
                    dimensionedScalar
                    (
                        "curlV",
                        dimVelocity/dimLength,
                        1e-18
                    )
               )
        )
    );

    volVectorField dirGradU
    (
        dirU^dirCurlU
        /
        (
            max(
                    mag(dirU^dirCurlU),
                    dimensionedScalar
                    (
                        "curlV",
                        dimless,
                        1e-18
                    )
               )
        )
    );


    return symm(
      
         lambda.x()*(dirU*dirU)
       + lambda.y()*(dirGradU*dirGradU)
       + lambda.z()*(dirCurlU*dirCurlU)
    );
}

Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::suspensionStressModel::Q
(
    const volVectorField& U,
    const volScalarField& lambda1,
    const volScalarField& lambda2,
    const volScalarField& lambda3
) const
{
    volVectorField dirU
    (
        U
        /
        (
            max(
                    mag(U),
                    dimensionedScalar
                    (
                        "V",
                        dimVelocity,
                        1e-18
                    )
               )
        )
    );

    volVectorField dirCurlU
    (
        fvc::curl(U)
        /
        (
            max(
                    mag(fvc::curl(U)),
                    dimensionedScalar
                    (
                        "curlV",
                        dimVelocity/dimLength,
                        1e-18
                    )
               )
        )
    );

    volVectorField dirGradU
    (
        dirU^dirCurlU
        /
        (
            max(
                    mag(dirU^dirCurlU),
                    dimensionedScalar
                    (
                        "curlV",
                        dimless,
                        1e-18
                    )
               )
        )
    );

    return symm
    (
         lambda1*(dirU*dirU)
       + lambda2*(dirGradU*dirGradU)
       + lambda3*(dirCurlU*dirCurlU)
    );
}


// ************************************************************************* //
