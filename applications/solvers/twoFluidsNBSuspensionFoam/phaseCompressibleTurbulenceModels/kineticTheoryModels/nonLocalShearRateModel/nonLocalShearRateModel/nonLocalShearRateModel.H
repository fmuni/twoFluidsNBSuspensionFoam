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

Class
    Foam::kineticTheoryModels::nonLocalShearRateModel

SourceFiles
    nonLocalShearRateModel.C

\*---------------------------------------------------------------------------*/

#ifndef nonLocalShearRateModel_H
#define nonLocalShearRateModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class nonLocalShearRateModel Declaration
\*---------------------------------------------------------------------------*/

class nonLocalShearRateModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        nonLocalShearRateModel(const nonLocalShearRateModel&);

        //- Disallow default bitwise assignment
        void operator=(const nonLocalShearRateModel&);


protected:

    // Protected data

        //- Reference to higher-level dictionary for re-read
        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("nonLocalShearRateModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        nonLocalShearRateModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        nonLocalShearRateModel(const dictionary& dict);


    // Selectors

        static autoPtr<nonLocalShearRateModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~nonLocalShearRateModel();


    // Member Functions


        virtual tmp<volScalarField> gammaDotEff
        (
            const phaseModel& phase,
            const volScalarField& Theta
        ) const = 0;


        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
