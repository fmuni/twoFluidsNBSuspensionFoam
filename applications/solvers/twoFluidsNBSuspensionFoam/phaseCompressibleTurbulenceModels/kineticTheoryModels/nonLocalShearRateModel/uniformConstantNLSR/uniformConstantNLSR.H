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

Class
    Foam::kineticTheoryModels::frictionalStressModels::Boyer

Description
    Uniform constant non local shear rate

SourceFiles
    uniformConstant.C

Author
    Federico Municchi, Purdue University (2018). All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef uniformConstantNLSR_H
#define uniformConstantNLSR_H

#include "nonLocalShearRateModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace nonLocalShearRateModels
{

class uniformConstantNLSR
:
    public nonLocalShearRateModel
{
    // Private data

        dictionary coeffDict_;

        dimensionedScalar gammaDotNL_;

        Switch  suspensionBased_;

public:

    //- Runtime type information
    TypeName("uniformConstant");


    // Constructors

        //- Construct from components
        uniformConstantNLSR(const dictionary& dict);


    //- Destructor
    virtual ~uniformConstantNLSR();


    // Member functions

        virtual tmp<volScalarField> gammaDotEff
        (
            const phaseModel& phase,
            const volScalarField& Theta
        ) const ;

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionalStressModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
