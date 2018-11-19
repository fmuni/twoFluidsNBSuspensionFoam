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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::kineticTheoryModels::suspensionStressModel>
Foam::kineticTheoryModels::suspensionStressModel::New
(
    const dictionary& dict
)
{
    word suspensionStressModelType(dict.lookup("suspensionStressModel"));

    Info<< "Selecting suspensionStressModel "
        << suspensionStressModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(suspensionStressModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "suspensionStressModel::New(const dictionary&) : " << endl
            << "    unknown suspensionStressModelType type "
            << suspensionStressModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid suspensionStressModelType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<suspensionStressModel>(cstrIter()(dict));
}


// ************************************************************************* //
