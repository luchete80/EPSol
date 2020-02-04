/************************************************************************
	
	Copyright 2012-2013 Luciano Buglioni
 
	Contact: luciano.buglioni@gmail.com

	This file is a part of FluxSol

	FluxSol is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    Free CFD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/
#ifndef products_H
#define products_H

#include "pTraits.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FluxSol
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<int rank>
class typeOfRank
{};

//El arg es el tipo al que le aplico el gradiente
//Es como el openfoam pero esta cambiada
template<class arg1,class arg2>
class outerProduct
{
public:

    typedef typename typeOfRank
    <
        int(pTraits<arg1>::rank)+ int(pTraits<arg2>::rank)
    >::type type;
};


//El arg es el tipo al que le aplico el gradiente
//Es como el openfoam pero esta cambiada
template<class arg1,class arg2>
class innerProduct
{
public:

    typedef typename typeOfRank
    <
        int(pTraits<arg1>::rank)+ int(pTraits<arg2>::rank) - 2
    >::type type;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
