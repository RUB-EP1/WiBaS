/**************************************************************
 *                                                            *
 *  WiBaS                                                     *
 *                                                            *
 *  Williams' Background Suppression                          *
 *                                                            *
 *  Author: Julian Pychy                                      *
 *   email: julian@ep1.rub.de                                 *
 *                                                            *
 *  Copyright (C) 2016  Julian Pychy                          *
 *                                                            *
 *                                                            *
 *  Description:                                              *
 *                                                            *
 *  License:                                                  *
 *                                                            *
 *  This file is part of WiBaS                                *
 *                                                            *
 *  WiBaS is free software: you can redistribute it and/or    *
 *  modify it under the terms of the GNU General Public       *
 *  License as published by the Free Software Foundation,     *
 *  either version 3 of the License, or (at your option) any  *
 *  later version.                                            *
 *                                                            *
 *  WiBaS is distributed in the hope that it will be useful,  *
 *  but WITHOUT ANY WARRANTY; without even the implied        *
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR   *
 *  PURPOSE. See the GNU General Public License for more      *
 *  details.                                                  *
 *                                                            *
 *  You should have received a copy of the GNU General        *
 *  Public License along with WiBaS (license.txt). If not,    *
 *  see <http://www.gnu.org/licenses/>.                       *
 *                                                            *
 *************************************************************/


#include "PhasespaceCoord.hh"


PhasespaceCoord::PhasespaceCoord() :
    _id(0),
    _isCircular(false),
    _norm(1)
{
}



PhasespaceCoord::PhasespaceCoord(unsigned short int id, double norm, bool isCircular) :
    _id(id),
    _isCircular(isCircular),
    _norm(norm)
{
}



void PhasespaceCoord::SetNorm(double norm){
    _norm = norm;
}



unsigned short int PhasespaceCoord::GetID() const
{
    return _id;
}



bool PhasespaceCoord::GetIsCircular() const
{
    return _isCircular;
}



double PhasespaceCoord::GetNorm() const
{
    return _norm;
}
