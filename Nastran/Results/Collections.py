"""
===============================================================================
    StrEngL.Nastran.Results.Collections
    Copyright (C) 2015  Benjamin E. Taylor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
===============================================================================
"""

class Results(dict):
    def __init__(self, *args):
        dict.__init__(self, *args)
    def sortXForces(self):
        raise NotImplementedError
    def sortYForces(self):
        raise NotImplementedError
    def sortZForces(self):
        raise NotImplementedError
    def sortXMoments(self):
        raise NotImplementedError
    def sortYMoments(self):
        raise NotImplementedError
    def sortZMoments(self):
        raise NotImplementedError

class ElementResults(Results):
    def __init__(self, *args):
        Results.__init__(self, *args)
    def sortStresses(self):
        raise NotImplementedError
    def sortStrains(self):
        raise NotImplementedError
        
class NodalResults(Results):
    def __init__(self, *args):
        Results.__init__(self, *args)
    def sortXDisplacements(self):
        raise NotImplementedError
    def sortYDisplacements(self):
        raise NotImplementedError
    def sortZDisplacements(self):
        raise NotImplementedError
    def sortXRotations(self):
        raise NotImplementedError
    def sortYRotations(self):
        raise NotImplementedError
    def sortZRotations(self):
        raise NotImplementedError