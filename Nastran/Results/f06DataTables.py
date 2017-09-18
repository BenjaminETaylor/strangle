"""
===============================================================================
    StrEngL.Nastran.Results.f06DataTables
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

import numpy as np
from StrEngL.Nastran.Results import BaseClasses
from StrEngL.Nastran.Results import Collections
# f06 results page data for each type

#==============================================================================
# page start lines {'type': dataStartLine}
#==============================================================================
startLines = {
    'QUAD4_FORCES': 9,
    'TRIA3_FORCES': 8,
    'CBUSH_FORCES': 7,
    'CBAR_FORCES': 7,
    'CBEAM_FORCES': 7}

#==============================================================================   
# general function tools
#==============================================================================
def formatLine(line, offset=0, remove=None):
    """ formats lines for line-by-line result tables
        line: line to format
        offset: index to begin (everything before discarded)
        remove: optional single item to remove"""
    data = line[offset:].strip().split()
    if remove: data.pop(remove)
    ID = int(data.pop(0))
    data = list(map(lambda x: float(x), data))
    return ID, data
    
def groupLines(page, count, stride):
    """ extracts 'count' number of lines at each 'stride'
        if 'stride' is greater than 'count', lines are skipped"""
    pageLineGroups = []
    for i in range(count):
        pageLineGroups.append(page[i::stride])
    groups = []
    while pageLineGroups[0]:
        groups.append(
            [lineGroup.pop(0) for lineGroup in pageLineGroups])
    return groups

#==============================================================================
# parser functions create results objects from BaseClasses
#==============================================================================
def parseQuad4Forces(page):
    results = Collections.ElementResults()
    for line in page:      
        ID, data = formatLine(line)
        forces = np.array(data[:3])
        moments = np.array(data[3:6])
        shears = np.array(data[6:8])
        results[ID] = BaseClasses.Element2D(ID, forces, moments, shears)
    return results
    
def parseQuad4ForcesBilin(page):
    results = Collections.ElementResults()
    groups = groupLines(page, 1, 5)
    for group in groups:      
        ID, data = formatLine(group[0], offset=1, remove=1)
        forces = np.array(data[:3])
        moments = np.array(data[3:6])
        shears = np.array(data[6:8])
        results[ID] = BaseClasses.Element2D(ID, forces, moments, shears)
    return results
    
def parseTria3Forces(page):
    return parseQuad4Forces(page)
    
def parseCbushForces(page):
    results = Collections.ElementResults()
    for line in page:
        ID, data = formatLine(line,offset=1)
        forces = np.array(data[0:3])
        moments = np.array(data[3:6])
        results[ID] = BaseClasses.Element0D(ID, forces, moments)
    return results

def parseCbarForces(page):
    results = Collections.ElementResults()
    for line in page:
        ID, data = formatLine(line)
        forces = np.array([data[6]]+data[4:6])
        momentsA = np.array([data[7]]+data[0:2])
        momentsB = np.array([data[7]]+data[2:4])
        results[ID] = BaseClasses.Element1D(ID, forces, momentsA, momentsB)
    return results
        
def parseCbeamForces(page):
    results = Collections.ElementResults()
    groups = groupLines(page, 3, 3)
    for group in groups:
        ID = int(group[0][1:].strip())
        g1, data1 = formatLine(group[1])
        g2, data2 = formatLine(group[2])
        forces = np.array([data1[5]]+data1[3:5])
        momentsA = np.array([data1[6]]+data1[1:3])
        momentsB = np.array([data2[6]]+data2[1:3])
        results[ID] = BaseClasses.Element1D(ID, forces, momentsA, momentsB)
    return results
        
#==============================================================================
# parser tools for each results title
#==============================================================================
parserTools = {
    'QUAD4_FORCES': lambda page: parseQuad4Forces(page),
    'QUAD4_FORCES_BILIN': lambda page: parseQuad4ForcesBilin(page),
    'TRIA3_FORCES': lambda page: parseTria3Forces(page),
    'CBUSH_FORCES': lambda page: parseCbushForces(page),
    'CBAR_FORCES': lambda page: parseCbarForces(page),
    'CBEAM_FORCES': lambda page: parseCbeamForces(page)}