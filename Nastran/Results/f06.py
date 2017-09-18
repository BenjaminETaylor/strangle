"""
===============================================================================
    StrEngL.Nastran.Results.f06
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

import time
#from itertools import ifilter
from StrEngL.Nastran.Results import f06DataTables

class f06File():
    
    def __init__(self, filename):
        # instance variables
        self.filename = filename
        self.lines = open(filename,'r').readlines()
        self.pages = []
        #self._scanFile()
        
    #def __getattr__(self, name):
        # for list attribute calls
        #return getattr(self.lines, name)
    
    def __getitem__(self, key):
        # for slicing operations
        return self.lines[key]
            
    def getPage(self, pageNum):
        return self.pages[pageNum]
       
    def getTitles(self):
        titles = set()
        for page in self.pages:
            if page.title != None:
                titles.add(page.title)
        return titles
        
    def getElementResults(self, titles):
        """        
        Accepts: titles = python list of result titles
        Returns: Dictionary
           results[title][ subcase ][ elementID ]
        Assumes: subcases do not span multiple f06 files 
                 (standard NATRAN convention)
        """
        for title in titles:
            if title not in self.getTitles():
                print("f06 does not contain %s!" % title)
                titles.remove(title)
        # initialize dictionaries to contain all data and headers 
        results = {}
        startTime = time.time()
        # filter the file for pages with indicated result titles and 
        # iterate through each page storing the data
        for title in titles:
            data = {}
            for page in self._filterPages(title):
                # grab parser function from f06DataTables module
                if page.option: parseTitle = title + '_' + page.option
                else: parseTitle = title
                if not parseTitle in f06DataTables.parserTools:
                    print("Could not parse %s!" % parseTitle)
                    break
                parserFunction = f06DataTables.parserTools[parseTitle]
                # check if subcase has been found yet, if not add subcase key
                # to dictionaries
                subcase = page.subcase
                if subcase not in data: 
                    data[subcase] = []
                # collect the page data
                data[subcase].extend(page.getDataList())
            for subcase in data:
                # initialize results collection object
                data[subcase] = parserFunction(data[subcase])
            for subcase in data:
                if subcase not in results:
                    results[subcase] = {}
                results[subcase].update(data[subcase])
        print("%s loaded in %.2f seconds" % 
              (self.filename, time.time() - startTime))
        return results
        
    def getHash(self):
        # generates a hash ID for future new file checks
        import hashlib
        fileObj = open(self.filename, 'rb')
        hashType = hashlib.sha256()
        hashType.update(fileObj.read(1000*128)) # approx. 1st 1000 lines
        fileObj.close()
        return hashType.digest()
        
    def scanHeaders(self):
        # scan the file and create f06Page objects
        startTime = time.time()
        print("scanning %s..." % self.filename)
        # generate pages
        self._generatePages()
        for page in self.pages:
            page.scanHeader()
        print("took %.2f seconds" % (time.time() - startTime,))
        
    def openFile(self):
        self.lines = open(self.filename, 'r').readlines()
    
    def closeFile(self):
        self.lines = []
        
    def _generatePages(self):
        # generates f06Page objects, stores each instance in self.pages
        startLines = self._getStartLines()
        for i in range(len(startLines) - 1):
            startLine = startLines[i]
            endLine = startLines[i+1]
            self.pages.append(f06Page(i, startLine, endLine, self))
        
    def _getStartLines(self):
        # returns line numbers for all start lines in f06 file
        # ends with last line in file
        lineNumbers = range(len(self.lines))
        linesDict = dict(zip(self.lines, lineNumbers))
        startLines = []
        filterFunc = lambda line: line.startswith('1') 
        for line in filter(filterFunc, self):
            lineNum = linesDict[line]
            startLines.append(lineNum)
        startLines.append(len(self.lines))
        return startLines
                
    def _filterPages(self, title):
        # accepts 'title', returns pages of type 'title'
        filterFunc = lambda page: page.title==title
        return filter(filterFunc, self.pages)

        
class f06Page():
    # class variables
    READ_HEADER_LENGTH = 10 # number of lines read during header scan
    
    def __init__(self, number, start, end, f06FileObj):
        # instance variables
        self.number = number
        self.start = start
        self.end = end
        self.f06File = f06FileObj
        self.title = None
        self.subcase = None
        self.option = None
       
    def __len__(self):
        # returns page length in bytes
        return int(self.end - self.start)
    
    def __iter__(self):
        # for iteration (may be redundant)
        return iter(self.getDataList())
        
    def __getitem__(self, key):
        # for slicing operations
        return self.getDataList()[key]
    
    def scanHeader(self):
        import re
        endLine = self.start + self.READ_HEADER_LENGTH
        scanLines = self.f06File[ self.start:endLine ]
        for line in scanLines:
            if 'SUBCASE ' in line: # found subcase
                # capture the subcase and add it to current page
                self.subcase = self._captureSubCase(line)
            elif re.search(r'\w\s\w\s\w\s.*\(.*\)',line): # found title
                # parse and capture the title, then add it to current page
                self.title = self._parseTitle(line)
    
    def getDataList(self):
        # returns: list of data lines from f06 page (excluding header)
        return self.f06File[ self.getDataStartLine() : self.end ]

    def getHeader(self):
        # returns header of f06 page
        return self.f06File[ self.start : self.getDataStartLine() ]
        
    def getDataStartLine(self):
        # sets data start line for page
        return self.start + f06DataTables.startLines[self.title]
        
    def _parseTitle(self, line):
        # converts f06 title to key variable format,
        # '<ElementTYPE>_<resultsTYPE>'
        title = [string.replace(' ','')
                    for string in line.strip().split('  ')]
        # remove empty elements
        while '' in title:
            title.remove('')
        if any([')' in word for word in title]):
            # check for options (CORNER,CUBIC,SGAGE,BILIN)
            if 'OPTION' in title[-1]: 
                self.option = title.pop()[7:]
            title = title[-1]+title[0]
            title = title.strip('(').replace(')','_')
        else: title = '_'.join(title)
        return title
        
    def _captureSubCase(self, line):
        # strips subcase ID
        return int(line.strip().split('SUBCASE ')[-1])
    
    # future method to capture load steps in dynamic analyses    
    #def _captureLoadStep(self, line):
    #    return float(line.strip().split('=')[-1])