"""
===============================================================================
    StrEngL.Nastran.Results.f06Databases
    Copyright (C) 2016  Benjamin E. Taylor

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

import shelve
from StrEngL.Nastran.Results.f06 import f06File

class f06Db():
    """
    This class creates objects that will control and monitor groups of f06
    files. Its attributes contain information about the files it monitors.
    Its methods run checks on those files, remove them, add them, and update 
    them if necessary. Moreover, some methods automate the action of 
    performing the same task on all files in the database. The database
    object can also save itself to disk.
    """
    def __init__(self, f06FilePaths=[]):
        # instance variables
        self.files = {}
        self.filename = None
        # f06FilePaths is List of f06 file paths
        self._buildCollection(f06FilePaths)
        
    def _buildCollection(self, f06FilePaths):
        # collects f06File objects in self
        # keyed by hash ID
        for filePath in f06FilePaths:
            f06 = f06File(filePath)
            f06.scanHeaders()
            f06.closeFile()
            self.files[f06.getHash()] = f06
            
    def getFileList(self):
        # returns a list with all the filenames in the database
        return self.files.values() #might have to be list comprehension
            
    def removeFile(self, filename):
        # removes item based on value
        for (hashKey, f06) in self.files.items():
            if f06.filename == filename:
                del self.files[hashKey] #change to pop()
                break
            
    def addFile(self, filename):
        # adds filename with associated hash ID as key
        f06 = f06File(filename)
        f06.scanHeaders()
        f06.closeFile()
        self.files[f06.getHash()] = f06
                    
    def checkForUpdates(self):
        # checks each file for same hash ID
        # if hash ID is different, the file is updated
        badFiles = []
        newFiles = []
        for (hashKey, f06) in self.files.items():
            try: f06Temp = f06File(f06.filename)
            except: 
                raise IOError("could not find: %s" % f06.filename)
                # highlight file row
                badFiles.append(f06.filename)
                continue
            if f06Temp.getHash() == hashKey: pass
            else: 
                newFiles.append(f06.filename)
                badFiles.append(hashKey)
        for filename in badFiles:
            self.removeFile(filename)
        for filename in newFiles:
            self.addFile(filename)
            print("updated: %s" % filename)
            
    def getAllElementResults(self, titles):
        """ 
        Returns results for each title in titles found in all files 
        of self. Assumes each f06 file contains a set of unique subcase IDs..
        """
        import time
        allResults = {}
        start = time.time()
        for f06 in self.getFileList():
            f06.openFile()
            results = f06.getElementResults(titles)
            allResults.update(results)
            f06.closeFile()
        print('All results read in %.2f seconds.' % (time.time() - start,))
        return allResults
        
    def save(self, filename=None):
        # saves self to disk at filename
        # if filename not given, assume a db has been loaded and
        # self.filename is not None
        if not filename: filename = self.filename
        db = shelve.open(filename)
        db.clear()
        db.update(self.files)
        db.close()
        print('%i files written to database' % (len(self.files),))
        
    def loadDb(self, filename):
        #erase any existing files
        #need to prompt user to save existing db if status is modified
        self.files = {}
        self.filename = filename
        self.files.update(shelve.open(filename))
