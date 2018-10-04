# -*- coding: utf-8 -*-
"""
Created on Mon Oct 03 17:58:12 2016

@author: Ben
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import vtk
from vtk import vtkQuad

class MeshViewer(object):
    
    def __init__(self, mesh):
        self.mesh = mesh
        
    def render(self):
        self._generate_points()
        self._build_quads()
        self._build_viewing_heirarchy()
        self.renderWindow.Render()
        self.interactor.Start()
        
    def _generate_points(self):
        self.points = vtk.vtkPoints()
        for node in self.mesh.nodes:
            self.points.InsertNextPoint(node.location)
            
    def _build_quads(self):
        
        self.cellArray = vtk.vtkCellArray()
        
        for element in self.mesh.elements:
            quad = vtkQuad()
            quad.GetPointIds().SetId(0, element.nodes[0].id)
            quad.GetPointIds().SetId(1, element.nodes[1].id)
            quad.GetPointIds().SetId(2, element.nodes[2].id)
            quad.GetPointIds().SetId(3, element.nodes[3].id)
            
            self.cellArray.InsertNextCell(quad)
            
    def _build_viewing_heirarchy(self):
        self.poly = vtk.vtkPolyData()
        self.poly.SetPoints(self.points)
        self.poly.SetPolys(self.cellArray)
        
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.poly)
        
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetInterpolationToFlat()
        self.actor.GetProperty().SetEdgeColor(1.0, 0.0, 0.0) #(R,G,B)
        self.actor.GetProperty().EdgeVisibilityOn()
        
        self.renderer = vtk.vtkRenderer()
        self.renderer.AddActor(self.actor)
        
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.AddRenderer(self.renderer)
        
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.renderWindow)