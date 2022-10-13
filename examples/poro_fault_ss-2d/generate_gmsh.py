#!/usr/bin/env nemesis
"""Generate a tri or quad mesh of a strike-slip fault using Gmsh, making use of the
built-in geometry engine.

We use the `gmsh_utils` module provided with PyLith. This module has helper functions
and classes for marking materials and boundaries compatible with PyLith. We also use the
`GenerateMesh` class from the module as a base class to our local `App` class. The
`GenerateMesh` class handles processing of command line options, initialization, and
finalizing the mesh.

Run `generate_gmsh.py --help` to see the command line options.
"""

import gmsh
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Block is 150km by 100km.

    -70.0 km <= x <= 75.0 km
    -50.0 km <= y <= 50.0 km

    The fault surface runs along the x-axis through the entire domain.

             p1--------------------p4
             |                      |
             |                      |
 fault ----> p5--------------------p6
             |                      |
             |                      |
    y=0 -->  p2--------------------p3
             ^
             |
             x=0
    """
    DOMAIN_X = 150.0e+3
    DOMAIN_Y = 100.0e+3

    def __init__(self):
        """Constructor.
        """
        # Set the cell choices available through command line options
        # with the default cell type `tri` matching the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
            }
        self.filename = "mesh_tri.msh"

    def create_geometry(self):
        """Create geometry.
        """
        # Set local variables for domain size and corner of the domain.
        lx = self.DOMAIN_X
        ly = self.DOMAIN_Y
        x1 = -0.5 * lx
        y1 = -0.5 * ly

        # Create points.
        p1 = gmsh.model.geo.add_point(x1, y1, 0.0)
        p2 = gmsh.model.geo.add_point(x1+lx, y1, 0.0)
        p3 = gmsh.model.geo.add_point(x1+lx, y1+ly, 0.0)
        p4 = gmsh.model.geo.add_point(x1, y1+ly, 0.0)

        p5 = gmsh.model.geo.add_point(x1+0.5*lx, y1, 0.0)
        p6 = gmsh.model.geo.add_point(x1+0.5*lx, y1+ly, 0.0)

        # Create curves. We store the curve tag as a data member
        # so that we can refer to them later.
        self.l_yneg1 = gmsh.model.geo.add_line(p1, p5)
        self.l_yneg2 = gmsh.model.geo.add_line(p5, p2)
        self.l_xpos = gmsh.model.geo.add_line(p2, p3)
        self.l_ypos2 = gmsh.model.geo.add_line(p3, p6)
        self.l_ypos1 = gmsh.model.geo.add_line(p6, p4)
        self.l_xneg = gmsh.model.geo.add_line(p4, p1)
        self.l_fault = gmsh.model.geo.add_line(p5, p6)

        # Create curve loops and surfaces from the curves.
        # We traverse the curves in a counter clock-wise direction.
        # If the curve is in the opporite direction, we use the negative tag.
        c0 = gmsh.model.geo.add_curve_loop([self.l_yneg1, self.l_fault, self.l_ypos1, self.l_xneg])
        self.s_xneg = gmsh.model.geo.add_plane_surface([c0])
        c1 = gmsh.model.geo.add_curve_loop([self.l_yneg2, self.l_xpos, self.l_ypos2, -self.l_fault])
        self.s_xpos = gmsh.model.geo.add_plane_surface([c1])

        gmsh.model.geo.synchronize()

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented.
        """
        # Create two materials, one for each side of the fault.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_xneg]),
            MaterialGroup(tag=2, entities=[self.s_xpos]),
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the fault.
        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.l_xneg]),
            VertexGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.l_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.l_yneg1, self.l_yneg2]),
            VertexGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.l_ypos1, self.l_ypos2]),
            VertexGroup(name="fault", tag=20, dim=1, entities=[self.l_fault]),
        )
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh.
        """
        # Set discretization size with geometric progression from distance to the fault.
        # We turn off the default sizing methods.
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # First, we setup a field `field_distance` with the distance from the fault.
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "CurvesList", [self.l_fault])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression` 
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=4.0e+3, bias=1.05)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        # Finally, we use the field `field_size` for the cell size of the mesh.
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        gmsh.option.setNumber("Mesh.Algorithm", 8)
        gmsh.model.mesh.generate(2)
        if cell == "quad":
            gmsh.model.mesh.recombine()
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()


# End of file
