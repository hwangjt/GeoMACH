GeoMACH Package
===============

Abstract
--------

This document presents a high-level description of
*GeoMACH: Geometry-centric MDO of Aircraft Configurations with High-fidelity*.
Anyone interested in GeoMACH is strongly encouraged to read this page and the
BSE and PGM pages in their entirety, as they contain a concise overview of the
objectives, features, approach, and implementation of GeoMACH.
The intended audience includes:

1. People considering using GeoMACH
2. Beginner users of GeoMACH
3. Advanced users / developers

Motivation
----------

GeoMACH's main purpose is parametric geometry modeling---it is a tool
for creating and manipulating aircraft geometries.
There are many existing software packages for this, but GeoMACH is unique in
simultaneously achieving the following objectives:

+---+-----------------+-----------------------------------------------------+
|   | Objective       | GeoMACH feature                                     |
+===+=================+=====================================================+
| 1 | High-fidelity   | - Fully differentiable parametrization              |
|   | geometries for  | - Sparse, analytic derivatives                      |
|   | large-scale     | - Efficient: :math:`\mathcal{O}` (1s) to initialize |
|   | optimization    |   and :math:`\mathcal{O}` (100ms) to update         |
|   |                 | - Accurate: continuous, watertight geometry         |
+---+-----------------+-----------------------------------------------------+
| 2 | Conceptual      | - Support for unconventional configurations         |
|   | design-level    |                                                     |
|   | parametrization |                                                     |
+---+-----------------+-----------------------------------------------------+
| 3 | Ease-of-use     | - Open-source software                              |
|   | and automation  | - Parametrization using high-level aircraft design  |
|   |                 |   parameters, or low-level shape variables, or      |
|   |                 |   anything in between                               |
+---+-----------------+-----------------------------------------------------+

The commercial computer-aided design (CAD) packages that are popular today use
constructive solid geometry (CSG), which defines geometries as a sequence of 
boolean operations on simple 3-D objects like cylinders.
There are many geometry engines specialized for aircraft design that compute
intersections between aircraft components to build a single, closed geometry.
Both CAD packages and aircraft geometry engines do not satisfy the first
objective because they are not differentiable, do not provide accurate
derivatives, and/or are inefficient to evaluate.

Approach
--------

GeoMACH achieves the above objectives by representing the aircraft
outer mold line (OML) in a unique way:

1. *The OML is a union of 4-sided patches*: The 4-sided patches facilitate
   smoothly interpolated junctions between components. These junctions
   are the key to enabling a differentiable parametrization.

   .. image:: Images/junction.jpg
      :width: 200 px

2. *Each patch is a B-spline surface*: In a B-spline surface, a 2-D 
   array of control points defines a smooth and continuous 4-sided patch.
   One advantage of using B-splines is that a surface mesh is defined 
   by a sparse linear mapping from the control points.

   .. image:: Images/bsplines.png
      :width: 200 px

3. *Decomposition of aircraft into components*: The control points are
   in turn defined by shape variables which are grouped by aircraft
   components.

   .. image:: Images/primitives.png
      :width: 200 px

Subpackages
-----------

GeoMACH is a Python package that uses Fortran subroutines for expensive
computations. It is divided into 3 subpackages:

1. BSE: B-spline Surface-modeling Engine - a general, standalone tool
   for modeling geometries as a union of B-spline surfaces
2. PGM: Parametric Geometry Modeler - an aircraft geometry modeling tool
   that maps user-defined shape variables to B-spline control points to
   be passed to BSE
3. PSM: Parametric Structural Modeler - aircraft structural modeling tool
   (in development)

.. toctree::

    GeoMACH.BSE
    GeoMACH.PGM
    GeoMACH.PSM

