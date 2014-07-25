GeoMACH
-------
GeoMACH (Geometry-centric MDAO of Aircraft Configurations with High fidelity) is an aircraft design tool suite featuring parametric modelers for the geometry and structure of the aircraft that are efficient and supply derivatives with respect to all parameters, enabling optimization with large numbers of design variables.

Installing and getting started (linux or Mac)
---------------------------------------------
[GeoMACH-top] = path to the top-level GeoMACH directory

1. Install GeoMACH:
   $ cd [GeoMACH-top]
   $ sudo python setup.py install

2. Compile the documentation:
   $ cd [GeoMACH-top]/docs
   $ make html

3. Read the documentation:
   (option 1, if firefox is installed)
   $ cd [GeoMACH-top]/docs/_build/html
   $ firefox index.html

   (option 2)
   Open a browser and enter the URL:
   [GeoMACH-top]/docs/_build/html/index.html

4. Run the example:
   $ cd [GeoMACH-top]/examples
   $ python conventional.py
