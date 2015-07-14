GeoMACH
-------
GeoMACH (Geometry-centric MDAO of Aircraft Configurations with High fidelity) is an aircraft design tool suite featuring parametric modelers for the geometry and structure of the aircraft that are efficient and supply derivatives with respect to all parameters, enabling optimization with large numbers of design variables.

Installing and getting started (linux or Mac)
---------------------------------------------
[GeoMACH-top] = path to the top-level GeoMACH directory

1. Install numpy, scipy, sphinx, and numpydoc

2. Install GeoMACH:
   $ cd [GeoMACH-top]
   $ sudo python setup.py develop

3. Compile the documentation:
   $ cd [GeoMACH-top]/docs
   $ make html

4. Read the documentation:
   (option 1, if firefox is installed)
   $ cd [GeoMACH-top]/docs/_build/html
   $ firefox index.html

   (option 2)
   Open a browser and enter the URL:
   [GeoMACH-top]/docs/_build/html/index.html

5. Run the example:
   $ cd [GeoMACH-top]/examples
   $ python conventional.py

Installing on a local directory (linux or Mac, for clusters)
---------------------------------------------
[GeoMACH-top] = path to the top-level GeoMACH directory

1. Install (or load) numpy, scipy

2. Install GeoMACH:
   $ cd [GeoMACH-top]
   $ python setup.py build --fcompiler=intelem (you may use another compiler)
   $ python setup.py install --user

3. You have to manually copy all airfoils to the installation folder
   $ cp -r GeoMACH/PGM/airfoils ~/.local/lib/python2.7/site-packages/GeoMACH-0.1-py2.7-linux-x86_64.egg/GeoMACH/PGM/

4. Run the example:
   $ cd examples
   $ python conventional.py
