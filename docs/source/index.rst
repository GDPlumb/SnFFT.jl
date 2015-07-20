SnFFT Package
===============

The *SnFFT* package is Julia package designed to facilitate harmonic analysis on the `symmetric group <http://en.wikipedia.org/wiki/Symmetric_group>`_ of order n, denoted |Sn|.   
Out of the box, *SnFFT* implements:

* Group operations and factorizations for |Sn|
* Functionality to set up functions over |Sn|
* The fast Fourier transform with additional options if the function is sparse or bandlimited
* The inverse fast Fourier transform with additional options if the function is bandlimited or the user is only interested in the result from the top few components
* The convolution and correlation of two Fourier transforms

**Contents:**

.. toctree::
   :maxdepth: 2

   starting.rst
   examples.rst
   groupoperations.rst
   snfunctions.rst
   yor.rst
   fft.rst
   ifft.rst
   corcon.rst
   miscellaneous.rst
   mlapplications.rst
   supplement.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |Sn| replace:: **S**\ :sub:`n` \
