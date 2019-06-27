
==============
Code Structure
==============

Class-Inheritance
-----------------

The main class in PROPOSAL is the Propagator class.

The code structure is described in detail in the `second PROPOSAL paper <https://arxiv.org/abs/1809.07740>`_.

The schematic code structure is displayed in the image below.

.. figure:: inheritance.png
   :scale: 50 %
   :align: center

The Propagator Class consists of a list of sectors and decides which sector propagates the particel.
Each sector then propagates the particle as shown above.

Propagation Routine
-------------------

The propagation routine is in the Propagator class and is the ``Propagator::propagate(particle)`` function. The steps which PROPOSAL does during the propagation are best explained in the figure below.

.. figure:: PROPOSAL_propagation_flow.png
   :scale: 30 %
   :align: center


Config files
----------------
The media and geometry used for the propagation are defined in a JSON file. In ``resources/config_docu.md`` the entries are thoroughly explained.
