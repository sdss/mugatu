.. _db-design-example:

FPS deisgns from targetdb
=========================

In order use ``mugatu`` to validate a design that exists in ``targetdb``, all that is needed is the primary key of the design in the ``targetdb.Design`` table and the Julian Date of the observation of this design. Additionally, you have to be connected to targetdb using the ``sdssdb`` package. Once this is done, a ``mugatu.FPSDesign`` object can be built:

	>>> from mugatu.fpsdesign import FPSDesign
	>>> des = FPSDesign(design_pk=3, obsTime=2459519.51)

Then, the design can be validated using ``kaiju``, where validation in this context means the design is checked for possible collisions and deadlocks of the robotic arms of the FPS. The design is built and validated by:

	>>> des.build_design_db()
	>>> des.validate_design()
	[WARNING]: Some targets could not be assigned to fiber (MugatuWarning)
	[WARNING]: Some targets removed from design due to collisions (MugatuWarning)

As can be seen, ``mugatu`` provides warnings for when targets either could not be assigned to the fiber because the assigned fiber could not reach that position or that the assinged fiber would cause a collision. To check which objects were removed based on these criteria you can do the following:

	>>> des.targets_unassigned
	[5162927115,
 	5162842743,
 	5160549269,
 	5160550675,
 	5162920152,
 	5162963877,
 	5139069961,
 	5163132766,
 	5163169243,
 	5163095154,
 	5163094799,
 	5163000011,
 	5139076236,
 	5139069478]
 	>>> des.targets_collided
 	[5162889296,
 	5162960038,
 	5163076735,
 	5163068285,
 	5162935471,
 	5162916769,
 	5162845224,
 	5162880987,
 	5139061379,
 	5139065658,
 	5163167280,
 	5163162899,
 	5163158588,
 	5163062228,
 	5162996112]

 Where the lists returned above are the catalogids of the targets that could not be assigned to the chosen fiber.

 With these removals, the final validated design is then stored in ``des.valid_design``.

.. _manual-design-example:

Build a design manually
========================