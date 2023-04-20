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

.. _validation-design-example:

Validate Designs with mugatu
============================

With ``mugatu`` you are able to validate designs in batches using the script ``bin/validate_designs_batches.py``. For validating designs, you can either validate designs in a user specified directory, all of a robostrategy run or a set of designs meant to replace some robostrategy designs.

User Specified Directory
------------------------

For designs in a user specified directory, one can validate designs like the following:

	>>> python bin/validate_designs_batches.py -t dir -l utah -d path_to_designs/

This will then output validation results in a file called:

::

	path_to_designs/design_validation_results.fits


Robostrategy Plan
-----------------

To validate an entire robostrategy run for one observatory, one would run the following:

	>>> python bin/validate_designs_batches.py -t rs -l utah -p plan_name -o apo -n 16

In the above call, we have also added ``-n 16`` to run the validaiton in parallel with 16 cores. This will output validation results in a file called:

::

	/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/sandbox/mugatu/rs_plan_validations/{plan_name}/rs_{plan_name}_apo_design_validation_results.fits


Robostrategy Replacement Designs
--------------------------------

Finally, you can validate designs that will replace an entire robostrategy field that has already been ingested into targetdb. This can be done via to following:

	>>> python bin/validate_designs_batches.py -t rs_replace -l utah -p plan_name -f 1000 1001 - n 16

In the above, ``-f 1000 10001`` indicated that we will validate the deisgns meant to replace the designs in fields 1000 and 1001. This will output validation results in a file called:

::

	/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/target/robostrategy_replacement/{plan_name}/{original_file_name}_validation.fits


.. _ingest-design-example:

Ingest Designs with mugatu
==========================

With ``mugatu`` you are able to ingest validated designs in targetdb. Similar to the above example, this can be done for designs in a user specified directory, for a robostrategy run or for a set of designs meant to replace some robostrategy designs.

User Specified Directory
------------------------

For designs in a user specified directory, one can ingest designs like the following:

	>>> python bin/load_manual_designs_batches.py -l utah -d path_to_designs/ -f path_to_designs/design_validation_results.fits.

In addition, this code will out put a file that associates the deisgn files with their design_id in targetdb:

::

	path_to_designs/design_ids_for_design_files.fits


Robostrategy Plan
-----------------

To load designs for an entire robostrategy, one can ingest all designs like the following:

	>>> python bin/RS_to_targetdb.py -p plan_name -o apo -t example_load

In the above, ``-t example_load`` tags the new version for the plan in targetdb.

Robostrategy Replacement Designs
--------------------------------

Finally, we can ingest designs that are meant to replace all designs in a field for a current robostrategy plan. This can be done like the following:

	>>> python bin/replace_RS_designs.py -l utah -p plan_name -f 1000 1001

The above will also create/update a change log for the new field_ids that may be created as a result of the design replacement process.


.. _full-example:

Full Example at Utah
====================

Here is a full worked example on how to validate, visualize the validation results and ingest the designs of a robostrategy run at Utah. First, the required software needs to be loaded via Modules:

	>>> module load miniconda/3.9 robostrategy tree

The above should load all the default versions of the required software. Next, we will validate the designs for each observatory. This is best done in parallel by submiting a slurm job script. Below is an example of a slurm job script called ``submit_validate.sh`` that will accomplish this:

::

	#!/bin/bash
	#SBATCH --job-name example_valid_apo
	#SBATCH --output=example_valid_apo.txt
	#SBATCH --account=sdss-np
	#SBATCH --partition=sdss-np
	#SBATCH --nodes=1
	#SBATCH --time=3-00:00:00
	srun -n 1 -c 8 python -u $MUGATU_DIR/bin/validate_designs_batches.py -t rs -l utah -p example_plan -o apo -n 8

The above will run the validation at APO for a robostrategy run called ``example_plan`` using 8 cores. The job is then submited by:

	>>> sbatch submit_validate.sh

The above then needs to be repeated for LCO by changing ``-o apo`` to ``-o lco`` in ``submit_validate.sh``.

Once completed, the results of the validations can be found here:

::

	/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/sandbox/mugatu/rs_plan_validations/{plan_name}/rs_{plan_name}_{obs}_design_validation_results.fits

Where in this example ``{plan_name} = example_plan`` and there will be two files, one with ``{obs} = apo`` and the other with ``{obs} = lco``. To better visualize the results, an HTML page can be generated with summary plots. This is done by running the command:

	>>> python $MUGATU_DIR/bin/rsValidation_html.py -p example_plan --kaiju_v 1.3.1 --coordio_v 1.5.2

In the above, the inputs for ``--kaiju_v`` and ``--coordio_v`` should be changed to the current versions of kaiju and coordio. These can be seen by using the command ``module list``. Once the above script is finished, the validation results can then be viewed at:

::

	https://data.sdss5.org/sas/sdsswork/sandbox/mugatu/rs_plan_validations/{plan_name}/{plan_name}_validation.html

On this page will be links to the various tests so one can assess the validation results and look for any major issues with the robostrategy run.

If the above validation is verified as acceptable, then the designs can be loaded into targetdb. This process can take some time, so it is best to run it in the background. This can be done via the following command:

	>>> nohup python -u  $MUGATU_DIR/bin/RS_to_targetdb.py -p example_plan -o apo -t example_plan > log_ingest_apo.txt 2>&1 &

In the above ``-p`` sets the plan name of a new entry ``targetdb.Version`` and ``-t`` sets the tag for that plan. In practice it is best if these are the same. Finally, the above command also needs to be run for LCO by changing ``-o apo`` to ``-o lco`` in the above. With the above completed, all designs should now be ingested into targetdb at Utah.
