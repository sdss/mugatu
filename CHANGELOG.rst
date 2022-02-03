.. _mugatu-changelog:

==========
Change Log
==========

This document records the main changes to the ``mugatu`` code.

* Update coordio/kaiju syntax to be compatible with coordio=1.2.1 and kaiju=1.2.2.
* Add function to calculate assignment_hash for designs. Used to identify identical designs.
* Check for reserved fieldids  in TargetFieldIDs object.
* Add code to create summary of validation results as HTML page.
* Add functionality to write validation results to new DesignModeCheckResults targetdb table.
* Remove restriction to only consider exclusion radii > 1" for bright neighbor check (so, consider all exclusion radii for bright sources).

* :release:`1.0.4 <2021-12-08>`
* In TargetdbFieldIDs class, account for gaps in fieldid when finding next available.
* Add values for designmode metrics to FPSDesign object design_errors dictionary.
* Specified additional outer joins when querying design in targetdb.
* Specified mugatu version and run_on date in design table when ingesting new designs.

* :release:`1.0.3 <2021-11-29>`
* Added some minor tweaks for database joins so designs are pulled from database completely.

* :release:`1.0.2 <2021-11-29>`
* Added some minor tweaks to the database column names used to pull designs from targetdb.

* :release:`1.0.1 <2021-11-29>`
* Changes have been made to make mugatu compatible with the new targetdb schema as of sdssdb=0.4.12.
* The bright neighbor check has been added to the verification of designs.
* Finalized versions of all designmode checks are included in the verification of designs.
* A new class has been added to check availability of field_id in targetdb.Field table.
