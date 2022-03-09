.. _mugatu-changelog:

==========
Change Log
==========

This document records the main changes to the ``mugatu`` code.

* :release:`1.1.3 <2022-03-09>`
* Add function to create design_status bit mask.
* Allow for replacement fields when creating new fields in targetdb.

* :release:`1.1.2 <2022-03-08>`
* Change designmode category accounting to use category string and not carton_pk.

* :release:`1.1.1 <2022-02-15>`
* For bright neighbor check, only fail check when assigned fibers are too near a bright source (i.e. no longer include unassigned fibers).
* In DesignMode outputs for bright neighbor check, include adjusted magnitudes for fibers near bright sources.

* :release:`1.1.0 <2022-02-03>`
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
