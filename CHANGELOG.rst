.. _mugatu-changelog:

==========
Change Log
==========

This document records the main changes to the ``mugatu`` code.

* :release:`1.3.1 <2022-08-01>`
* Write design_version_pk for new design entries.
* :release:`1.3.0 <2022-07-19>`
* Update mugatu.designmode.build_brigh_neigh_query to include proper motions, only query on the current version of catalogdb and remove Tycho/Gaia duplicates.
* Add optional argument to FPSDesign for RS_VERSION (needed to be used with new DesignToField schema).
* Make changes to queries to include new DesignToField table.
* Add field_exposure to designs added to targetdb.
* Make obsTime optional argument in mugatu.fspdesign.FPSDesign so it is calculated at LST when not provided.

* :release:`1.2.0 <2022-03-23>`
* Update code to accommodate new length of magnitude array to N = 10. Now all magnitude inputs need to be of form: [g, r, i, z, BP, G, RP, J, H, K]. Older versions of mugatu will not be compatible with this new format, so some functions will be broken if trying to use older versions.

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
