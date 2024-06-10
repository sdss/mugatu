.. _mugatu-changelog:

==========
Change Log
==========

This document records the main changes to the ``mugatu`` code.

`2.3.1 <https://github.com/sdss/mugatu/compare/2.3.0...2.3.1>`_ - 2024-06-10

* Add option to validation script to skip dark_rm designs.
* Fix syntax issue in html script.

`2.3.0 <https://github.com/sdss/mugatu/compare/2.2.2...2.3.0>`_ - 2024-05-31

* Add validation to check designid_status coomapred to info in opsdb.

`2.2.2 <https://github.com/sdss/mugatu/compare/2.2.1...2.2.2>`_ - 2023-10-03

* rsValidation_html script now can generate summary plots for type = dir, rs, rs_catchup.
* validation script saves validation file to new subdirectory for type = dir.
* Validation script now writes to environmental variable called MUGATU_DATA.

`2.2.1 <https://github.com/sdss/mugatu/compare/2.2.0...2.2.1>`_ - 2023-08-25

* Below patches do not effect core functions of mugatu.
* Change to replacement of designs in targetdb with custom designs to recognize previously defined custom fields.
* Add designid_status column to results form validation script so designs created in a previous robostrategy run can be ignored when analyzing validation results.

`2.2.0 <https://github.com/sdss/mugatu/compare/2.1.0...2.2.0>`_ - 2023-07-26

* Add option in validation to used cached bright neighbor queries from robostrategy.
* Updated offsetting function calls for function versions that use all magnitudes and are on a per observatory basis.

`2.1.0 <https://github.com/sdss/mugatu/compare/2.0.0...2.1.0>`_ - 2023-04-24

* Changed designmode.build_brigh_neigh_query to default to v1 and use Gaia DR3. This means mugatu versions from here on need to use v1 targets and older versions should be used if using v0.5 targeting.

`2.0.0 <https://github.com/sdss/mugatu/compare/2.0.0-alpha.0...2.0.0>`_ - 2023-01-05

* Added stable versions for coordio and kaiju in dependencies.

`2.0.0-alpha.0 <https://github.com/sdss/mugatu/compare/1.3.2...2.0.0-alpha.0>`_ - 2022-11-29

* Added functions to allow for offsetting of targets in design based on offset function from coordio.
* Added optional parameter to mugatu.designmode.build_brigh_neigh_query to specify observatory and change search radius.
* Fixed bug in function that outputs next available field_id
* Script to create HTML page of RS validation results now does not require both observatories
* Validation script now validates fields in parallel
* Can pass observatory to bright neighbor query to change search radius

`1.3.2 <https://github.com/sdss/mugatu/compare/1.3.1...1.3.2>`_ - 2022-08-05

* Change bright safety query to exclude Tycho catalogids in Gaia

`1.3.1 <https://github.com/sdss/mugatu/compare/1.3.0...1.3.1>`_ - 2022-08-01

* Write design_version_pk for new design entries.

`1.3.0 <https://github.com/sdss/mugatu/compare/1.2.0...1.3.0>`_ - 2022-07-19

* Update mugatu.designmode.build_brigh_neigh_query to include proper motions, only query on the current version of catalogdb and remove Tycho/Gaia duplicates.
* Add optional argument to FPSDesign for RS_VERSION (needed to be used with new DesignToField schema).
* Make changes to queries to include new DesignToField table.
* Add field_exposure to designs added to targetdb.
* Make obsTime optional argument in mugatu.fspdesign.FPSDesign so it is calculated at LST when not provided.

`1.2.0 <https://github.com/sdss/mugatu/compare/1.1.3...1.2.0>`_ - 2022-03-23

* Update code to accommodate new length of magnitude array to N = 10. Now all magnitude inputs need to be of form: [g, r, i, z, BP, G, RP, J, H, K]. Older versions of mugatu will not be compatible with this new format, so some functions will be broken if trying to use older versions.

`1.1.3 <https://github.com/sdss/mugatu/compare/1.1.2...1.1.3>`_ - 2022-03-09

* Add function to create design_status bit mask.
* Allow for replacement fields when creating new fields in targetdb.

`1.1.2 <https://github.com/sdss/mugatu/compare/1.1.1...1.1.2>`_ - 2022-03-08

* Change designmode category accounting to use category string and not carton_pk.

`1.1.1 <https://github.com/sdss/mugatu/compare/1.1.0...1.1.1>`_ - 2022-02-15

* For bright neighbor check, only fail check when assigned fibers are too near a bright source (i.e. no longer include unassigned fibers).
* In DesignMode outputs for bright neighbor check, include adjusted magnitudes for fibers near bright sources.

`1.1.0 <https://github.com/sdss/mugatu/compare/1.0.4...1.1.0>`_ - 2022-02-03

* Update coordio/kaiju syntax to be compatible with coordio=1.2.1 and kaiju=1.2.2.
* Add function to calculate assignment_hash for designs. Used to identify identical designs.
* Check for reserved fieldids  in TargetFieldIDs object.
* Add code to create summary of validation results as HTML page.
* Add functionality to write validation results to new DesignModeCheckResults targetdb table.
* Remove restriction to only consider exclusion radii > 1" for bright neighbor check (so, consider all exclusion radii for bright sources).

`1.0.4 <https://github.com/sdss/mugatu/compare/1.0.3...1.0.4>`_ - 2021-12-08

* In TargetdbFieldIDs class, account for gaps in fieldid when finding next available.
* Add values for designmode metrics to FPSDesign object design_errors dictionary.
* Specified additional outer joins when querying design in targetdb.
* Specified mugatu version and run_on date in design table when ingesting new designs.

`1.0.3 <https://github.com/sdss/mugatu/compare/1.0.2...1.0.3>`_ - 2021-11-29

* Added some minor tweaks for database joins so designs are pulled from database completely.

`1.0.2 <https://github.com/sdss/mugatu/compare/1.0.1...1.0.2>`_ - 2021-11-29

* Added some minor tweaks to the database column names used to pull designs from targetdb.

`1.0.1 <https://github.com/sdss/mugatu/compare/1.0.0...1.0.1>`_ - 2021-11-29

* Changes have been made to make mugatu compatible with the new targetdb schema as of sdssdb=0.4.12.
* The bright neighbor check has been added to the verification of designs.
* Finalized versions of all designmode checks are included in the verification of designs.
* A new class has been added to check availability of field_id in targetdb.Field table.
