#!/usr/bin/env python
"""Convert all datasets in a Galaxy library to a specific organism.

This coverts the dbkey for datasets in a library to a new organism, and is
useful for bulk updates of pre-loaded libraries.

Usage:
    covert_library_dbkey.py <config ini file> <library name> <organism dbkey>
"""
import sys
import os
import tempfile
import ConfigParser

def main(ini_file, library_name, dbkey):
    sys.path.append(os.path.join(os.getcwd(), "lib"))
    app = get_galaxy_app(ini_file)

    #for library in app.model.Library.query():
    #    print library.name, library.deleted

    library = app.model.Library.query().filter_by(name=library_name,
            deleted=False).first()
    app.model.session.begin()
    for dataset in library.root_folder.datasets:
        print 'Assigning', dataset.library_dataset_dataset_association.name, \
                'to', dbkey
        #print dataset.library_dataset_dataset_association.dbkey
        dataset.library_dataset_dataset_association.dbkey = dbkey
    app.model.session.commit()

def get_galaxy_app(ini_file):
    import galaxy.app

    conf_parser = ConfigParser.ConfigParser({'here':os.getcwd()})
    conf_parser.read(ini_file)
    configuration = {}
    for key, value in conf_parser.items("app:main"):
        configuration[key] = value
    # If we don't load the tools, the app will startup much faster
    empty_xml = tempfile.NamedTemporaryFile()
    empty_xml.write( "<root/>" )
    empty_xml.flush()
    configuration['tool_config_file'] = empty_xml.name
    configuration['enable_job_running'] = False
    configuration['database_create_tables'] = False
    app = galaxy.app.UniverseApplication( global_conf = ini_file, **configuration )
    return app

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print __doc__
        sys.exit()
    main(*sys.argv[1:4])
