#!/usr/bin/env python
"""Upload a set of next-gen sequencing data files to a data library in Galaxy.

Usage:
    upload_to_galaxy.py <config file> <flowcell directory> <analysis output dir>

The configuration file is in YAML format with the following key/value pairs:

galaxy_url: Base URL of Galaxy for uploading.
galaxy_api_key: Developer's API key.
galaxy_config: Path to Galaxy's universe_wsgi.ini file. This is required so
we know where to organize directories for upload based on library_import_dir.
"""
import sys
import os
import glob
import shutil
import ConfigParser
import urllib
import urllib2
import json

import yaml

from Mgh.Solexa.Flowcell import get_flowcell_info, get_fastq_dir

def main(config_file, fc_dir, analysis_dir):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    galaxy_api = GalaxyApiAccess(config["galaxy_url"],
        config["galaxy_api_key"])

    fc_name, fc_date = get_flowcell_info(fc_dir)
    folder_name = "%s_%s" % (fc_date, fc_name)
    run_info = lims_run_details(galaxy_api, fc_name)
    for dl_folder, access_role, dbkey, lane, name, description in run_info:
        print folder_name, lane, name, description, dl_folder
        library_id = get_galaxy_library(dl_folder, galaxy_api)
        folder, cur_galaxy_files = get_galaxy_folder(library_id, folder_name,
                                                     name, description,
                                                     galaxy_api)
        print "Creating storage directory"
        store_dir = move_to_storage(lane, folder_name,
                select_upload_files(lane, fc_dir, analysis_dir),
                cur_galaxy_files, config)
        if store_dir:
            print "Uploading directory of files to Galaxy"
            print galaxy_api.upload_directory(library_id, folder['id'],
                                              store_dir, dbkey, access_role)

# LIMS specific code for retrieving information on what to upload from
# the Galaxy NGLIMs.
# Also includes function for selecting files to upload from flow cell and
# analysis directories.
# These should be editing to match a local workflow if adjusting this.

def lims_run_details(galaxy_api, fc_name):
    """Retrieve run infomation on a flow cell from Next Gen LIMS.
    """
    run_info = galaxy_api.run_details(fc_name)
    for lane_info in (l for l in run_info["details"] if l.has_key("researcher")):
        description = "%s: %s" % (lane_info["researcher"],
                lane_info["description"])
        libname, role = _get_galaxy_libname(lane_info["private_libs"],
                                            lane_info["lab_association"])
        yield (libname, role, lane_info["genome_build"],
                lane_info["lane"], lane_info["name"], description)

def _get_galaxy_libname(private_libs, lab_association):
    # simple case -- one private library. Upload there
    if len(private_libs) == 1:
        return private_libs[0]
    # no private libraries -- use the lab association
    elif len(private_libs) == 0:
        return lab_association, ""
    # multiple libraries -- find the one that matches the lab association
    else:
        check_libs = [l.lower() for (l, _) in private_libs]
        try:
            i = check_libs.index(lab_association.lower())
            return private_libs[i]
        # can't find the lab association, give us the first library
        except IndexError:
            return private_libs[0]

def select_upload_files(lane, fc_dir, analysis_dir):
    """Select fastq, bam alignment and summary files for upload to Galaxy.
    """
    for fname in glob.glob(os.path.join(get_fastq_dir(fc_dir),
            "%s_*_fastq.txt" % lane)):
        yield fname
    for summary_file in glob.glob(os.path.join(analysis_dir,
            "%s_*-summary.pdf" % lane)):
        yield summary_file
    for bam_file in glob.glob(os.path.join(analysis_dir,
            "%s_*-sort-dup.bam" % lane)):
        yield bam_file

# General functionality for interacting with Galaxy via the Library API

def get_galaxy_folder(library_id, folder_name, lane, description, galaxy_api):
    """Return or create a folder within the given library.

    Creates or retrieves a top level directory for a run, and then creates
    a lane specific directory within this run.
    """
    items = galaxy_api.library_contents(library_id)
    root = _folders_by_name('/', items)[0]
    run_folders = _folders_by_name("/%s" % folder_name, items)
    if len(run_folders) == 0:
        run_folders = galaxy_api.create_folder(library_id,
                root['id'], folder_name)
    lane_folders = _folders_by_name("/%s/%s" % (folder_name, lane), items)
    if len(lane_folders) == 0:
        lane_folders = galaxy_api.create_folder(library_id,
                run_folders[0]['id'], str(lane), description)
        cur_files = []
    else:
        cur_files = [f for f in items if f['type'] == 'file'
                and f['name'].startswith("/%s/%s" % (folder_name, lane))]
    return lane_folders[0], cur_files

def _folders_by_name(name, items):
    return [f for f in items if f['type'] == 'folder' and
                                f['name'] == name]

def move_to_storage(lane, fc_dir, select_files, cur_galaxy_files, config):
    """Create directory for long term storage before linking to Galaxy.
    """
    galaxy_conf = ConfigParser.SafeConfigParser({'here' : ''})
    galaxy_conf.read(config["galaxy_config"])
    try:
        lib_import_dir = galaxy_conf.get("app:main", "library_import_dir")
    except ConfigParser.NoOptionError:
        raise ValueError("Galaxy config %s needs library_import_dir to be set."
                % config["galaxy_config"])
    storage_dir = _get_storage_dir(fc_dir, lane, os.path.join(lib_import_dir,
                                   "storage"))
    existing_files = [os.path.basename(f['name']) for f in cur_galaxy_files]
    upload_files = []
    for to_upload in select_files:
        if os.path.basename(to_upload) in existing_files:
            upload_files = []
            break
        else:
            if not os.path.exists(os.path.join(storage_dir,
                                               os.path.basename(to_upload))):
                shutil.copy(to_upload, storage_dir)
            upload_files.append(to_upload)
    return (storage_dir if len(upload_files) > 0 else None)

def _get_storage_dir(cur_folder, lane, storage_base):
    store_dir = os.path.join(storage_base, cur_folder, str(lane))
    if not os.path.exists(store_dir):
        os.makedirs(store_dir)
    return store_dir

def get_galaxy_library(lab_association, galaxy_api):
    ret_info = None
    for lib_info in galaxy_api.get_libraries():
        if lib_info["name"].find(lab_association) >= 0:
            ret_info = lib_info
            break
    # need to add a new library
    if ret_info is None:
        ret_info = galaxy_api.create_library(lab_association)[0]
    return ret_info["id"]

class GalaxyApiAccess:
    """Simple front end for accessing Galaxy's REST API.
    """
    def __init__(self, galaxy_url, api_key):
        self._base_url = galaxy_url
        self._key = api_key

    def _make_url(self, rel_url, params=None):
        if not params:
            params = dict()
        params['key'] = self._key
        vals = urllib.urlencode(params)
        return ("%s%s" % (self._base_url, rel_url), vals)

    def _get(self, url, params=None):
        url, params = self._make_url(url, params)
        response = urllib2.urlopen("%s?%s" % (url, params))
        return json.loads(response.read())

    def _post(self, url, data, params=None):
        url, params = self._make_url(url, params)
        request = urllib2.Request("%s?%s" % (url, params),
                headers = {'Content-Type' : 'application/json'},
                data = json.dumps(data))
        response = urllib2.urlopen(request)
        return json.loads(response.read())

    def get_libraries(self):
        return self._get("/api/libraries")

    def show_library(self, library_id):
        return self._get("/api/libraries/%s" % library_id)

    def create_library(self, name, descr="", synopsis=""):
        return self._post("/api/libraries", data = dict(name=name,
            description=descr, synopsis=synopsis))

    def library_contents(self, library_id):
        return self._get("/api/libraries/%s/contents" % library_id)

    def create_folder(self, library_id, parent_folder_id, name, descr=""):
        return self._post("/api/libraries/%s/contents" % library_id,
                data=dict(create_type="folder", folder_id=parent_folder_id,
                          name=name, description=descr))

    def show_folder(self, library_id, folder_id):
        return self._get("/api/libraries/%s/contents/%s" % (library_id,
            folder_id))

    def upload_directory(self, library_id, folder_id, directory, dbkey,
            access_role='', file_type='auto', link_data_only=True):
        """Upload a directory of files with a specific type to Galaxy.
        """
        return self._post("/api/libraries/%s/contents" % library_id,
                data=dict(create_type='file', upload_option='upload_directory',
                    folder_id=folder_id, server_dir=directory,
                    dbkey=dbkey, roles=str(access_role),
                    file_type=file_type, link_data_only=str(link_data_only)))

    def run_details(self, run):
        """Next Gen LIMS specific API functionality.
        """
        return self._get("/nglims/api_run_details", dict(run=run))

if __name__ == "__main__":
    main(*sys.argv[1:])
