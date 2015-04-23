#info
import h5py, os
from events_tools_pm import is_basecalled
import poreminion

def get_minknow_version(f5):
    return f5['/UniqueGlobalKey/tracking_id'].attrs["version_name"]

def get_basecaller_events_limits(f5):
    return f5['/Analyses/Basecall_2D_000/Configuration/general'].attrs['min_events'], f5['/Analyses/Basecall_2D_000/Configuration/general'].attrs['max_events']

def get_min_max_ratio(f5):
    return f5['/Analyses/Basecall_2D_000/Configuration/hairpin_align'].attrs['min_ratio'], f5['/Analyses/Basecall_2D_000/Configuration/hairpin_align'].attrs['max_ratio']

def get_model_type(f5):
    return f5['/Analyses/Basecall_2D_000/Configuration/general'].attrs['model_type']


def run_get_attributes(f5file):
    script = os.path.join(poreminion.scripts.__path__[0], 'getAttributes.sh')
    cmd = script + " " + f5file
    os.system(cmd)

## TODO: add more get functions
## TODO: add option allowing user to request specific attributes

def run(parser, args):
    if args.basic:
        f5 = h5py.File(args.fast5)
        basecalled = is_basecalled(f5)
        print get_minknow_version(f5)
        if basecalled:
            print get_basecaller_events_limits(f5)
            print get_min_max_ratio(f5)
            print get_model_type(f5)
    if args.all:
        run_get_attributes(args.fast5)
