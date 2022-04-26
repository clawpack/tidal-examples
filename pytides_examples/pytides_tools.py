"""
Functions to facilitate using PyTides.
"""

import sys, os

# Temporary fix to import local submodule version of pytides:
CLAW = os.environ['CLAW']  # path to Clawpack files
pathstr = os.path.join(CLAW, 'tidal-examples/pytides')
assert os.path.isdir(pathstr), '*** clawpack/tidal-examples/pytides ***'
print('pytides_tools is using Pytides from: %s' % pathstr)
if pathstr not in sys.path:
    sys.path.insert(0,pathstr)

from pytides.tide import Tide
import numpy as np

def new_tide_instance_from_existing(constit_list,existing_tide_instance):
    """
    constit_list is the list of constituents to be used in the
    new_tide_instance.
    The values of the amplitudes and phases for each of them is to be
    pulled from an existing_tide_instance.  If no such constituent is in
    the existing_tide_instance, an error message is printed.
    """
    existing_constits = existing_tide_instance.model['constituent']
    existing_amps = existing_tide_instance.model['amplitude']
    existing_phases = existing_tide_instance.model['phase']
    len_existing = len(existing_constits)
    new_model = np.zeros(len(constit_list), dtype = Tide.dtype)
    new_constits=[]; new_amps=[]; new_phases=[];
    for ic in constit_list:
        success = False
        for j in range(len_existing):
            ie = existing_constits[j]
            if (ie.name == ic.name):    #grab it
                success = True
                new_constits.append(ie)
                new_amps.append(existing_amps[j])
                new_phases.append(existing_phases[j])
        if (success == False):
            print ('Did not find consituent name: ',ic.name,\
                   'in existing tide instance')
    new_model['constituent'] = new_constits
    new_model['amplitude']   = new_amps
    new_model['phase']       = new_phases

    # The new_model is now complete, so make a tide instance called
    #called new_tide from it.
    new_tide = Tide(model = new_model, radians = False)
    return new_tide
