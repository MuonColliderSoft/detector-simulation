#!/usr/bin/env python
"""This script converts a MARS file to an SLCIO file with MCParticles that can be used by ILCSoft"""

import os
import argparse


parser = argparse.ArgumentParser(description='Convert MARS text file to SLCIO file with MCParticles')
parser.add_argument('file_in', metavar='FILE_IN', help='Path to the input MARS15 file')
parser.add_argument('file_out', metavar='FILE_OUT.slcio', help='Path to the output SLCIO file')
parser.add_argument('-d', '--decays', metavar='D',  help='Number of decays simulated in the input sample (default: 5e5)', type=float, default=5e5)
parser.add_argument('-f', '--format',  help='Input format: [mars15 fluka]', type=str, default='mars15')
parser.add_argument('-i', '--invert_pdgs', help='Invert PDG IDs selection', default=False, action='store_true')
parser.add_argument('-l', '--lines_event', metavar='L',  help='Number of lines to put in a single LCIO event (default: 1000)', type=int, default=1000)
parser.add_argument('-m', '--max_lines', metavar='M',  help='Maximum number of lines to process', type=int, default=None)
parser.add_argument('-n', '--ne_min', metavar='E',  help='Minimum energy of accepted neutrons [GeV]', type=float, default=None)
parser.add_argument('-o', '--overwrite',  help='Overwrite existing output file', action='store_true', default=False)
parser.add_argument('-p', '--pdgs', metavar='ID',  help='PDG IDs of particles to be included', type=int, default=None, nargs='+')
parser.add_argument('-t', '--t_max', metavar='T',  help='Maximum time of accepted particles [ns]', type=float, default=None)
parser.add_argument('-w', '--weight', metavar='W',  help='Weight nominator in the input sample (default: 427780)', type=float, default=427780)

args = parser.parse_args()

if not args.overwrite and os.path.isfile(args.file_out):
    raise FileExistsError('Output file already exists: {0:s}'.format(args.file_out))

print("Converting {0:s} data from file\n  {1:s}\nto SLCIO file\n  {2:s}\nwith N decays: {3:.0f}".format(args.format, args.file_in, args.file_out, args.decays))
print("Splitting by {0:d} lines/event".format(args.lines_event));
if args.pdgs is not None:
    flag = 'in' if not args.invert_pdgs else 'NOT in'
    print("Will only use particles with PDG IDs {0:s} the list: {1}".format(flag, args.pdgs))


from math import sqrt
from pdb import set_trace as br
from array import array
from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL

import random
import math

from bib_pdgs import MARS_PIDS, FLUKA_PIDS, PDG_PROPS

def line_to_numbers(iL, line, line_format):
    """Converts the text line to a list of numbers according to the input format"""
    # Splitting the values into words and converting values to numbers
    ss = line.split()
    if line_format == 'm':
        n_words = 18
    elif line_format == 'f':
        # n_words = 17
        n_words = 20
    else:
        raise RuntimeError('Wrong line format: ' + line_format)

    if len(ss) != n_words:
        print('WARNING: wrong number of values on line {0:d}: {1:d}'.format(iL, len(ss)))
        return None

    # Extracting the PDGID
    if line_format == 'f':
        pdg = int(ss[0])
    elif line_format == 'm':
        pid = int(ss[1])
        if pid not in MARS_PIDS:
            print('WARNING! No PDG available for particle with MARS ID: {0:d}'.format(pid))
            print('         Skipping the particle...')
            return None
        pdg = MARS_PIDS[pid]
    # Skipping if this pdgId should be excluded
    if args.pdgs is not None:
        if args.invert_pdgs and pdg in args.pdgs:
            return None
        elif not args.invert_pdgs and pdg not in args.pdgs:
            return None

    if line_format == 'm':
        x,y,z,px,py,pz,toff,w = map(float, ss[2:10])
        toff *= 1e9  # s -> ns
    elif line_format == 'f':
        x,y,z,px,py,pz,toff,w = map(float, [ss[i] for i in [6,7,8, 2,3,4, 9, 5]])

    # Converting cm -> mm
    x *= 10
    y *= 10
    z *= 10

    return pdg, x,y,z, px,py,pz, toff, w




wrt = IOIMPL.LCFactory.getInstance().createLCWriter()
wrt.open( args.file_out, EVENT.LCIO.WRITE_NEW )

# Write a RunHeader
run = IMPL.LCRunHeaderImpl()
run.setRunNumber(0)
run.parameters().setValue("InputPath", args.file_in)
run.parameters().setValue("InputFormat", args.format)
if args.format == 'mars15':
    run.parameters().setValue("NDecays", args.decays)
    run.parameters().setValue("WNominator", args.weight)
if args.t_max:
    run.parameters().setValue("Time_max", args.t_max)
if args.ne_min:
    run.parameters().setValue("NE_min", args.ne_min)
if args.pdgs:
    run.parameters().setValue("PdgIds", str(args.pdgs))
    run.parameters().setValue("PdgIdsInverted", args.invert_pdgs)
wrt.writeRunHeader(run)

# Calculating the correction to particle weights (from the original MARS2ROOT conversion script)
w_corr = 1.0
if args.format == 'mars15':
    w_corr = args.weight/args.decays

random.seed()
nEventLines = 0
nEvents = 0
col = None
evt = None

line_format = args.format[0]

with open(args.file_in) as file_in:
    for iL, line in enumerate(file_in):
        if nEventLines == 0:
            col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE)
            # Creating the LCIO event and collection
            evt = IMPL.LCEventImpl()

            evt.setEventNumber(nEvents)
            evt.addCollection(col, "MCParticle")
        if args.max_lines and iL >= args.max_lines:
            break

        # Extracting relevant values from the line
        pdg, x,y,z, px,py,pz, toff, w = line_to_numbers(iL, line, line_format)

        # Skipping if particle's time is greater than allowed
        if args.t_max is not None and toff > args.t_max:
            continue

        # Skipping if it's a neutron with too low kinetic energy
        if args.ne_min is not None and abs(pdg) == 2112 and sqrt(px*px + py*py + pz*pz) < args.ne_min:
            continue

        # Correcting the weight
        w *= w_corr

        # Getting the charge and mass of the particle
        if pdg not in PDG_PROPS:
            print('WARNING! No properties defined for PDG ID: {0:d}'.format(pdg))
            print('         Skpping the particle...')
            continue
        charge, mass = PDG_PROPS[pdg]

        # print(iL, pdg, toff, w)

        # Calculating how many copies of the particle to create according to the weight
        nW = math.floor(w)
        if random.random() < math.modf(w)[0]:
            nW += 1

        # Creating the particle with original parameters
        particle = IMPL.MCParticleImpl()
        particle.setPDG(pdg)
        particle.setGeneratorStatus(1)
        particle.setTime(toff)
        particle.setMass(mass)
        particle.setCharge(charge)
        pos = array('d', [x, y, z])
        mom = array('f', [px, py, pz])

        # Creating the particle copies with random Phi rotation
        for iP in range(nW):
            p = IMPL.MCParticleImpl(particle)
            # Rotating position and momentum by random angle in Phi
            dPhi = random.random() * math.pi * 2
            co = math.cos(dPhi)
            si = math.sin(dPhi)
            pos[0] = co * x - si * y
            pos[1] = si * x + co * y
            mom[0] = co * px - si * py
            mom[1] = si * px + co * py
            p.setVertex(pos)
            p.setMomentum(mom)
            # Adding particle to the collection
            col.addElement(p)

        # Updating counters
        nEventLines += 1
        if nEventLines >= args.lines_event:
            nEvents +=1
            nEventLines = 0
            wrt.writeEvent(evt)
            print('Wrote event: {0:d}  (line {1:d})'.format(nEvents, iL+1))

    # Writing the last event
    if nEventLines > 0:
        wrt.writeEvent(evt)
        nEvents += 1
        print('Wrote event: {0:d}'.format(nEvents))

print('Wrote {0:d} events to file: {1:s}'.format(nEvents, args.file_out))
wrt.close()
