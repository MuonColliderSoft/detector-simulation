#!/usr/bin/env python
"""This script converts a FLUKA binary file to an SLCIO file with LCIO::MCParticle instances"""

import os
import argparse
import numpy as np


parser = argparse.ArgumentParser(description='Convert FLUKA binary file to SLCIO file with MCParticles')
parser.add_argument('file_in', metavar='FILE_IN', help='Path to the input FLUKA file')
parser.add_argument('file_out', metavar='FILE_OUT.slcio', help='Path to the output SLCIO file')
parser.add_argument('-b', '--bx_time', metavar='TIME',  help='Time of the bunch crossing [ns]', type=float, default=117.78143152396451)
parser.add_argument('-n', '--normalization', metavar='N',  help='Normalization of the generated sample', type=float, default=1.0)
parser.add_argument('-l', '--lines_event', metavar='L',  help='Number of lines to put in a single LCIO event (default: 1000)', type=int, default=1000)
parser.add_argument('-m', '--max_lines', metavar='M',  help='Maximum number of lines to process', type=int, default=None)
parser.add_argument('-o', '--overwrite',  help='Overwrite existing output file', action='store_true', default=False)
parser.add_argument('--pdgs', metavar='ID',  help='PDG IDs of particles to be included', type=int, default=None, nargs='+')
parser.add_argument('--nopdgs', metavar='ID',  help='PDG IDs of particles to be excluded', type=int, default=None, nargs='+')
parser.add_argument('--ne_min', metavar='E',  help='Minimum energy of accepted neutrons [GeV]', type=float, default=None)
parser.add_argument('--t_max', metavar='T',  help='Maximum time of accepted particles [ns]', type=float, default=None)

args = parser.parse_args()

if not args.overwrite and os.path.isfile(args.file_out):
	raise FileExistsError(f'Output file already exists: {args.file_out:s}')


from math import sqrt
from pdb import set_trace as br
from array import array
from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL

import random
import math

from bib_pdgs import FLUKA_PIDS, PDG_PROPS

def bytes_from_file(filename):
	with open(filename, 'rb') as f:
		while True:
			chunk = np.fromfile(f, dtype=line_dt, count=1)
			if not len(chunk):
				return
			yield chunk

# Binary format of a single entry
line_dt=np.dtype([
	('fid',  np.int32),
	('fid_mo',  np.int32),
	('E', np.float64),
	('x', np.float64),
	('y', np.float64),
	('z', np.float64),
	('cx', np.float64),
	('cy', np.float64),
	('cz', np.float64),
	('age', np.float64),
	('age_mu', np.float64),
	('x_mu', np.float64),
	('y_mu', np.float64),
	('z_mu', np.float64),
	('x_mo', np.float64),
	('y_mo', np.float64),
	('z_mo', np.float64),
	('px_mo', np.float64),
	('py_mo', np.float64),
	('pz_mo', np.float64),
	('age_mo', np.float64)
])

######################################## Start of the processing
print(f'Converting data from file\n  {args.file_in:s}\nto SLCIO file\n  {args.file_out:s}\nwith normalization: {args.normalization:.1f}')
print(f'Splitting by {args.lines_event:d} lines/event');
if args.pdgs is not None:
	print(f'Will only use particles with PDG IDs: {args.pdgs}')

# Initialize the LCIO file writer
wrt = IOIMPL.LCFactory.getInstance().createLCWriter()
wrt.open(args.file_out, EVENT.LCIO.WRITE_NEW)

# Write a RunHeader
run = IMPL.LCRunHeaderImpl()
run.setRunNumber(0)
run.parameters().setValue('InputPath', args.file_in)
run.parameters().setValue('Normalization', args.normalization)
if args.t_max:
	run.parameters().setValue('Time_max', args.t_max)
if args.ne_min:
	run.parameters().setValue('NeutronEnergy_min', args.ne_min)
if args.pdgs:
	run.parameters().setValue('PdgIds', str(args.pdgs))
if args.nopdgs:
	run.parameters().setValue('NoPdgIds', args.nopdgs)
wrt.writeRunHeader(run)

# Bookkeeping variables
random.seed()
nEventLines = 0
nEvents = 0
col = None
evt = None

# Reading the the file one entry at a time
for iL, data in enumerate(bytes_from_file(args.file_in)):
	if nEventLines == 0:
		col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE)
		# Creating the LCIO event and collection
		evt = IMPL.LCEventImpl()

		evt.setEventNumber(nEvents)
		evt.addCollection(col, 'MCParticle')
	if args.max_lines and iL >= args.max_lines:
		break

	# Extracting relevant values from the line
	fid,e, x,y,z, cx,cy,cz, toff,toff_mo = (data[n][0] for n in [
		'fid', 'E',
		'x','y','z',
		'cx', 'cy', 'cz',
		'age', 'age_mo'
	])

	# Converting FLUKA ID to PDG ID
	try:
		pdg = FLUKA_PIDS[fid]
	except KeyError:
		print(f'WARNING: Unknown PDG ID for FLUKA ID: {pid}')
		continue

	# Calculating the absolute time of the particle [ns]
	t = toff - toff_mo - args.bx_time

	# Skipping if particle's time is greater than allowed
	if args.t_max is not None and t > args.t_max:
		continue

	# Calculating the components of the momentum vector
	mom = np.array([cx, cy, cz], dtype=np.float32)
	mom *= e

	# Skipping if it's a neutron with too low kinetic energy
	if args.ne_min is not None and abs(pdg) == 2112 and np.linalg.norm(mom) < args.ne_min:
		continue

	# Getting the charge and mass of the particle
	if pdg not in PDG_PROPS:
		print('WARNING! No properties defined for PDG ID: {0:d}'.format(pdg))
		print('         Skpping the particle...')
		continue
	charge, mass = PDG_PROPS[pdg]

	# Calculating how many random copies of the particle to create according to the weight
	nP_frac, nP = math.modf(args.normalization)
	if nP_frac > 0 and random.random() < nP_frac:
		nP += 1
	nP = int(nP)

	# Creating the particle with original parameters
	particle = IMPL.MCParticleImpl()
	particle.setPDG(pdg)
	particle.setGeneratorStatus(1)
	particle.setTime(t)
	particle.setMass(mass)
	particle.setCharge(charge)
	pos = np.array([x, y, z], dtype=np.float64)

	# Creating the particle copies with random Phi rotation
	px, py, pz = mom
	for i, iP in enumerate(range(nP)):
		p = IMPL.MCParticleImpl(particle)
		# Rotating position and momentum of the copies by a random angle in Phi
		if i > 0:
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
