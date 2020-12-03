#!/usr/bin/env python3

import sys
import re


def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    print('Reading ploidy output from %s and SMN_caller output from %s' % (file1, file2))

    positions = [70_940_000, 70_950_000]
    ploidies = [ {}, {} ]
    with open(file1) as inp:
        for line in inp:
            if line.startswith('#'):
                continue

            line = line.strip().split('\t')
            reg_start = int(line[1])
            reg_end = int(line[2])
            reg_name = line[3]
            assert reg_name == 'SMN1'
            sample = line[4]

            for i, pos in enumerate(positions):
                if reg_start <= pos < reg_end:
                    ploidies[i][sample] = line
    
    with open(file2) as inp:
        next(inp)
        for line in inp:
            smnc_line = line.strip().split('\t')
            sample = smnc_line[0].split('.')[0]

            det_line1 = ploidies[0].get(sample)
            det_line2 = ploidies[1].get(sample)
            if det_line1 is None and det_line2 is None:
                continue

            try:
                smnc_x, smnc_y, smnc_z = tuple(map(int, smnc_line[3:6]))
                det_ploidy1 = int(det_line1[7])
                det_ploidy2 = tuple(map(int, det_line2[10].split(',')))
                match = True
            except ValueError:
                match = False

            if match:
                match &= det_ploidy1 == smnc_x + smnc_y + smnc_z
                match &= det_ploidy2 == (smnc_x, smnc_y)

            partial_mism = False
            if match:
                try:
                    det_par_ploidy1 = tuple(map(int, det_line1[10].split(',')))
                    partial_mism = det_par_ploidy1 != (smnc_x, smnc_y + smnc_z)
                except ValueError:
                    pass

            if not match or partial_mism:
                print('{}Mismatch for sample {}:'.format('PARTIAL ' if partial_mism else '', sample))
                print('\t{}'.format('\t'.join(smnc_line)))
                print('\t{}'.format('\t'.join(det_line1) if det_line1 else '*'))
                print('\t{}'.format('\t'.join(det_line2) if det_line2 else '*'))
                print()

if __name__ == '__main__':
    main()
