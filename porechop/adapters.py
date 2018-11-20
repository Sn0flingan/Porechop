"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the class and sequences for known adapters used in Oxford Nanopore library
preparation kits.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""


class Adapter(object):

    def __init__(self, name, start_sequence=None, end_sequence=None, both_ends_sequence=None):
        self.name = name
        self.start_sequence = start_sequence if start_sequence else []
        self.end_sequence = end_sequence if end_sequence else []
        if both_ends_sequence:
            self.start_sequence = both_ends_sequence
            self.end_sequence = both_ends_sequence
        self.best_start_score, self.best_end_score = 0.0, 0.0

    def best_start_or_end_score(self):
        return max(self.best_start_score, self.best_end_score)

    def is_barcode(self):
        if 'Barcode' in self.name:
            return True
        else:
            return False

    def barcode_direction(self):
        if '_rev' in self.start_sequence[0]:
            return 'reverse'
        else:
            return 'forward'

    def get_barcode_name(self):
        """
        Gets the barcode name for the output files. We want a concise name, so it looks at all
        options and chooses the shortest.
        """
        '''
        possible_names = [self.name, self.start_sequence[0]]
        if self.end_sequence:
            possible_names.append(self.end_sequence[0])
        barcode_name = sorted(possible_names, key=lambda x: len(x))[0]
        return barcode_name.replace(' ', '_')
        '''
        return self.name
            
    def get_name(self):
        return self.name

    def get_start_seq(self):
        return self.start_sequence

    def get_end_seq(self):
        return self.end_sequence

    def toString(self):
        return "--------------\n" + self.name + "\n" + self.start_sequence[1]  + "\n" + self.end_sequence[1] + "\n"


##Motor protein binding part of top adapter sequence was removed
kit_adapters = {"SQK-LSK109": Adapter('SQK-LSK109',
                    start_sequence=('SQK-NSK007_Y_Top_trunk', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
                    end_sequence=('SQK-NSK007_Y_Bottom', 'GCAATACGTAACTGAACGAAGT')),
                "EXP-PBC096-SQK-LSK108": Adapter('EXP-PBC096-SQK-LSK108',
                    start_sequence=('SQK-NSK007_Y_Top_trunk-EXP-PBC096', 'AATGTACTTCGTTCAGTTACGTATTGCTGGTGCTG'),
                    end_sequence=('SQK-NSK007_Y_Bottom-EXP-PBC096', 'GCAATACGTAACTGAACGAAGT')),
                "SQK-LSK108": Adapter('SQK-LSK108',
                    start_sequence=('SQK-NSK007_Y_Top_trunk', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
                    end_sequence=('SQK-NSK007_Y_Bottom', 'GCAATACGTAACTGAACGAAGT')),
                "SQK-LSK308": Adapter('SQK-LSK308',
                    start_sequence=('SQK-LSK308_1D2_Top', 'GTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT'),
                    end_sequence=('SQK-LSK308_1D2_Bottom','GGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT')),
                "SQK-RAD004": Adapter('SQK-RAD004',
                      start_sequence=('SQK-RAD004_Y_Top_trunk', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
                                      end_sequence=('SQK-RAD004_Y_Bottom','GCAATACGTAACTGAACGAAGT'))}

def make_full_native_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (reverse)'][0]
    start_barcode_seq = barcode.start_sequence[1]
    end_barcode_seq = barcode.end_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA' + start_barcode_seq + 'CAGCACCT'
    end_full_seq = 'AGGTGCTG' + end_barcode_seq + 'TTAACCTTAGCAATACGTAACTGAACGAAGT'

    return Adapter('Native barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('NB' + '%02d' % barcode_num + '_start', start_full_seq),
                   end_sequence=('NB' + '%02d' % barcode_num + '_end', end_full_seq))


def make_full_rapid_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCT' + start_barcode_seq + 'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))
