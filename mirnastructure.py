from __future__ import print_function, absolute_import, division, unicode_literals

import logging, re


def diagramFromStructure(sequence, structure):
    """
    Convert a sequence and a vienna structure to a diagram. It assumes the 
    structure contains only one loop (the stem)

    eg for sequence
CUCAAGGUCCUGAUGAUGUGGAGCCACUCCUAUGGGGGUAGCUGUGGCUUUAACUUUGGGGGCAACUUUGAG
    and structure
(((((((((((.......(((((((((..((((....))))..))))))))).......)))..))))))))
    return

        --   gaugaug         uc    g
cucaaggu  ccu       uggagccac  cuau g
||||||||  |||       |||||||||  ||||
gaguuuca  ggg       auuucggug  gaug g
        ac   gguuuca         uc    g

    """

    loop_finder = re.compile(r'\([^()]*\)')
    loops = [loop for loop in loop_finder.finditer(structure)]
    if len(loops) != 1:
        raise ValueError("We expected only one loop but we got this: %s" % structure)

    # initialise each of the lines. We keep them in list form until the end since strings are immutable
    unmatched_upstream = []
    matched_upstream = []
    linkages = []
    matched_downstream = []
    unmatched_downstream = []
    sequence = sequence.lower()

    loop = loops[0]
    loop_start = loop.start() + 1
    loop_end = loop.end() - 2

    # draw the loop
    loop_length = loop_end - loop_start + 1

    def add_character(unmatched_up, matched_up, link, matched_down, unmatched_down):
        """
        Add a column to the drawing
        """
        unmatched_upstream.append(unmatched_up)
        matched_upstream.append(matched_up)
        linkages.append(link)
        matched_downstream.append(matched_down)
        unmatched_downstream.append(unmatched_down)


    centre_point = loop_start + loop_length // 2
    upstream_structure = structure[:centre_point][::-1] # reversed
    upstream_sequence = sequence[:centre_point][::-1]

    if (loop_length % 2) == 0:
        # standard case
        downstream_structure = structure[centre_point:]
        downstream_sequence = sequence[centre_point:]
    else:
        # odd-length. hmm, assign an extra base unmatched to the top
        downstream_structure = structure[centre_point+1:]
        downstream_sequence = sequence[centre_point+1:]

        # draw the outstanding base
        add_character(' ', sequence[centre_point], ' ', '-', ' ')

    up_index = down_index = 0 

    # draw the end of loop special-cased since we draw it a bit different
    if upstream_structure[up_index] == '.':
        add_character(' ', upstream_sequence[up_index], ' ', downstream_sequence[down_index], ' ')
        up_index += 1
        down_index += 1

    up_length = len(upstream_structure)
    down_length = len(downstream_structure)

    while up_index < up_length and down_index < down_length:
        if upstream_structure[up_index] == '(' and downstream_structure[down_index] == ')':
            # matched character
            add_character(' ', upstream_sequence[up_index], '|', downstream_sequence[down_index], ' ')
            up_index += 1
            down_index += 1
        elif upstream_structure[up_index] == '.' and downstream_structure[down_index] == '.':
            # bilateral bulge
            add_character(upstream_sequence[up_index], ' ', ' ', ' ', downstream_sequence[down_index])
            up_index += 1
            down_index += 1
        elif upstream_structure[up_index] == '.' and downstream_structure[down_index] == ')':
            # upstream bulge
            add_character(upstream_sequence[up_index], ' ', ' ', ' ', '-')
            up_index += 1
        elif upstream_structure[up_index] == '(' and downstream_structure[down_index] == '.':
            # downstream bulge
            add_character('-',' ', ' ', ' ', downstream_sequence[down_index])
            down_index += 1
        else:
            # shouldn't be anything else
            raise ValueError("Strange pairing. Upstream: %s, downstream: %s" % (upstream_structure[up_index], downstream_structure[down_index]))


    # Maybe some unmatched downstream characters left
    while down_index < down_length:
        add_character('-',' ',' ',' ', downstream_sequence[down_index])
        down_index += 1

    # and the same for upstream
    while up_index < up_length:
        add_character(upstream_sequence[up_index],' ',' ',' ', '-')
        up_index += 1

    # convert to strings but remember to reverse them before since we were drawing them backwards
    return ''.join(unmatched_upstream[::-1]), ''.join(matched_upstream[::-1]), ''.join(linkages[::-1]), ''.join(matched_downstream[::-1]), ''.join(unmatched_downstream[::-1])


if __name__ == '__main__':
    # some testing

    up_unmatched, up_matched, linkage, down_matched, down_unmatched = diagramFromStructure('CUCAAGGUCCUGAUGAUGUGGAGCCACUCCUAUGGGGGUAGCUGUGGCUUUAACUUUGGGGGCAACUUUGAG', '(((((((((((.......(((((((((..((((....))))..))))))))).......)))..))))))))')

    assert up_unmatched == '        --   gaugaug         uc    g '
    assert up_matched == 'cucaaggu  ccu       uggagccac  cuau g'
    assert linkage == '||||||||  |||       |||||||||  ||||  '
    assert down_matched == 'gaguuuca  ggg       auuucggug  gaug g'
    assert down_unmatched == '        ac   gguuuca         uc    g '


    up_unmatched, up_matched, linkage, down_matched, down_unmatched = diagramFromStructure('UUUGGAAAUUUGAUGAUCUCAAGGUCCUGAUGAUGUGGAGCCACUCCUAUGGGGGUAGCUGUGGCUUUAACUUUGGGGGCAACUUUGAGGGAAUAAAGUCUCAAAA', '....((.((((.....((((((((((((.......(((((((((..((((....))))..))))))))).......)))..))))))))).....)))).))....')

    assert up_unmatched == 'uuug  a    gauga         --   gaugaug         uc    g '
    assert up_matched == '    ga auuu     ucucaaggu  ccu       uggagccac  cuau g'
    assert linkage == '    || ||||     |||||||||  |||       |||||||||  ||||  '
    assert down_matched == '    cu ugaa     ggaguuuca  ggg       auuucggug  gaug g'
    assert down_unmatched == 'aaaa  c    auaag         ac   gguuuca         uc    g '

