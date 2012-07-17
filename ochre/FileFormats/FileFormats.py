ext_to_type = {'fa': 'fa',
               'fasta': 'fa',
               'fq': 'fq',
               'fastq': 'fq',
               'gb': 'gb',
               'genbank': 'gb',
               'embl': 'em',
               'sff': 'sff'}
file_types = set(ext_to_type.values())


def guess_filetype(ext, magic_byte=None):
    ftype = ''
    ftype = ext_to_type.get(ext, '')
    if ftype == '' and magic_byte is not None:
        if magic_byte[0] == 0x3E:
            ftype = 'fa'
        elif magic_byte[0] == 0x40:
            ftype = 'fq'
    return ftype
