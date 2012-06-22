file_types = ['fa','fq','gb','em']

def guess_filetype(ext, magic_byte=None):
    ftype = ''
    ext_to_type = {'fa':'fa',
                   'fasta':'fa',
                   'fq':'fq',
                   'fastq':'fq',
                   'gb':'gb',
                   'genbank':'gb',
                   'embl':'em'}
    ftype = ext_to_type.get(ext,'')
    if ftype == '' and magic_byte is not None:
        if magic_byte[0] == 0x3E:
            ftype = 'fa'
        elif magic_byte[0] == 0x40:
            ftype = 'fq'
    return ftype
