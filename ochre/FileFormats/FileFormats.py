def guess_filetype(ext, magic_byte=None):
    ext_to_type = dict((a, plg.file_type) for plg in get_plugins() \
      for a in plg.abbrevs)

    ftype = ''
    ftype = ext_to_type.get(ext, '')
    #TODO: move magic byte stuff into modules
    if ftype == '' and magic_byte is not None:
        if magic_byte == 0x3E:
            ftype = 'fa'
        elif magic_byte == 0x40:
            ftype = 'fq'
    return ftype


def file_reader(ftype, *args, **kwargs):
    for pl in get_plugins():
        if ftype in pl.abbrevs:
            return pl().read(*args, **kwargs)


def file_writer(ftype, *args, **kwargs):
    for pl in get_plugins():
        if ftype in pl.abbrevs:
            return pl().write(*args, **kwargs)


def get_plugins():
    from ochre.FileFormats.FASTA import FASTA
    from ochre.FileFormats.FASTQ import FASTQ

    return [FASTA, FASTQ]

#    import types
#    import inspect
#    import ochre.FileFormats
#
#    plugins = []
#    for mdle in ochre.FileFormats.__dict__.values():
#        if type(mdle) is types.ModuleType:
#            for obj in mdle.__dict__.values():
#                if inspect.isclass(obj):
#                    if issubclass(obj, SeqFile) and obj is not SeqFile:
#                        plugins.append(obj)
#    return plugins


class SeqFile(object):
    pass
