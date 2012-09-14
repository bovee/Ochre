""" Utility commands for simplifying modules that deal with external
    programs."""
import os.path as op
import pipes
import subprocess
import ochre


def get_cfg(header, key):
    try:
        import configparser
    except ImportError:
        import ConfigParser as configparser
    cfg = configparser.RawConfigParser()
    cfg.read([op.join(ochre._path, 'ochre.cfg')])
    return cfg.get(header, key)


def app(app_name, subprog_name):
    """ Resolves the directory for app_name and
    appends the program named subprog_name to it."""
    apppath = get_cfg('Locales', app_name).replace('$OCHRE', ochre._path)
    return op.join(apppath, subprog_name)


def temp(*ident):
    """ Appends ident, if given, to the temporary directory. """
    return op.join(get_cfg('Locales', 'TEMP'), op.join('', *ident))


def run(cmds, lsf_opts=None):
    """ Given a list of commands (also as lists & suitable for
    subprocess), execute them either through subprocess or through
    a batch manager such as LSF. """
    if lsf_opts is None:
        for c in cmds:
            p = subprocess.Popen(c)
            p.wait()
            if p.returncode != 0:
                raise OSError('Line "' + ' '.join(c) + \
                  '"\n failed with error ' + str(p.returncode))
    else:
        lsf_opts = {'-u': 'bovee@fas.harvard.edu',
                '-J': 'process_illumina',
                '-q': 'bigmem',
                '-o': '/n/pearsonfs1/SequenceData/MahoneyLake7M/output.txt',
                '-n': '8',
                '-R': '"span[ptile=8];rusage[mem=50000]"'}

        cmd_str = '#!/bin/sh\n'
        for opt in lsf_opts:
            cmd_str += '#BSUB ' + opt + ' ' + lsf_opts[opt] + '\n'
        cmd_str += '\n'
        cmd_str += '\n'.join(' '.join(pipes.quote(i) for i in c) for c in cmds)

        p = subprocess.Popen(['bsub', '-K'], stdin=subprocess.PIPE)
        p.communicate(input=bytes(cmd_str, encoding='UTF-8'))
        p.wait()
