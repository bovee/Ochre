def run_cmds(cmds, lsf_opts=None):
    import pipes
    import subprocess

    if lsf_opts is None:
        for c in cmds:
            p =  subprocess.Popen(c)
            p.wait()
    else:
        lsf_opts = {'-u':'bovee@fas.harvard.edu', '-J':'process_illumina', '-q':'bigmem', '-o':'/n/pearsonfs1/SequenceData/MahoneyLake7M/output.txt', '-n':'8', '-R':'"span[ptile=8];rusage[mem=50000]"'}

        cmd_str = '#!/bin/sh\n'
        for opt in lsf_opts:
            cmd_str += '#BSUB '+opt+' '+lsf_opts[opt]+'\n'
        cmd_str += '\n'
        cmd_str += '\n'.join(' '.join(pipes.quote(i) for i in c) for c in cmds)

        p = subprocess.Popen(['bsub','-K'], stdin=subprocess.PIPE)
        p.communicate(input=bytes(cmd_str,encoding='UTF-8'))
        p.wait()
