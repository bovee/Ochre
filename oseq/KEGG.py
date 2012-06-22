
def get_genes(ko_num='01676'):
    from SOAPpy import WSDL

    with open('K'+fo_num+'.fna') as outfile:
        serv = WSDL.Proxy('http://soap.genome.jp/KEGG.wsdl')
        genes = serv.get_genes_by_ko('ko:K' + ko_num) 
        for gene in genes:
            fasta = serv.bget('-f -n n ' + gene.entry_id)
            outfile.write(fasta)
