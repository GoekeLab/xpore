import sys

from .dataprep import dataprep
from .diffmod import diffmod
from .postprocessing import postprocessing

def parse_options(argv):

    """Parses options from the command line """

    from argparse import ArgumentParser
    from xpore import __version__

    parser = ArgumentParser(prog='xpore')
    subparsers = parser.add_subparsers(help='Running modes', metavar='{dataprep, diffmod, postprocessing}')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))

    ### RUN MODE "DATAPREP"
    parser_dataprep = subparsers.add_parser('dataprep', help='run mode for preprocessing nanopolish eventalign.txt before differential modification analysis')
    optional_dataprep = parser_dataprep._action_groups.pop()
    required_dataprep = parser_dataprep.add_argument_group('required arguments')
    # Required arguments
    required_dataprep.add_argument('--eventalign', dest='eventalign', help='eventalign filepath, the output from nanopolish.',required=True)
    ##required.add_argument('--summary', dest='summary', help='eventalign summary filepath, the output from nanopolish.',required=True)
    required_dataprep.add_argument('--out_dir', dest='out_dir', help='output directory.',required=True)
    optional_dataprep.add_argument('--gtf_path_or_url', dest='gtf_path_or_url', help='GTF file path or url.',type=str)
    optional_dataprep.add_argument('--transcript_fasta_paths_or_urls', dest='transcript_fasta_paths_or_urls', help='transcript FASTA paths or urls.',type=str)
    # Optional arguments
    optional_dataprep.add_argument('--skip_eventalign_indexing', dest='skip_eventalign_indexing', help='skip indexing the eventalign nanopolish output.',default=False,action='store_true')
    # parser.add_argument('--features', dest='features', help='Signal features to extract.',type=list,default=['norm_mean'])
    optional_dataprep.add_argument('--genome', dest='genome', help='to run on Genomic coordinates. Without this argument, the program will run on transcriptomic coordinates',default=False,action='store_true') 
    optional_dataprep.add_argument('--n_processes', dest='n_processes', help='number of processes to run.',type=int, default=1)
    optional_dataprep.add_argument('--chunk_size', dest='chunk_size', help='number of lines from nanopolish eventalign.txt for processing.',type=int, default=1000000)
    optional_dataprep.add_argument('--readcount_min', dest='readcount_min', help='minimum read counts per gene.',type=int, default=1)
    optional_dataprep.add_argument('--readcount_max', dest='readcount_max', help='maximum read counts per gene.',type=int, default=1000)
    optional_dataprep.add_argument('--resume', dest='resume', help='with this argument, the program will resume from the previous run.',default=False,action='store_true') #todo
    parser_dataprep._action_groups.append(optional_dataprep)
    parser_dataprep.set_defaults(func=dataprep)

    ### RUN MODE "DIFFMOD"
    parser_diffmod = subparsers.add_parser('diffmod', help='run mode for performing differential modification analysis')
    optional_diffmod = parser_diffmod._action_groups.pop()
    required_diffmod = parser_diffmod.add_argument_group('required arguments')
    # Required arguments
    required_diffmod.add_argument('--config', dest='config', help='YAML configuraion filepath.',required=True)
    # Optional arguments
    optional_diffmod.add_argument('--n_processes', dest='n_processes', help='number of processes to run.',type=int,default=1)
    optional_diffmod.add_argument('--save_models', dest='save_models', help='with this argument, the program will save the model parameters for each id.',default=False,action='store_true') # todo
    optional_diffmod.add_argument('--resume', dest='resume', help='with this argument, the program will resume from the previous run.',default=False,action='store_true') 
    optional_diffmod.add_argument('--ids', dest='ids', help='gene or transcript ids to model.',default=[],nargs='*')
    parser_diffmod._action_groups.append(optional_diffmod)
    parser_diffmod.set_defaults(func=diffmod)

    ### RUN MODE "POSTPROCESSING"
    parser_postprocessing = subparsers.add_parser('postprocessing', help='run mode for post-processing diffmod.table, the result table from differential modification analysis.')
    required_postprocessing = parser_postprocessing.add_argument_group('required arguments')
    # Required arguments
    required_postprocessing.add_argument('--diffmod_dir', dest='diffmod_dir', help='diffmod directory path, the output from xpore-diffmod.',required=True)
    parser_postprocessing.set_defaults(func=postprocessing)

    return parser.parse_args(argv[1:])

def main(argv=sys.argv):

    ### get command line options
    options = parse_options(argv)
    options.func(options)

if __name__ == "__main__":
    main(sys.argv)
