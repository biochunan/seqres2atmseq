# main.py
import os
import json
import shutil
import argparse
from pathlib import Path
from loguru import logger
from pprint import pprint
from typing import Dict, List, Union, Tuple
# Biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqIO import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline
import tempfile

# ==================== Configuration ====================
PathLike = Union[str, Path]


# ==================== Function ====================
# align seq using ClustalOmega
def run_align_clustalomega(clustal_omega_executable: str,
                           seq1: str = None, seq2: str = None,
                           seqs: List[str] = None) -> List[SeqRecord]:
    """

    Args:
        seq1: sequence of a chain e.g. seqres sequence 
        seq2: sequence of a chain e.g. atmseq sequence 
        or you can provide a list of strings using seqs 
        seqs: e.g. ["seq1", "seq2", ...]
        clustal_omega_executable: (str) path to clustal omega executable
            e.g. "/usr/local/bin/clustal-omega"
    Returns:
        aln_seq_records: (List)
    """
    # assert input
    if seqs is None and (seq1 is None or seq2 is None):
        raise NotImplemented("Provide either List of seqs as `seqs` OR a pair of seqs as `seq1` and `seq2`.")

    # generate seq_recs
    seq_rec = [None]
    if seqs:
        seq_rec = [SeqRecord(id=f"seq{i + 1}", seq=Seq(seqs[i]), description="")
                   for i in range(len(seqs))]
    elif seq1 is not None and seq2 is not None:
        seq_rec = [SeqRecord(id=f"seq{1}", seq=Seq(seq1), description=""),
                   SeqRecord(id=f"seq{2}", seq=Seq(seq2), description="")]

    with tempfile.TemporaryDirectory() as tmpdir:
        # executable
        cmd = clustal_omega_executable

        # create input seq fasta file and output file for clustal-omega
        in_file = os.path.join(tmpdir, "seq.fasta")
        out_file = os.path.join(tmpdir, f"aln.fasta")
        with open(in_file, "w") as f:
            SeqIO.write(seq_rec, f, "fasta")
        # create Clustal-Omega commands
        clustalomega_cline = ClustalOmegaCommandline(cmd=cmd, infile=in_file, outfile=out_file, verbose=True, auto=True)

        # run Clustal-Omega
        stdout, stderr = clustalomega_cline()

        # read aln
        aln_seq_records = []
        with open(out_file, "r") as f:
            for record in AlignIO.read(f, "fasta"):
                aln_seq_records.append(record)

        return aln_seq_records

def gen_seqres_to_atmseq_mask(seqres: str, 
                              atmseq: str,
                              clustal_omega_executable: Union[Path, str]
    ) -> Tuple[List[int], List[SeqRecord]]:
    aln = run_align_clustalomega(clustal_omega_executable=clustal_omega_executable,
                                 seq1=seqres,
                                 seq2=atmseq)
    # assert no dash in seqres
    try:
        assert '-' not in str(aln[0].seq)
    except AssertionError:
        print(f'seqres: {seqres}')
        print(f'atmseq: {atmseq}')
        print(f'aln[0].seq: {aln[0].seq}')
        print(f'aln[1].seq: {aln[1].seq}')
        raise AssertionError('seqres contains dash, please check your input seqres and atmseq')
    
    aln1 = str(aln[1].seq)  # => this is the atmseq in aln, may contain "-"
    seqres_to_atmseq_mask=[1 if i != '-' else 0 for i in aln1]  # 1 => in atmseq; 0 => not in atmseq
    
    return seqres_to_atmseq_mask, aln

def parse_fasta_file(fasta_fp: PathLike) -> Dict[str, str]:
    """
    Parse fasta file. Add id (sequence header) and sequence to 
    the dictionary metadata field.

    Args:
        fasta_fp (str): fasta file path

    Returns:
        Dict[str, str]: dictionary of sequence headers and sequences
    """
    seq_dict = {'headers': [], 'seqs': []}
    with open(fasta_fp, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq_dict['headers'].append(record.id)
            seq_dict['seqs'].append(str(record.seq))
    
    return seq_dict

def validate_sequence_headers(seq_dict: Dict[str, str]) -> None:
    """
    Each chain sequence header should follow the format
    <structure name>|<chain id>|<seqres/atmseq>)
    e.g. 1a2y_0P|H|seqres
    e.g. 1a2y_0P|L|atmseq
    And each chain should have both seqres and atmseq.

    Args:
        seq_dict (Dict[str, str]): a dictionary returned by parse_fasta_file. 
            It contains sequence headers and sequences.
    """
    # get all seqres/atmseq pairs 
    seqres_atmseq_pairs = {}
    for header in seq_dict['headers']:
        seq_name, seq_type = [x.strip() for x in header.rsplit('|', 1)]
        if seq_name not in seqres_atmseq_pairs:
            seqres_atmseq_pairs[seq_name] = {'seqres': False, 'atmseq': False}
        seqres_atmseq_pairs[seq_name][seq_type] = True
        
    # check if each chain has both seqres and atmseq
    for seq_name, seqres_atmseq in seqres_atmseq_pairs.items():
        try:
            assert seqres_atmseq['seqres'] and seqres_atmseq['atmseq']
        except AssertionError:
            raise AssertionError(
                f"Each chain should have both seqres and atmseq. {seq_name} only has {seqres_atmseq}"
            )
    
    return True

def process_seq_dict(seq_dict: Dict[str, str]) -> Dict[str, str]:
    """
    Reorganize seq_dict to a dictionary of 
    {
        f"{structure name}|{chain id}": {
            'seqres': str,
            'atmseq': str,
        }
    }

    Args:
        seq_dict (Dict[str, str]): a dictionary returned by parse_fasta_file. 
            It contains sequence headers and sequences.

    Returns:
        Dict[str, str]: See above
    """
    # get all seqres/atmseq pairs 
    seqres_atmseq_pairs = {}
    
    for h, s in zip(seq_dict['headers'], seq_dict['seqs']):
        seq_name, seq_type = [x.strip() for x in h.rsplit('|', 1)]
        if seq_name not in seqres_atmseq_pairs:
            seqres_atmseq_pairs[seq_name] = {seq_type: s}
        else:
            seqres_atmseq_pairs[seq_name][seq_type] = s
    
    return seqres_atmseq_pairs

def cli():
    parser = argparse.ArgumentParser(description="Process sequence file.")
    
    # Define arguments
    parser.add_argument('-i', '--fasta', type=Path, help='Input FASTA file path', default=None)
    parser.add_argument('-s', '--seqres', type=str, help='SEQRES sequence', default=None)
    parser.add_argument('-a', '--atmseq', type=str, help='ATMSEQ sequence', default=None)
    parser.add_argument('-n', '--seq_name', type=str, help='Sequence name', default='seq')
    parser.add_argument('-c', '--clustal_omega_executable', 
                        type=str, default=None,
                        help='Path to clustal omega executable')
    # add output dir or file path 
    parser.add_argument('-o', '--output', type=Path, default=Path('./mask.json'),
                        help='Output directory or file path')
    # verbose 
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='Verbose mode')
    
    args = parser.parse_args()
    
    # ---------------------------------------------------------------
    # Check for correct combination of arguments, either 
    # 1. --fasta or 
    # 2. both --seqres and --atmseq must be provided.
    # ---------------------------------------------------------------
    if args.fasta is None and (args.seqres is None or args.atmseq is None):
        parser.error("Either --fasta or both --seqres and --atmseq must be provided.")
    elif args.fasta is not None and (args.seqres is not None or args.atmseq is not None):
        parser.error("Provide either --fasta or both --seqres and --atmseq, not both.")

    # use default clustal omega executable if not provided
    if args.clustal_omega_executable is None:
        args.clustal_omega_executable = shutil.which('clustalo')
    # raise error if clustal omega executable not found
    if args.clustal_omega_executable is None:
        raise FileNotFoundError('clustal omega executable not found')
    
    def _set_clustalo_path_by_os_type() -> None:
        """
        Set clustalo path based on OS type
        if OS is Mac, use clustaloMac -> clustal-omega-1.2.3-macosx
        if OS is Linux, use clustalo -> clustalo-1.2.4-Ubuntu-x86_64
        
        Raises:
            NotImplementedError: OS type not supported
        Returns:
            None
        """
        # determine OS type
        import platform
        os_type = platform.system()
        if os_type == 'Darwin':
            args.clustal_omega_executable = Path(__file__).parent/'assets'/'clustaloMac'
        elif os_type == 'Linux':
            args.clustal_omega_executable = Path(__file__).parent/'assets'/'clustalo'
        else:
            raise NotImplementedError(f'OS type {os_type} not supported.')
    
    # set clustalo path by OS type if not provided
    if args.clustal_omega_executable is None:
        _set_clustalo_path_by_os_type()
    
    return args

def main():
    args = cli()
    
    # process inputs 
    if args.fasta: 
        try:
            assert args.fasta.exists()
        except AssertionError:
            raise FileNotFoundError(f'{args.fasta} not found')
        seq_dict = parse_fasta_file(fasta_fp=args.fasta)
    else: 
        seq_dict = {'headers': [f'{args.seq_name}|seqres', f'{args.seq_name}|atmseq'], 
                    'seqs': [args.seqres, args.atmseq]}
    
    # validate sequence headers
    validate_sequence_headers(seq_dict=seq_dict)
    
    # reorganize seq_dict
    seqres_atmseq_pairs = process_seq_dict(seq_dict=seq_dict)
    
    # iterate over seqres/atmseq pairs and 
    # run clustal omega and generate seqres to atmseq mask
    mask_dict = {}
    for k in seqres_atmseq_pairs:
        mask, aln = gen_seqres_to_atmseq_mask(seqres=seqres_atmseq_pairs[k]['seqres'],
                                              atmseq=seqres_atmseq_pairs[k]['atmseq'],
                                              clustal_omega_executable=args.clustal_omega_executable)
        mask_dict[k] = {'seqres': str(aln[0].seq), 
                        'mask': mask, 
                        'atmseq': str(aln[1].seq)}
        if args.verbose:
            info = '\n'
            info += f"{k} seqres: {str(aln[0].seq)}\n"
            info += f"{k} atmseq: {str(aln[1].seq)}\n"
            info += f"{k} mask  : {''.join(map(str, mask))}\n"
            logger.debug(info)
    
    # save mask as json 
    if args.output.is_dir():
        args.output.mkdir(parents=True, exist_ok=True)
        args.output = args.output/'mask.json'
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        # convert list to str
        for k in mask_dict:
            m = mask_dict[k]['mask']
            mask_dict[k]['mask'] = ''.join(map(str, m))
        json.dump(mask_dict, f, indent=4)


# ==================== Main ====================
if __name__ == "__main__":
    main()
