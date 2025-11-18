from typing import Tuple
import sqlite3
from pydantic import BaseModel
from typing import List
import logging
import pickle
from .utils import cmd_out
from tqdm import tqdm
import filelock
import numpy as np

class Link(BaseModel):
    source: str
    target: str
    distance: float
    positions: List[int]
    missing: int

def extract_variant_set(vcf_file: str) -> Tuple[set,set]:
    ref_diffs = set()
    missing = set()
    for l in cmd_out(f"bcftools view {vcf_file} | bcftools query -f '%POS[\t%GT]\n'"):
        if l[0]=="#": continue
        row = l.strip().split()
        pos = int(row[0])
        gt = row[1]
        if gt==".": 
            missing.add(pos)
            continue
        elif gt=="1":
            ref_diffs.add(int(pos))
        else:
            raise Exception("Unknown GT: %s" % gt)

    return ref_diffs, missing



class SnpDistDB:
    """
    Class for storing and searching for SNP differences between samples
    
    Arguments
    ---------
    filename : str
        Filename of the sqlite database to use
    
    Attributes
    ----------
    filename : str
        Filename of the sqlite database to use
    conn : sqlite3.Connection
        Connection to the sqlite database
    c : sqlite3.Cursor
        Cursor to the sqlite database
        
    Methods
    -------
    store(result,vcf_file)
        Store the SNP differences from a sample in the database
    search(result,vcf_file,cutoff=20)
        Search for samples with similar SNP differences in the database
    """
    def __init__(self, filename: str) -> None:
        with filelock.SoftFileLock(f"{filename}.lock"):
            self.filename = filename
            self.conn = sqlite3.connect(filename)
            self.c = self.conn.cursor()
            self.c.execute('''CREATE TABLE IF NOT EXISTS samples (sample text, taxa text, diffs binary, missing binary)''')
            self.c.execute('''CREATE TABLE IF NOT EXISTS links (source text, target text, snps binary, dist integer, missing integer)''')
    
    def store(self, sample_name:str, vcf_file: str, taxa: str, cutoff: int = 10) -> List[Link]:
        """
        Store the SNP differences from a sample in the database.
        
        Parameters
        ----------
        sample_name : str
            Name of the sample
        vcf_file : str
            VCF file containing the SNP differences
        taxa : str
            Taxonomic classification of the sample
        cutoff : int, optional
            Maximum SNP distance to consider a sample as similar, by default 10

        Returns
        -------
        List[Link]
            List of Link objects representing similar samples
        """
        new_diffs, new_missing = extract_variant_set(vcf_file)
        res = self.c.execute("SELECT sample FROM samples WHERE sample=?",(sample_name,)).fetchone()
        if not res:
            self.c.execute("INSERT INTO samples VALUES (?,?,?,?)",(sample_name, taxa, pickle.dumps(new_diffs), pickle.dumps(new_missing)))
        

        logging.info("Searching for close samples in %s" % self.filename)

        self.c.execute("SELECT sample, diffs, missing FROM samples WHERE taxa=?", (taxa,))
        
        sample_dists = []
        for s,d,m in tqdm(self.c.fetchall(),desc="Searching for close samples"):
            dist = new_diffs.symmetric_difference(pickle.loads(d))
            dist -= new_missing
            dist -= pickle.loads(m) 
            logging.debug(f"Sample {sample_name} vs {s}: {len(dist)} SNPs different")
            if (ld:=len(dist))<cutoff:
                sample_dists.append(
                    Link(
                        source = sample_name,
                        target = s,
                        distance = ld,
                        positions = list(dist),
                        missing = len(new_missing.union(pickle.loads(m)))
                    )
                )

        self.c.executemany("INSERT INTO links VALUES (?,?,?,?,?)", [
            (l.source, l.target, pickle.dumps(l.positions), l.distance, l.missing)
            for l in sample_dists if l.source != l.target
        ])
        self.conn.commit()
        return sample_dists
    
    def extract_matrix(self) -> Tuple[List[str], List[List[int]]]:
        """
        Extract a distance matrix from the database.

        Returns
        -------
        Tuple[List[str], List[List[int]]]
            A tuple containing a list of sample names and a 2D list representing the distance matrix.
        """
        samples = []
        sample_idx = {}
        self.c.execute("SELECT sample FROM samples")
        for i, (s,) in enumerate(self.c.fetchall()):
            samples.append(s)
            sample_idx[s] = i
        n = len(samples)
        matrix = [[np.nan for _ in range(n)] for _ in range(n)]
        self.c.execute("SELECT source, target, dist FROM links")
        for s,t,d in self.c.fetchall():
            i = sample_idx[s]
            j = sample_idx[t]
            matrix[i][j] = d
            matrix[j][i] = d
        for i in range(n):
            matrix[i][i] = 0
        return samples, matrix
    def _sample_present(self, sample_name: str) -> bool:
        res = self.c.execute("SELECT sample FROM samples WHERE sample=?",(sample_name,)).fetchone()
        return res is not None
    def inspect_link(self, source: str, target: str) -> Link:
        """
        Inspect a specific link between two samples.

        Parameters
        ----------
        source : str
            Name of the source sample
        target : str
            Name of the target sample

        Returns
        -------
        Link
            A Link object representing the SNP differences between the two samples.
        """
        if not self._sample_present(source):
            raise Exception(f"Sample {source} not found in database")
        if not self._sample_present(target):
            raise Exception(f"Sample {target} not found in database")
        res1 = self.c.execute("SELECT snps, dist, missing FROM links WHERE source=? AND target=?", (source, target)).fetchone()
        res2 = self.c.execute("SELECT snps, dist, missing FROM links WHERE source=? AND target=?", (target, source)).fetchone()
        if not res1 and not res2:
            raise Exception(f"No link found between {source} and {target}")
        snps, dist, missing = res1 if res1 else res2
        return Link(
            source=source,
            target=target,
            distance=dist,
            positions=pickle.loads(snps),
            missing=missing
        )